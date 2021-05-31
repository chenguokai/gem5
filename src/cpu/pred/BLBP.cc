//
// Created by xim on 5/16/21.
//
/*
 * Design details:
 * PC bits used for prediction: bit 11 to 2,
 *
*/


#include "BLBP.hh"

#include "base/intmath.hh"
#include "debug/Indirect.hh"

BLBP::BLBP(const ITTAGEParams &params)
           : IndirectPredictor(params),
             numPredictors(params.numPredictors),
             ghrNumBits(params.indirectGHRBits),
             numTageBits(params.indirectTageBits),
             pathLength(params.indirectPathLength),
             ghrMask((1 << params.indirectGHRBits) - 1)
{
    threadInfo.resize(params.numThreads);

    WeightTable.resize(params.numThreads);
    for (unsigned i = 0; i < params.numThreads; ++i) {
        WeightTable[i].resize(numPredictors);
        for (unsigned j = 0; j < numPredictors; ++j) {
            WeightTable[i][j].resize(numWeightEntry);
        }
    }
    IBTB.resize(params.numThreads);
    for (unsigned i = 0; i < params.numThreads; ++i) {
        IBTB[i].resize(numBTBEntry);
        // no need to resize priority queue
    }

    intervalTable.push_back(make_pair(0, 13));
    intervalTable.push_back(make_pair(1, 33));
    intervalTable.push_back(make_pair(23, 49));
    intervalTable.push_back(make_pair(44, 85));
    intervalTable.push_back(make_pair(77, 149));
    intervalTable.push_back(make_pair(159, 270));
    intervalTable.push_back(make_pair(252, 630));
}

void BLBP::genIndirectInfo(ThreadID tid, void* &indirect_history)  {
    // record the GHR as it was before this prediction
    // It will be used to recovery the history in case this prediction is wrong or belongs to bad path
    indirect_history = new ghrEntry(threadInfo[tid].ghr);
}

void BLBP::shiftGhr(ghrEntry &ghr) {
    // shift GHR by 1 bit
    for (unsigned i = 0; i < 10; ++i) {
        ghr.ghr[i] <<= 1;
        if (i  > 0) {
            ghr[i] |= (ghr[i - 1] >> 63); // get carry out bit
        }
    }
}

void BLBP::shiftGhr_right(ghrEntry &ghr) {
    // shift ghr 1 bit right
    uint64_t higher = 0;
    for (unsigned i = 9; i >= 0; --i) {
        ghr.ghr[i] |= (higher << 63);
        higher = ghr.ghr[i] & 1;
        ghr.ghr[i] >>= 1;
    }
}

void BLBP::updateDirectionInfo(ThreadID tid, bool actually_taken) {
    // we need to move all bits one bit left
    // storage design:
    shiftGhr(threadInfo[tid].ghr);
    threadInfo[tid].ghr[0] |= actually_taken;
    // no need to mask
}

void BLBP::changeDirectionPrediction(ThreadID tid, void * indirect_history, bool actually_taken) {
    ghrEntry *previousGhr = static_cast<ghrEntry *>(indirect_history);
    threadInfo[tid].ghr = *previousGhr;
    shiftGhr(threadInfo[tid].ghr);
    threadInfo[tid].ghr[0] |= actually_taken;
}

bool BLBP::lookup_helper(Addr br_addr, TheISA::PCState& target, ThreadID tid) {
    // we need to lookup through

    unsigned btb_entry = addr_fold_btb(br_addr);
    std::vector<std::pair<int, TheISA::PCState> >ans; // first: the total sum second: the target address


    // get the weight table
    std::vector<int> weight_sum;
    weight_sum.resize(numWeightBits);
    // get weight from the first table
    unsigned entry = local_to_entry(br_addr);
    for (int i = 0; i < numWeightBits; i++) {
        weight_sum[i] = WeightTable[tid][0][entry].weight[i];
    }
    // get weight from the following tables
    for (int i = 1; i < numPredictors; ++i) {
        entry = history_to_entry(tid, i, threadInfo[tid].ghr);
        for (int j = 0; j < numWeightBits; ++j) {
            weight_sum[j] += WeightTable[tid][i][entry].weight[j];
        }
    }

    // get the BTB entry

    BTBEntry it;
    uint64_t available = 0xffffffffffffffffUL;
    while (!tmp_quque.empty()) tmp_quque.pop(); // we expect an empty tmp queue
    while (!IBTB[tid][btb_entry].empty()) {
        it = IBTB[tid][btb_entry].top();
        IBTB[tid][btb_entry].pop();
        available &= IBTB[tid][btb_entry].targetAddr;
        tmp_quque.push(it);
    }

    available >>= 2; // skip the lowest two bits
    // now we have common bits kept 1 and non-common ones 0

    int max_val = -1, current_val;
    TheISA::PCState max_target, current_target;
    Addr max_target_addr = 0;
    while (!tmp_quque.empty()) {
        current_val = 0;
        // check every target and get the maximum value
        current_target = tmp_quque.front().target;
        Addr target = current_target.instAddr();
        tmp_quque.pop();
        target >>= 2; // skip the lowest two bits
        for (int i = 0; i < numWeightBits; ++i) {
            if ((available & (1UL << i)) == 0) {
                // non-common bit, calculate the bits
                // iterate through all given entries
                current_val += (target & 1) * weight_sum[i];
            }
            target >>= 1;
        }
        if (current_val > max_val) {
            max_val = current_val;
            max_target_addr = current_target;
            max_target = current_target;
        }
    }
    if (!max_target_addr) {
        // we have found the target
        target = max_target;
        return true;
    }
    return false;
}

bool BLBP::lookup(Addr br_addr, TheISA::PCState& target, ThreadID tid) {
    return lookup_helper(br_addr, target, tid);
}

void BLBP::recordIndirect(Addr br_addr, Addr tgt_addr, InstSeqNum seq_num, ThreadID tid) {
    DPRINTF(Indirect, "Recording %x seq:%d\n", br_addr, seq_num);
    HistoryEntry entry(br_addr, tgt_addr, seq_num);
    threadInfo[tid].pathHist.push_back(entry);
}

void BLBP::commit(InstSeqNum seq_num, ThreadID tid, void * indirect_history) {
    DPRINTF(Indirect, "Committing seq:%d\n", seq_num);
    ThreadInfo &t_info = threadInfo[tid];

    // we do not need to recover the GHR, so delete the information
    unsigned * previousGhr = static_cast<unsigned *>(indirect_history);
    delete previousGhr;

    if (t_info.pathHist.empty()) return;

    if (t_info.headHistEntry < t_info.pathHist.size() &&
        t_info.pathHist[t_info.headHistEntry].seqNum <= seq_num) {
        if (t_info.headHistEntry >= pathLength) {
            t_info.pathHist.pop_front();
        } else {
            ++t_info.headHistEntry;
        }
    }
}

void BLBP::squash(InstSeqNum seq_num, ThreadID tid) {
    DPRINTF(Indirect, "Squashing seq:%d\n", seq_num);
    ThreadInfo &t_info = threadInfo[tid];
    auto squash_itr = t_info.pathHist.begin();
    int valid_count = 0;
    while (squash_itr != t_info.pathHist.end()) {
        if (squash_itr->seqNum > seq_num) {
            break;
        }
        ++squash_itr;
        ++valid_count;
    }
    if (squash_itr != t_info.pathHist.end()) {
        DPRINTF(Indirect, "Squashing series starting with sn:%d\n",
                squash_itr->seqNum);
    }
    int queue_size = t_info.pathHist.size();
    for (int i = 0; i < queue_size - valid_count; ++i) {
        t_info.ghr >>=1;
    }
    t_info.pathHist.erase(squash_itr, t_info.pathHist.end());

    // todo: maintain local history here
}

void
BLBP::deleteIndirectInfo(ThreadID tid, void * indirect_history)
{
    ghrEntry * previousGhr = static_cast<ghrEntry *>(indirect_history);
    threadInfo[tid].ghr = *previousGhr;

    delete previousGhr;
}

void BLBP::recordTarget(InstSeqNum seq_num, void * indirect_history, const TheISA::PCState& target, ThreadID tid) {
    // here ghr was appended one more
    int ghr_last = threadInfo[tid].ghr[0] | 1;
    shiftGhr_right(threadInfo[tid].ghr);

    // mainly associate target
    ThreadInfo &t_info = threadInfo[tid];
    auto hist_entry = *(t_info.pathHist.rbegin());

    // Temporarily pop it off the history
    // TODO: check if essential

    // first get the btb_entry
    unsigned btb_entry = addr_fold_btb(hist_entry.pcAddr);
    BTBEntry it;
    uint64_t available = 0xffffffffffffffffUL;
    while (!tmp_quque.empty()) tmp_quque.pop(); // we expect an empty tmp queue
    while (!IBTB[tid][btb_entry].empty()) {
        it = IBTB[tid][btb_entry].top();
        IBTB[tid][btb_entry].pop();
        available &= IBTB[tid][btb_entry].targetAddr;
        tmp_quque.push(it);
    }

    available >>= 2; // skip the lowest two bits
    // now we have common bits kept 1 and non-common ones 0
    

    // get common bits from btb_entry

    // update non-common bits for each bits
}

unsigned BLBP::transferFunc(unsigned num) {
    return transferTable[num];
}

std::pair<int, int> BLBP::intervalFunc(unsigned num) {
    return intervalTable[num];
}

// we need to fold br_addr to 6 bits
unsigned BLBP::addr_fold_btb(Addr br_addr) {
    unsigned ret = 0;
    // fold 64 bit to 6 bit
    for (int i = 0; i < 11; ++i) {
        ret = ret ^ (br_addr & 0x3f);
        br_addr >>= 6;
    }
    return ret & 0x3f;
}

unsigned BLBP::history_to_entry(int predictor, ghrEntry &ghr) {
    unsigned ret = 0;
    // expected predictor ranging from 1 to 7
    std::pair<int, int> interval = intervalFunc(predictor - 1);
    // we need to fold values inside this interval to numLocalHistBits
    int i = 0;
    ghrEntry tmp_ghr = ghr;
    for (; i < interval.first; ++i) {
        shiftGhr_right(tmp_ghr); // skip bits that we do not need
    }
    unsigned tmp;
    for (; i + numLocalHistBits < interval.second; ++i) {
        tmp = 0;
        for (int j = 0; j < numLocalHistBits; ++j) {
            tmp |= tmp_ghr.ghr[0] & 1;
            shiftGhr_right(tmp_ghr);
        }
        ret ^= tmp;
    }
    tmp = 0;
    for (;i < interval.second; ++i) {
        tmp |= tmp_ghr.ghr[0] & 1;
        shiftGhr_right(tmp_ghr);
    }
    ret ^= tmp;
    return ret;
}