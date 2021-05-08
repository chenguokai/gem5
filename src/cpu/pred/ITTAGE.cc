//
// Created by xim on 5/8/21.
//

#include "cpu/pred/ITTAGE.hh"

#include "base/intmath.hh"
#include "debug/Indirect.hh"

ITTAGE::ITTAGE(
        const ITTAGEParams &params)
        : IndirectPredictor(params),
          numPredictors(params.numPredictors),
          ghrNumBits(params.indirectGHRBits),
          ghrMask((1 << params.indirectGHRBits)-1)
{
    threadInfo.resize(params.numThreads);

    targetCache.resize(params.numThreads);
    for (unsigned i = 0; i < numSets; i++) {
        targetCache[i].resize(numPredictors);
    }
    use_alt = 8;

    fatal_if(ghrNumBits > (sizeof(ThreadInfo::ghr)*8), "ghr_size is too big");
}

void
ITTAGE::genIndirectInfo(ThreadID tid,
                                         void* & indirect_history)
{
    // record the GHR as it was before this prediction
    // It will be used to recover the history in case this prediction is
    // wrong or belongs to bad path
    indirect_history = new unsigned(threadInfo[tid].ghr);
}

void
ITTAGE::updateDirectionInfo(
        ThreadID tid, bool actually_taken)
{
    threadInfo[tid].ghr <<= 1;
    threadInfo[tid].ghr |= actually_taken;
    threadInfo[tid].ghr &= ghrMask;
}

void
ITTAGE::changeDirectionPrediction(ThreadID tid,
                                                   void * indirect_history, bool actually_taken)
{
    unsigned * previousGhr = static_cast<unsigned *>(indirect_history);
    threadInfo[tid].ghr = ((*previousGhr) << 1) + actually_taken;
    threadInfo[tid].ghr &= ghrMask;
    // maybe we should update hash here?
    // No: CSRs are calculated at use-time
}

bool
ITTAGE::lookup_helper(Addr br_addr, TheISA::PCState& target, ThreadID tid, int &predictor, int &predictor_index)
{
    // todo: adjust according to ITTAGE
    DPRINTF(Indirect, "Looking up %x (set:%d)\n", br_addr, set_index);


    unsigned index = getAddrFold(br_addr);
    int pred_counts = 0;
    TheISA::PCState target1, target2;
    int predictor_1, predictor_2, predictor_index_1, predictor_index_2;

    for (int i = numPredictors - 1; i >= 0; --i) {
        unsigned csr1 = getCSR1(threadInfo[tid].ghr, i);
        unsigned csr2 = getCSR2(threadInfo[tid].ghr, i);
        unsigned tmp_index = index ^ csr1;
        unsigned tmp_tag = (br_addr & 0xff) ^ csr1 ^ (csr2 << 1);
        if (targetCache[tid][i][tmp_index].tag == tmp_tag)) {
            if (pred_counts == 0) {
                target_1 = targetCache[tid][i][tmp_index].target;
                predictor_1 = i;
                predictor_index_1 = tmp_index;
                ++pred_counts;
            }
            if (pred_counts == 1) {
                target_2 = targetCache[tid][i][tmp_index].target;
                predictor_2 = i;
                predictor_index_2 = tmp_index;
                break;
            }
        }
    }
    // decide whether use altpred or not
    if (pred_counts > 0) {
        if (use_alt > 7 && targetCache[tid][predictor_1][predicotr_index_1].counter == 1 &&
                targetCache[tid][predictor_1][predicotr_index_1].useful == 0 && pred_counts == 2 &&
                targetCache[tid][predictor_2][predictor_index_2].counter > 0
        ) {
            target = target_2;
            predictor = predictor_2;
            predictor_index = predictor_index_2;
        } else {
            target = target_1;
            predictor = predictor_1;
            predictor_index = predictor_index_1;
        }
        DPRINTF(Indirect, "Hit %x (target:%s)\n", br_addr, target);
        return true;
    }

    DPRINTF(Indirect, "Miss %x\n", br_addr);
    return false;
}

bool ITTAGE::lookup(Addr br_addr, TheISA::PCState& target, ThreadID tid) {
    int predictor, predictor_index; // no use
    return lookup_helper(br_addr, target, tid, predictor, predictor_index);
}

void
ITTAGE::recordIndirect(Addr br_addr, Addr tgt_addr,
                                        InstSeqNum seq_num, ThreadID tid)
{
    DPRINTF(Indirect, "Recording %x seq:%d\n", br_addr, seq_num);
    HistoryEntry entry(br_addr, tgt_addr, seq_num);
    threadInfo[tid].pathHist.push_back(entry);
}

void
ITTAGE::commit(InstSeqNum seq_num, ThreadID tid,
                                void * indirect_history)
{
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

void
ITTAGE::squash(InstSeqNum seq_num, ThreadID tid)
{
    DPRINTF(Indirect, "Squashing seq:%d\n", seq_num);
    ThreadInfo &t_info = threadInfo[tid];
    auto squash_itr = t_info.pathHist.begin();
    while (squash_itr != t_info.pathHist.end()) {
        if (squash_itr->seqNum > seq_num) {
            break;
        }
        ++squash_itr;
    }
    if (squash_itr != t_info.pathHist.end()) {
        DPRINTF(Indirect, "Squashing series starting with sn:%d\n",
                squash_itr->seqNum);
    }
    t_info.pathHist.erase(squash_itr, t_info.pathHist.end());
}

void
ITTAGE::deleteIndirectInfo(ThreadID tid, void * indirect_history)
{
    unsigned * previousGhr = static_cast<unsigned *>(indirect_history);
    threadInfo[tid].ghr = *previousGhr;

    delete previousGhr;
}

void
ITTAGE::recordTarget(
        InstSeqNum seq_num, void * indirect_history, const TheISA::PCState& target,
        ThreadID tid)
{
    // todo: adjust according to ITTAGE
    ThreadInfo &t_info = threadInfo[tid];

    // Should have just squashed so this branch should be the oldest
    auto hist_entry = *(t_info.pathHist.rbegin());
    // Temporarily pop it off the history so we can calculate the set
    t_info.pathHist.pop_back();

    // we have lost the original lookup info, so we need to lookup again
    int predictor, predictor_index;
    TheISA::PCState dontcare_tgt;
    bool predictor_found = lookup_helper(hist_entry.pcAddr, dontcare_tgt, tid, predicto, predictor_index);

    // update global history anyway
    hist_entry.targetAddr = target.instAddr();
    t_info.pathHist.push_back(hist_entry);

    bool allcate_values = true;

    if (predictor_found && dontcare_tgt == target) {
        // the prediction was from predictor tables and correct
        // increment the counter
        if (targetCache[tid][predictor][predictor_index].counter <= 2) {
            ++targetCache[tid][predictor][predictor_index].counter;
        }
    } else {
        bool bimodal_correct; // todo: add bimodal algo here

        // either a misprediction or no prediction
        if (bimodal_correct) {
            targetCache[tid][predictor][predictor_index].useful = 1;
        }
    }

    // unsigned * ghr = static_cast<unsigned *>(indirect_history);
    // Addr set_index = getSetIndex(hist_entry.pcAddr, *ghr, tid);
    // Addr tag = getTag(hist_entry.pcAddr);


    // assert(set_index < numSets);

    // auto &iset = targetCache[set_index];
    /*
    for (auto way = iset.begin(); way != iset.end(); ++way) {
        if (way->tag == tag) {
            DPRINTF(Indirect, "Updating Target (seq: %d br:%x set:%d target:"
                              "%s)\n", seq_num, hist_entry.pcAddr, set_index, target);
            way->target = target;
            return;
        }
    }

    DPRINTF(Indirect, "Allocating Target (seq: %d br:%x set:%d target:%s)\n",
            seq_num, hist_entry.pcAddr, set_index, target);
    // Did not find entry, random replacement
    auto &way = iset[rand() % numWays];
    way.tag = tag;
    way.target = target;
    */
}

/*
inline Addr
SimpleIndirectPredictor::getSetIndex(Addr br_addr, unsigned ghr, ThreadID tid)
{
    ThreadInfo &t_info = threadInfo[tid];

    Addr hash = br_addr >> instShift;
    if (hashGHR) {
        hash ^= ghr;
    }
    if (hashTargets) {
        unsigned hash_shift = floorLog2(numSets) / pathLength;
        for (int i = t_info.pathHist.size()-1, p = 0;
             i >= 0 && p < pathLength; i--, p++) {
            hash ^= (t_info.pathHist[i].targetAddr >>
                                                   (instShift + p*hash_shift));
        }
    }
    return hash & (numSets-1);
}

inline Addr
SimpleIndirectPredictor::getTag(Addr br_addr)
{
    return (br_addr >> instShift) & ((0x1<<tagBits)-1);
}
*/
inlien int ITTAGE::getTableGhrLen(int table) {
    return 8 << (table - 1);
}

inline unsigned ITTAGE::getCSR1(unsigned ghr, int table) {
    int ghrLen = getTableGhrLen(table);
    unsigned ret = 0, mask = 0x7f;
    int i = 0;
    while (i + 7 < ghrLen) {
        ret = ret ^ (ghr & mask);
        ghr >>= 7;
        i += 7;
    }
    ret = ret ^ (ghr & mask);
    return ret & mask;
}

inline unsigned ITTAGE::getCSR2(unsigned ghr, int table) {
    int ghrLen = getTableGhrLen(table);
    unsigned ret = 0, mask = 0xff;
    int i = 0;
    while (i + 8 < ghrLen) {
        ret = ret ^ (ghr & mask);
        ghr >>= 8;
        i += 8;
    }
    ret = ret ^ (ghr & mask);
    return ret & mask;
}

uint8_t ITTAGE::getAddrFold(int address) {
    uint8_t folded_address, k;
    folded_address = 0;
    for (k = 0; k < 3; k++) {
        folded_address ^= ((address % (1 << ((k + 1) * 8))) / (1 << (k * 8)));
    }
    folded_address ^= address / (1 << (24));
    return folded_address;
}