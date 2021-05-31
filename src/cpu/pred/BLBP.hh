//
// Created by xim on 5/16/21.
//

/*
 * This file provides a basic implementation of BLBP
 * Note that we may have different parameters for BLBP and ITTAGE (and VPC if applicable)
 *
 * */

#ifndef __CPU_PRED_BLBP_HH__
#define __CPU_PRED_BLBP_HH__

#include <deque>
#include <queue>
#include "config/the_isa.hh"
#include "cpu/inst_seq.hh"
#include "cpu/pred/indirect.hh"
#include "params/BLBP.hh"

class BLBP : public IndirectPredictor {
public:
    BLBP(const BLBPParams &params);

    bool lookup(Addr br_addr, TheISA::PCState& br_target, ThreadID tid);
    void recordIndirect(Addr br_addr, Addr tgt_addr, InstSeqNum seq_num, ThreadID tid);
    void commit(InstSeqNum seq_num, ThreadID tid, void * indirect_history);
    void squash(InstSeqNum seq_num, ThreadID tid);
    void recordTarget(InstSeqNum seq_num, void * indirect_history, const TheISA::PCState &target, ThreadID tid);
    void genIndirectInfo(ThreadID tid, void* & indirect_history);
    void updateDirectionInfo(ThreadID tid, bool actually_taken);
    void deleteIndirectInfo(ThreadID tid, void * indirect_history);
    void changeDirectionPrediction(ThreadID tid, void * indirect_history, bool actually_taken);

private:
    bool lookup_helper(Addr, RiscvISA::PCState&, RiscvISA::PCState&, ThreadID, int&, int&, int&, int&, int&, bool&);
    uint8_t getAddrFold(int address);
    const unsigned ghrMask;
    const unsigned pathLength;
    const unsigned numPredictors;
    const unsigned ghrNumBits;
    const unsigned numTageBits;
    const unsigned numWeightBits = 10; // a little longer than paper, the predicted PC length in bit.
    const unsigned numLocalHistBits = 10; // according to the local histroy description
    const unsigned numWeightEntry = (1 << numWeightBits);
    const unsigned numBTBAssoc = 64;
    const unsigned numBTBEntry = 64;

    unsigned transferFunc(unsigned);
    std::pair<int, int> intervalFunc(unsigned);

    uint64_t global_counter = 0; // the global counter for LRU

    struct BTBEntry {
        BTBEntry() : tag(0), target(0), counter(0), useful(0) {}
        Addr tag;
        TheISA::PCState target;
        uint64_t counter; // used for LRU
        bool operator <(const BTBEntry &a) const { return counter > a.counter; } 
        // priority_queue provide the max value by default, here we want to get the min one
    };

    std::vector<std::vector<std::priority_queue<BTBEntry> > >  IBTB; // indirect BTB, level 1: thread, level 2: address by hist level 3 priority_queue: help retire
    std::queue<BTBEntry> tmp_quque; // priority queue does not provide iterator, iterate by ourselves.
    struct HistoryEntry {
        HistoryEntry(Addr br_addr, Addr tgt_addr, InstSeqNum seq_num) : pcAddr(br_addr), targetAddr(tgt_addr), seqNum(seq_num) {}
        Addr pcAddr;
        Addr targetAddr;
        InstSeqNum seqNum;
    };

    struct ghrEntry {
        ghrEntry() :ghr({0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) {}
        uint64_t ghr[10];
    }

    struct ThreadInfo {
        ThreadInfo() : headHistEntry(0), ghr(0) {}

        std::deque<HistoryEntry> pathHist;
        unsigned headHistEntry;
        ghrEntry ghr;
    };
    struct Weight{
        uint8_t weight[numWeightBits];
    }

    std::vector<ThreadInfo> threadInfo;
    std::vector<std::vector<std::vector<Weight> > >WeightTable; // level 1: thread level 2: predictor level 3: entry
    // std::vector<>
    const unsigned transferTable = [2, 4, 6, 8, 12, 14, 17, 24]; // approximated value from paper
    std::vector<std::pair<int, int> > intervalTable;

    // local history: record the indirect local history
    // maybe it refers to both indirect and conditional here, but we have no way to determin for now.
    // 10 bit
    const unsigned numLocalHistoryEntryBits = 8;
    const unsigned numLocalHistoryEntry = (1 << numLocalHistoryEntryBits);

    // note that we do not have any efficient way to recover, so we actually reserve more than needed history here
    std::vector<std::vector<uint64_t> > LocalHistory; // level 1: thread level 2: 256 different entries, mappped by pc


#endif __CPU_PRED_BLBP_HH__
