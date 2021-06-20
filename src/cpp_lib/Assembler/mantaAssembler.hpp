/*
 * assembly.hpp
 *
 *  Created on: 2020年4月24日
 *      Author: fenghe
 */

#pragma once

#include <iosfwd>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include "../../clib/utils.h"

typedef int32_t pos_t;

struct AddReadAction{
	AddReadAction(	int position_, int read_ID_, bool isAdd_){
		setBasic(position_, read_ID_, isAdd_);
	}
	void setBasic(int position_, int read_ID_, bool isAdd_){
		position_in_contig = position_;
		read_ID = read_ID_;
		 isAdd = isAdd_;
	}
	int position_in_contig;
	int read_ID;
	bool isAdd;//true for adding read action; false for delete read action

	//used in contig alignment
	int position_read;
	int suggest_contig_offset_in_ref;
	int search_len;//debugcode
	int wrong_base;

	void set_read_pos(std::string &read_seq, std::string &contig_seq, int ass_begin_offset_in_contig, int word_len, std::set<int> &remove_read_set, int read_in_ref_offset){

		if(isAdd == false){
			remove_read_set.emplace(read_ID);
			return;
		}

		position_read = -1;
		int read_kmer_num = read_seq.size() - word_len + 1;
		search_len = 0;
		if(position_in_contig < 0){
			for(int i = read_kmer_num - 1; i >= 0; i--){
				search_len++;
				int cmp = read_seq.compare(i, word_len, contig_seq, position_in_contig - ass_begin_offset_in_contig, word_len);
				if(cmp == 0){
					position_read = i;
					break;
				}
			}
		}else{
			for(int i = 0; i < read_kmer_num; i++){
				search_len++;
				int cmp = read_seq.compare(i, word_len, contig_seq, position_in_contig - ass_begin_offset_in_contig, word_len);
				if(cmp == 0){
					position_read = i;
					break;
				}
			}
		}

		suggest_contig_offset_in_ref = read_in_ref_offset - (position_in_contig - ass_begin_offset_in_contig - position_read);

		xassert(position_read != -1, "");
		//position_ref -=	ref_begin_pos;
	}
	void print(FILE *output, bool is_from_main_SV) const {
		fprintf(output,
			"[CP: %d RID:%d %c from_MAIN: %s ",
			position_in_contig, read_ID, isAdd?'+':'-', is_from_main_SV?"YES":"NO");

		if(isAdd)
			fprintf(output,
				"RP: %d SL: %d SUG: %d WB %d ",
				position_read, search_len, suggest_contig_offset_in_ref, wrong_base);
		fprintf(output,	"]\n");
	}
};

struct AssembledContig {
	std::string seq;  ///< contigsequence
	unsigned seedReadCount = 0;  ///< no of reads containing the seeding kmer
	std::set<unsigned> supportReads;
	std::set<unsigned> rejectReads;
	pos_t conservativeRange_bgn = 0; ///< subsection of the contig with conservative coverage
	pos_t conservativeRange_end = 0;
	int ending_reason[2]; //index: 0 for left and 1 for right; data : 0 for reason 0[maxBaseCount < o.minCoverage] and 1 for reason 1[ after seeing one repeat word]
	int new_support_read = 0;
	unsigned wordLength;
	std::vector<AddReadAction> actions;
	int ass_begin_offset_in_contig;// the offset of begin of contig compare with the start assembly kmer

	void debug_print(FILE *output) const {

		fprintf(output,
				"CONTIG size: [%ld] seedCount: [%d] supportReads: [%ld] "
				"ending_reason: [%d %d]"
				"seq:\n",
				seq.size(), seedReadCount, supportReads.size(), ending_reason[0], ending_reason[1]);

		const char *seq_ = seq.c_str();
		static const unsigned rowSize(200);
		static const unsigned sectionSize(100);
		assert(nullptr != seq_);
		const unsigned seqLen(strlen(seq_));
		for (unsigned i(0); i < seqLen; ++i) {
			if (i) {
				if (0 == (i % rowSize))
					fprintf(output, "\n");
				else if (0 == (i % sectionSize))
					fprintf(output, " ");
			}
			fprintf(output, "%c", seq_[i]);
		}

		fprintf(output, "\nsupportReads: ");
		for (auto &r : supportReads)
			fprintf(output, "%d ", r);
		fprintf(output, "\n");

		fprintf(output, "\nrejectReads: ");
		for (auto &r : rejectReads)
			fprintf(output, "%d ", r);
		fprintf(output, "\n");
	}
};

/**************************************************************/
/// Information added to each read in the process of assembly
struct AssemblyReadInfo {
	AssemblyReadInfo(bool isPseudo_ = false) :
			isPseudo(isPseudo_) {
	}
	bool isUsed = false;
	/// If true, the read was 'used' but filtered out, so there is no meaningful contig id association
	bool isFiltered = false;
	/// If true, the read was an assembled contig
	bool isPseudo = false;
	/// Index of the contigs that this read is used in
	std::vector<unsigned> contigIds;
};

typedef std::vector<AssembledContig> Assembly;
typedef std::vector<std::string> AssemblyReadInput;
typedef std::vector<AssemblyReadInfo> AssemblyReadOutput;

/// Input parameters for IterativeAssembler
///
struct IterativeAssemblerOptions {
	IterativeAssemblerOptions() {
	}

	/// the symbol set used during assembly
	std::string alphabet = "ACGT";
	/// minimum basecall quality for assembly input
	int minQval = 5;
	/// initial word (kmer) length
	unsigned minWordLength = 26;
	unsigned maxWordLength = 106;
	unsigned maxWordLength_without_enough_read = 66;
	unsigned wordStepSize = 10;
	unsigned minContigLength = 45;
	/// min. coverage required for contig extension
	unsigned minCoverage = 2;
	/// coverage required for conservative contig sub-range
	unsigned minConservativeCoverage = 2;
	/// max error rates allowed during contig extension
	double maxError = 0.35;
	/// min. number of unused reads to enable search for more contigs
	unsigned minUnusedReads = 3;
	/// min. number of reads required to start assembly
	unsigned minSupportReads = 3;
	/// Max. number of assembly returned for a given set of reads
	unsigned maxAssemblyCount = 20;
};

typedef std::unordered_map<std::string, unsigned> str_uint_map_t;
// maps kmers to supporting reads
typedef std::unordered_map<std::string, std::set<unsigned>> str_set_uint_map_t;
typedef std::unordered_map<std::string, std::pair<unsigned, unsigned>> str_pair_uint_map_t;

struct AssemblyManager {

	std::vector<std::string> reads;
	//results
	std::vector<AssembledContig> contigs;
	void addRead(const char *seq) {
		reads.push_back(std::string(seq));
	}
	void clear() {
		reads.clear();
		readInfo.clear();
	}
	std::vector<AssembledContig>& getResults() {
		return contigs;
	}
	void assembley();

private:
	//************************get kmer count***************************************/
	void getKmerCounts(const unsigned wordLength, str_uint_map_t &wordCount,
			str_set_uint_map_t &wordSupportReads);
	//buffs for getKmerCounts
	std::set<std::string> getKmerCounts_readWords_BUFF;
	//****************************getRepeatKmers*******************************************/
	void getRepeatKmers(const str_uint_map_t &wordCount,
			std::set<std::string> &repeatWords);
	//buffs
	str_pair_uint_map_t getRepeatKmers_wordIndices;
	std::vector<std::string> getRepeatKmers_wordStack;

	//**************************buildContigs********************************/
	bool buildContigs(const unsigned wordLength, unsigned & global_maxWordCount);
	//buffs
	str_uint_map_t buildContigs_wordCount;
	// records the supporting reads for each kmer
	str_set_uint_map_t buildContigs_wordSupportReads;
	// identify repeat kmers (i.e. circles from the de bruijn graph)
	std::set<std::string> buildContigs_repeatWords;
	std::set<std::string> buildContigs_unusedWords;
	/************************selectContigs*******************************/
	void selectContigs(const unsigned normalReadCount);
	//buffs
	std::set<unsigned> selectContigs_usedReads;
	std::set<unsigned> selectContigs_usedPseudoReads;
	std::set<unsigned> selectContigs_contigs2Remove;
	std::set<unsigned> selectContigs_newSupportReads;

	IterativeAssemblerOptions o;
	std::vector<AssemblyReadInfo> readInfo;
	//buffs
	std::vector<AssembledContig> tmpContigs;

	/***********************************Walk**************************************/
	bool walk(const std::string &seed, const unsigned wordLength,
			const str_uint_map_t &wordCount,
			const str_set_uint_map_t &wordReads,
			const std::set<std::string> &repeatWords,
			std::set<std::string> &unusedWords, AssembledContig &contig);
	//buffs for walk
	std::set<std::string> walk_wordsInContig;
	std::set<unsigned> walk_maxWordReads;
	std::set<unsigned> walk_maxContigWordReads;
	std::set<unsigned> walk_previousWordReads;
	std::set<unsigned> walk_supportReads2Remove;
	std::set<unsigned> walk_rejectReads2Add;
	std::set<unsigned> walk_contigWordReads;
	std::set<unsigned> walk_sharedReads;
	std::set<unsigned> walk_toRemove;
	std::set<unsigned> walk_toAdd;
	std::set<unsigned> walk_sharedReads_alleles;
	std::set<unsigned> walk_toUpdate;
};

//demo:
/**
 * void main()
 * {
 * 		AssemblyManager am;
 * 		while(1)
 * 		{
 * 			if(END_of_all_ass_event)
 * 				break;
 * 			am.clear();
 * 			for(read in read_list)
 * 				am.addRead(read);
 * 			am.assembley();
 * 			auto &contigs = am.getResults();
 * 			for(auto &contig:contigs)
 * 				std::cout << contig;
 * 		}
 * }
 *
 */

void assembly_test();
