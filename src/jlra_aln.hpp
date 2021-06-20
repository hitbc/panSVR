/*
 * jlra_aln.hpp
 *
 *  Created on: 2021年3月10日
 *      Author: fenghe
 */

#ifndef JLRA_ALN_HPP_
#define JLRA_ALN_HPP_

#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cstring>
#include <ctype.h>

#include <vector>
#include <map>
#include "deBGA_index.hpp"
#include "cpp_lib/get_option_cpp.hpp"
extern "C"
{
#include "kswlib/kalloc.h"
#include "kswlib/ksw2.h"

#include "clib/bam_file.h"
#include "clib/desc.h"
#include "clib/binarys_qsort.h"
}

#define SEED_OFFSET (K_T - LEN_KMER)
#define LEN_KMER_LEFT (32 - LEN_KMER) //LEN_KMER_LEFT = 32 - LEN_KMER
#define SEED_STEP 5
#define UNI_POS_N_MAX 32

#define SEED_STEP 5
#define MIN_CHAIN_SCORE 20

#define GAP_OPEN 16
#define GAP_EXT 1
#define GAP_OPEN2 32
#define GAP_EXT2 0
#define MATCH_SCORE 2
#define MISMATCH_SCORE 12

#define ZDROP_SCORE 400
#define BANDWIDTH 500

#define OUTPUT "./aln.sam"

struct MAP_PARA{
	//basic option
	int thread_n;

	//score option
	int match_D;
	int mismatch_D;
	int gap_open_D;
	int gap_ex_D;
	int gap_open2_D;
	int gap_ex2_D;
	int bw;
	int zdrop_D; //for DNA zdrop = 400, 200 for RNA

	int max_use_read;
	//alignment option

	//output option
	bool NOT_output_original_rst;//when set to be true, NOT output original result
	char *sam_path;//output file name
	bool output_sam;

	//input option
	char * indexDir;
	char * read_fastq1;
	char * ori_header_fn;

	int get_option(int argc, char *argv[]){
		options_list l;

	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  aln  [Options] <IndexDir> [ReadFiles.fa][ori_header_fn.sam]>\n");
	    l.add_title_string("  Basic:   \n");
	    l.add_title_string("    <IndexDir>      FOLDER   the directory contains index\n");
	    l.add_title_string("    [ReadFiles.fa]  FILES    reads files, FASTQ(A)(or fa.gz/fq.gz) format only one file accepted, \n");
	    l.add_title_string("                             for pair end NGS read, read 1 and 2 in a pair should store together.\n");
	    l.add_title_string("                             Using [signal] command to generate this type of file\n");
	    l.add_title_string("    [ori_header.sam]  FILES  Header file of original BAM/CRAM file, using [signal] command or [samtools view -H] to generate it\n");

		//thread number
		l.add_option("thread", 			't', "Number of threads", true, 4); l.set_arg_pointer_back((void *)&thread_n);
		//score
		l.add_option("gap-open1", 		'O', "Gap open penalty 1", true, GAP_OPEN); l.set_arg_pointer_back((void *)&gap_open_D);
			l.l.back().add_help_msg("a k-long gap costs min{O+k*E,P+k*F}..");
		l.add_option("gap-open2", 		'P', "Gap open penalty 2.", true, GAP_OPEN2); l.set_arg_pointer_back((void *)&gap_open2_D);
		l.add_option("gap-extension1", 	'E', "Gap extension penalty 1.", true, GAP_EXT);l.set_arg_pointer_back((void *)&gap_ex_D);
		l.add_option("gap-extension2", 	'F', "Gap extension penalty 2.", true, GAP_EXT2);l.set_arg_pointer_back((void *)&gap_ex2_D);
		l.add_option("match-score", 	'M', "Match score for SW-alignment.", true, MATCH_SCORE);l.set_arg_pointer_back((void *)&match_D);
		l.add_option("mis-score", 		'm', "Mismatch score for SW-alignment.", true, MISMATCH_SCORE);l.set_arg_pointer_back((void *)&mismatch_D);
		l.add_option("zdrop", 			'z', "Z-drop score for splice/non-splice alignment.", true, ZDROP_SCORE);l.set_arg_pointer_back((void *)&zdrop_D);
		l.add_option("band-width", 		'w', "Bandwidth used in chaining and DP-based alignment.", true, BANDWIDTH);l.set_arg_pointer_back((void *)&bw);
		//output
		l.add_option("output", 			'o', "Output file (BAM format)", true, "./output.bam");l.set_arg_pointer_back((void *)&sam_path);
		l.add_option("not-ori",			'Q', "when set to be true, NOT output original result when score of ORI is bigger");l.set_arg_pointer_back((void *)&NOT_output_original_rst);
		l.add_option("SAM", 			'S', "Output as SAM instead of (BAM format)"); l.set_arg_pointer_back((void *)&output_sam);
		l.add_option("max_use_read", 	'R', "Max number of read to alignment, used for debug or test PG.", true, MAX_int32t); l.set_arg_pointer_back((void *)&max_use_read);

		if(l.default_option_handler(argc, argv)) return 1;

		if (argc - optind < 3)
			return l.output_usage();
		xassert((thread_n >= 1) && (thread_n <= 48), "Input error: thread_n cannot be less than 1 or more than 48\n");

		indexDir = strdup(argv[optind]);
		if (indexDir[strlen(indexDir) - 1] != '/') strcat(indexDir, "/");
		read_fastq1 = strdup(argv[optind + 1]);
		ori_header_fn = strdup(argv[optind + 2]);
		return 0;
	}
};

#define OUT_BAM_CIGAR_STR   "MIDNSHP=XB"
struct CIGAR_PATH{
	CIGAR_PATH(char type_, uint16_t size_){
		switch(type_){
		case 'M': type = 0; break;
		case 'I': type = 1; break;
		case 'D': type = 2; break;
		case 'N': type = 3; break;
		case 'S': type = 4; break;
		case 'H': type = 5; break;
		case 'P': type = 6; break;
		case '=': type = 7; break;
		case 'X': type = 8; break;
		case 'B': type = 9; break;
		default : type = 0; break;
		}
		size = size_;
	}
	CIGAR_PATH(uint32_t bin_cigar){
		type = (int)((bin_cigar & 0xf));
		size = (bin_cigar >> 0x4);
	}
	CIGAR_PATH(const CIGAR_PATH &cp){
		type = cp.type;
		size = cp.size;
	}

	//return true when merge to the previous node; false when not merge and need to new a node
	bool try_merge(const CIGAR_PATH &cp){
		if(cp.size < 0){
			xassert(cp.type == 2, ""); //deletion
			if(type == 0){
				size += cp.size;
				xassert(size > 0, "");
				return true;
			}else if(type == 2){
				size -= cp.size;
				xassert(size > 0, "");
				return true;
			}else{
				xassert(0, "");
			}
		}else if(type == cp.type || cp.size == 0){
			size += cp.size;
			return true;
		}
		return false;
	}

	uint8_t type;
	int16_t size;
};

struct KSW_ALN_handler{

#define KSW_ALN_left_extend 0
#define KSW_ALN_right_extend 1
#define KSW_ALN_end_to_end 2
//#define KSW_ALN_left_extend_forward_cigar 3

	//read info
	uint8_t *read_str;
	//uint32_t read_cigar[100];
	//uint32_t n_cigar;
	int32_t read_score;
	uint32_t total_q_len;

	//sequence buff
	uint8_t *qseq;
	uint8_t *tseq;
	uint8_t *qseq_rev;
	uint32_t qlen;
	uint32_t tlen;
	uint8_t type;

	//result
	bool is_simple_aln;
	uint32_t simple_NM;
	ksw_extz_t ez;

	//kmalloc
	void *km;

	//mapping options
	int8_t mata_D[25];

	int8_t match_D;
	int8_t mismatch_D;
	int8_t gap_open_D;
	int8_t gap_ex_D;
	int8_t gap_open2_D;
	int8_t gap_ex2_D;
	uint16_t zdrop_D; //for DNA zdrop = 400, 200 for RNA
	int bandwith;

	int flag;
	std::vector<CIGAR_PATH> cigar_tmp;

	deBGA_INDEX * idx;

	void copy_option(MAP_PARA 	*o);
	void ksw_gen_mat_D();
	void init(MAP_PARA 	*o, deBGA_INDEX * idx_);
	void free_memory();
	void setRead(uint8_t *read_str_);
	void alignment(int read_st_, int read_ed_, int ref_st_, int ref_ed_, int type_);
	void align_non_splice();
	int get_misMatch(int read_st_, int read_ed_, int ref_st_, int ref_ed_);

};

//used to store temporary alignment results
struct MAX_IDX_OUTPUT{
	//alignment scores
	uint32_t align_score;
	uint32_t chain_score;

	//position in read
	uint32_t max_index;
	uint32_t read_bg;

	//sv information that result belong to
	//uint32_t sv_ID;//for a new alignment, which SV region it belonged to
	SV_chr_info *sv_info_p;

	//mapq
	uint8_t mapq;

	//mate information
	bool has_mate;
	uint32_t mate_chrID;
	uint32_t mate_ref_bg;
	SV_chr_info *mate_sv_info_p;

	//position in reference
	bool is_ori;
	uint32_t chrID;
	uint32_t ref_bg;
	int direction;

	//alignment cigar
	std::vector<CIGAR_PATH> cigar;

	//index in results
	int rst_idx;

    int reverseGIGAR(std::vector<CIGAR_PATH> &cigar_tmp, int read_len){
		cigar.clear();
		cigar.emplace_back(cigar_tmp.back());
		for(int i = cigar_tmp.size() - 2; i >= 0 ; i--){
				if(cigar.back().try_merge(cigar_tmp[i]) == false)
						cigar.emplace_back(cigar_tmp[i]);
		}
		if(!cigar.empty() && cigar[0].size == 0)
				cigar.erase(cigar.begin());
		//check:
		int total_read_len = 0;
		for(auto & ci: cigar){
			if(  ci.type == 0 || ci.type == 1 || ci.type == 3 || ci.type == 4){	total_read_len += ci.size;	}
		}
		if(total_read_len != read_len){
			fprintf(stderr, "ERROR cigar: read_len: %d\t", read_len);
			print(stderr, '\t');
			for(auto & c : cigar){
					fprintf(stderr, "%d%c", c.size, OUT_BAM_CIGAR_STR[c.type]);
			}
			fprintf(stderr, "\n");
			return false;
		}
		return true;
    }

	static int cmp_chain_score(const void*a_,const void*b_)	{
		MAX_IDX_OUTPUT *a = (MAX_IDX_OUTPUT *)a_;
		MAX_IDX_OUTPUT *b = (MAX_IDX_OUTPUT *)b_;
		if(a->chain_score != b->chain_score)	return a->chain_score < b->chain_score;
		else									return a->max_index > b->max_index;
	}

	static int cmp_align_score(const void*a_,const void*b_)	{
		MAX_IDX_OUTPUT *a = (MAX_IDX_OUTPUT *)a_;
		MAX_IDX_OUTPUT *b = (MAX_IDX_OUTPUT *)b_;
		if(a->align_score != b->align_score)	return a->align_score < b->align_score;
		else									return a->max_index > b->max_index;
	}

	void print(FILE * f, char endl){ fprintf(f, "chain_score:[%d], align_score:[%d], read_bg:[%d], chrID:[%d], ref_bg:[%d], %c", chain_score, align_score, read_bg, chrID, ref_bg, endl); }

};

//used for alignment of single end read
#define MAX_READ_LEN 1600
#define MAX_OUTPUT_NUMBER 6
struct single_end_handler{
	random_data rand_buff;
	//result
	int result_num;
	MAX_IDX_OUTPUT result[MAX_OUTPUT_NUMBER + MAX_OUTPUT_NUMBER];

	MAX_IDX_OUTPUT * primary_result;
	MAX_IDX_OUTPUT * secondary_result;

	//original result
	MAX_IDX_OUTPUT ori;
	bool ORI_is_UNMAPPED;

	void init(MAP_PARA 	*o, deBGA_INDEX * idx_){
		//random
		char * c_statebuf = (char *)xcalloc(128, 1);
		initstate_r(rand(), c_statebuf, 128, &rand_buff);

		kswh.init(o, idx_);
		s.s = sam_string;
		s.m = 2048;
		NOT_output_original_rst = o->NOT_output_original_rst;
		idx = idx_;
	}

	void destory(){	kswh.free_memory(); }
	// step 1:
	void read_register(kseq_t * c_read_){
		c_read = c_read_;
		read_l = c_read->seq.l;
		read_seq_char = c_read->seq.s;
	}

	// step2 :aligning read into reference and store all result in [result]，total result number stored in [result_num], and
	// primary ID and second ID stored in 	[int primary_index;] and [int secondary_index;]
	void align();
	void output_BAM(bam1_t *b, bool is_first_read, bool read_proper_mate, char * comment_str);//step3: store result in bam1_t *b

	//private function and variable:
private:

	uint64_t read_l;
	char * read_seq_char;
	kseq_t * c_read;
	uint8_t bin_read[2][MAX_READ_LEN];

	//buff:
	KSW_ALN_handler kswh;
	Graph_handler g[2];
	std::vector<vertex_MEM> vertexm_v;
	std::vector<vertex_U> vertexu_v;
	std::vector<UNI_SEED> uniseed_v[2];//for read 1 and read 2; for forward and reverse
	uint8_t repeat_seed_info[MAX_READ_LEN];
	std::map<uint64_t, int> kmer_set_buff;
	bool readIsSTR;

	//output buff
	char sam_string[2048];
	kstring_t s;

	//flags
	bool NOT_output_original_rst;

	//index
	deBGA_INDEX * idx;

	//ori_rst[0] is the original result of read 1 and ori_rst[1] is the second read
	void parse_ori_mapping_rst(char * readComment, int read_len){
		//set result is original
		ori.is_ori = true;
		//set basic informations
		char *multy_safty_strtok_ptr = NULL;
		int str_len_comment = strlen(readComment);

		char * token = strtok_r(readComment, "_", &multy_safty_strtok_ptr);		ori.chrID = atoi(token); //tid
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	ori.ref_bg = atoi(token); //ori_rst.ref_ed = ori_rst.ref_bg;
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	ori.read_bg = atoi(token); //soft_left[i]
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	ori.align_score = atoi(token); // score[i]
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	ori.mapq = atoi(token); // mapq[i]

		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	// mate mapq
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	// XA_number[i]
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); 	// mate XA_number[i]
		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); //isize

		token = strtok_r(NULL, "_", &multy_safty_strtok_ptr); //signal flags
		ori.direction = (token[0] == 'F')?FORWARD:REVERSE; //FORWARD?REVERSE
		ORI_is_UNMAPPED = (token[1] == 'Y')?true:false; // unmapped?mapped
		//store cigar
		ori.cigar.clear();
		if(ori.read_bg > 0)
			ori.cigar.emplace_back('S', ori.read_bg);
		ori.cigar.emplace_back('M', read_len - ori.read_bg);
		//store sv_info
		ori.sv_info_p = NULL;
		ori.has_mate = false;

		if(ori.ref_bg >= MAX_int32t)
			ori.ref_bg = 1;
		//restore comment
		for(int i = 0; i < str_len_comment - 1; i++)
			if(readComment[i] == 0)
				readComment[i] = ',';

	}
	void chainning_one_read(int read_is_reverse);
	void binary_read_2_bit();
};

#define K_T 22 		//kmer length used for building deBGA index

//data used for each ALIGN-PROCESSIONG thread
struct Classify_buff_pool
{
	single_end_handler SE_h[2];
};

//data used for PIPELINE each threads
struct CLASSIFY_THREAD_DATA
{
	kseq_t	*	seqs1; // first  read in a pair
	kseq_t	*	seqs2; // second read in a pair
	//OUTPUT_RST * rst; 	//mapping result for read1 and read2
	bam1_t 	*b1;
	bam1_t 	*b2;
	int 		readNum;
	void 	*	share_data_pointer;// register for the shared data
};

//data SHARED by all threads
struct CLASSIFY_SHARE_DATA
{
	//shared
	deBGA_INDEX * idx = NULL;//deBGA index
	kstream_t 	*_fp1 = NULL; //fastq files 1

	MAP_PARA 	*o; 		//mapping option
	Classify_buff_pool 	*buff = NULL;//data used for each classify thread
	CLASSIFY_THREAD_DATA *data = NULL;//for each pipeline thread

	htsFile *output_file = NULL;//result output file
};

struct CLASSIFY_MAIN{

	CLASSIFY_SHARE_DATA *share = NULL;

	void init_run(int argc, char *argv[]);
private:
	static void *classify_pipeline(void *shared, int step, int tid, void *_data);									//pipeline
	static int load_reads(kstream_t *_fp1, kseq_t *_seqs1, kseq_t *_seqs2, int n_needed, int MAX_load_read_number);	//function in step 1
	static void inline worker_for(void *_data, long data_index, int thread_index); 									//function in step 2
	static void output_results(int readNum, deBGA_INDEX * idx, htsFile *output_file, bam1_t *b1, bam1_t *b2);		//function in step 3
};

#endif /* JLRA_ALN_HPP_ */
