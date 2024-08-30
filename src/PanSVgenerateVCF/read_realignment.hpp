
#ifndef READ_REALIGNMENT_HPP_
#define READ_REALIGNMENT_HPP_

#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <cstring>
#include <ctype.h>

#include <vector>
#include <map>
#include "deBGA_index.hpp"
#include "../cpp_lib/get_option_cpp.hpp"
extern "C"
{
#include "../kswlib/kalloc.h"
#include "../kswlib/ksw2.h"

#include "../clib/bam_file.h"
#include "../clib/desc.h"
#include "../clib/binarys_qsort.h"
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
	char *sam_path_signal_ori;//output file name, used to store reads signals that not aligned by aligner nor pan-genome based realigner

	bool output_sam;

	//input option
	char * indexDir;
	char * read_fastq1;
	char * ori_header_fn;

	//read status data
	bool read_status_options_already_set;// == false
	int ISIZE_MIN;
	int ISIZE_MID;
	int ISIZE_MAX;
	int normal_read_length;

	int get_option(int argc, char *argv[]){
		options_list l;
		l.show_command(stderr, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  fc_aln  [Options] <IndexDir> [ReadFiles.fa][ori_header_fn.sam]>\n");
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
		l.add_option("output_signal_ori",'p', "Reads not fully alignment by aligner nor re-aligner(BAM format)", true, "./output_ori.bam");l.set_arg_pointer_back((void *)&sam_path_signal_ori);

		l.add_option("not-ori",			'Q', "when set to be true, NOT output original result when score of ORI is bigger");l.set_arg_pointer_back((void *)&NOT_output_original_rst);
		l.add_option("SAM", 			'S', "Output as SAM instead of (BAM format)"); l.set_arg_pointer_back((void *)&output_sam);
		l.add_option("max_use_read", 	'R', "Max number of read to alignment, used for debug or test PG.", true, MAX_int32t); l.set_arg_pointer_back((void *)&max_use_read);
		if(l.default_option_handler(argc, argv)) {
			l.show_c_value(stderr);
			return 1;
		}
		l.show_c_value(stderr);

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
		s.m = 204800;
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
	void output_BAM(bam1_t *b, bool is_first_read, int ABS_isize, char * comment_str);//step3: store result in bam1_t *b
	void output_ori_bam(bam1_t *b, int max_score);

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
	char sam_string[204800];
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

struct PE_score{
	//results part
	int max_same;
	int max_score;
	//properly mated
	bool read_pair_is_proper_mated;
	int cur_isize; // current isize

	//bool: this value indicated whether the pan-genome gain better results compared with original results
	bool pan_genome_gain_better_result;
	//bool: this value indicated whether read pair nearly full-match-aligned to the original genome or pan-genome
	//bool read_pair_full_match; // when all reads in a pair has 5 mis-matches or other mutations, this value will be set to true.

	MAX_IDX_OUTPUT *max_1;
	MAX_IDX_OUTPUT *max_2;

	//bam/cram status
	int max_isize; //the max ISZIE the that can be accepted, 99% of reads have ISIZEs shorter than that value
	int min_isize; //the min ISZIE the that can be accepted, 99% of reads have ISIZEs longer than that value
	int normal_read_len;
	int min_filter_score;

	void init(int max_isize_, int min_isize_, int normal_read_len_, int min_filter_score_){
		max_isize = max_isize_ + 200;
		min_isize = min_isize_ - 200;
		min_isize = MAX(0, min_isize);
		normal_read_len = normal_read_len_;
		min_filter_score = min_filter_score_;
		clear();
	}

	void clear(){
		max_same = 1;
		max_score = 0;
		max_1 = NULL;
		max_2 = NULL;
		cur_isize = 0;
		read_pair_is_proper_mated = false;
		//read_pair_full_match = false;
		pan_genome_gain_better_result = false;
	}

	void read_get_best_pairing_results(single_end_handler *SE_h){
		clear();
		//for both original
		int se0_rst_num = SE_h[0].result_num; int se1_rst_num = SE_h[1].result_num;
		if(!SE_h[0].ORI_is_UNMAPPED) se0_rst_num++;
		if(!SE_h[1].ORI_is_UNMAPPED) se1_rst_num++;

		//if(se0_rst_num == 0 )
		//when using only second end
		for(int i = 0; i < se0_rst_num; i++)
			store_pair_end_score(&((i < SE_h[0].result_num)?SE_h[0].result[i]:SE_h[0].ori), NULL);
		//when using only first end
		for(int j = 0; j < se1_rst_num; j++)
			store_pair_end_score(NULL, &((j < SE_h[1].result_num)?SE_h[1].result[j]:SE_h[1].ori));
		//when using both end
		for(int i = 0; i < se0_rst_num; i++){
			MAX_IDX_OUTPUT & se1 = (i < SE_h[0].result_num)?SE_h[0].result[i]:SE_h[0].ori;
			for(int j = 0; j < se1_rst_num; j++)
				store_pair_end_score(&se1, &((j < SE_h[1].result_num)?SE_h[1].result[j]:SE_h[1].ori));
		}
		//discard result when both best single end results are original(or NULL) OR max score is 0
		pan_genome_gain_better_result = pangenome_gain_better_result();
		//if(!pan_genome_gain_better_result) { max_1 = NULL; max_2 = NULL; }
	}

	void set_primary_secondary_mate(single_end_handler *SE_h){
		for(int i = 0; i < 2; i++){
			MAX_IDX_OUTPUT *c_rst = get_cur_output(i);
			if(c_rst == NULL)
				continue;
			single_end_handler * c_SE_h = SE_h + i;
			//set primary results
			c_SE_h->primary_result = c_rst;
			//set secondary
			c_SE_h->secondary_result = NULL;
			if(c_rst->is_ori && c_SE_h->result_num > 0){
				c_SE_h->secondary_result = &(c_SE_h->result[0]);
			}
			else if(c_SE_h->result_num > 1){
				c_SE_h->secondary_result = (c_rst->rst_idx == 0)?&(c_SE_h->result[1]):&(c_SE_h->result[0]);
			}
			//set mate
			MAX_IDX_OUTPUT *mate_rst = get_cur_output(1 - i);
			if(mate_rst != NULL && mate_rst->chrID != MAX_uint32_t){
				c_rst->has_mate = true;
				c_rst->mate_chrID = mate_rst->chrID;
				c_rst->mate_ref_bg = mate_rst->ref_bg;
				c_rst->mate_sv_info_p = mate_rst->sv_info_p;
				//when one end is original and mate read is new alignment, set the original has same SV info as new alignment
				if(c_rst->is_ori)
					c_rst->sv_info_p = c_rst->mate_sv_info_p;
			}
			else{
				c_rst->has_mate = false;
				c_rst->mate_chrID = 0;
				c_rst->mate_sv_info_p = NULL;
			}
		}
	}

private:

	MAX_IDX_OUTPUT *get_cur_output(int read_ID){ return (read_ID == 0)?max_1:max_2; }

	bool pangenome_gain_better_result(){
		if(max_score <= 0)					return false;
		if(max_1 != NULL && !max_1->is_ori) return true;
		if(max_2 != NULL && !max_2->is_ori)	return true;
		return false;
	}

	//a sub-function used in "store_pair_end_score"
	void set_score(int new_score, MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2, int ISIZE){
		bool store_new_best = true;
		if(new_score > max_score) max_same = 1;
		else if(new_score == max_score){
			max_same ++;
			if(rand() % max_same != 0) store_new_best = false;
		}
		if(store_new_best){
			max_1 = se1; max_2 = se2; max_score = new_score; cur_isize = ISIZE; read_pair_is_proper_mated = (cur_isize > 0);
		}
	}

	//return true isize (> 0); when not proper paired, return 0;
	int get_isize(int p1, int p2, int d1, int d2){
		if(d1 == d2) return 0;
		int isize = normal_read_len + ((d1 == FORWARD)?(p2 - p1):(p1 - p2));
		return (isize < max_isize && isize > min_isize)?isize: 0;
	}
	//it is a sub-function for "store_pair_end_score"
	//return isize (>0) when reads are properly mapped; otherwise return 0
	int read_is_proper_mated(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2){
		if(se1 == NULL || se2 == NULL || se1->chrID != se2->chrID){ return 0; }
		//s1_p1 is the position of read 1 in first SV region, s1_p2 is the position in the second region
		int s1_p1 = se1->ref_bg; int s1_p2 = s1_p1 + ((se1->is_ori)?0:se1->sv_info_p->end_offset);
		int s2_p1 = se2->ref_bg; int s2_p2 = s2_p1 + ((se2->is_ori)?0:se2->sv_info_p->end_offset);
		int isize = get_isize(s1_p1, s2_p1, se1->direction, se2->direction); if(isize > 0) return isize;
			isize = get_isize(s1_p1, s2_p2, se1->direction, se2->direction); if(isize > 0) return isize;
			isize = get_isize(s1_p2, s2_p1, se1->direction, se2->direction); if(isize > 0) return isize;
			isize = get_isize(s1_p2, s2_p2, se1->direction, se2->direction); if(isize > 0) return isize;
		return 0;
	}

	//it is a sub-function for "store_pair_end_score"
//	void read_is_proper_mated_Backup(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2, bool *read_is_nearby, bool *direction_right, int * isize){
//		if(se1 != NULL && se2 != NULL && se1->chrID == se2->chrID){
//			//when one read is ORI and its mapq is 0, read nearby is false
//			if((se1->is_ori && se1->mapq == 0) || (se2->is_ori && se2->mapq == 0)){
//				*read_is_nearby = false;
//			}
//			else{
//				//s1_p1 is the position of read 1 in first SV region, s1_p2 is the position in the second region
//				int s1_p1 = se1->ref_bg; int s1_p2 = s1_p1 + ((se1->is_ori)?0:se1->sv_info_p->end_offset);
//				int s2_p1 = se2->ref_bg; int s2_p2 = s2_p1 + ((se2->is_ori)?0:se2->sv_info_p->end_offset);
//				//first round check, for distance
//				if(		(ABS_U(s1_p1, s2_p1) < 2500) || //todo::///
//						(ABS_U(s1_p1, s2_p2) < 2500) ||
//						(ABS_U(s1_p2, s2_p1) < 2500) ||
//						(ABS_U(s1_p2, s2_p2) < 2500)){
//					*isize = (s2_p1 - s1_p1) + normal_read_len; //todo::~!~!
//					*read_is_nearby = true;
//				}
//				//second round check, for direction, the position of forward read must be smaller, the position of reverse must be bigger
//				if(*read_is_nearby && se1->direction != se2->direction){
//					int min_forward = (se1->direction == FORWARD)?(MIN(s1_p1, s1_p2)):(MIN(s2_p1, s2_p2));
//					int max_reverse = (se1->direction == REVERSE)?(MAX(s1_p1, s1_p2)):(MAX(s2_p1, s2_p2));
//					if(min_forward <= max_reverse)
//						*direction_right = true;
//				}
//			}
//		}
//	}

	inline bool one_read_is_new(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2){ return ((se1 != NULL && !se1->is_ori) || (se2 != NULL && !se2->is_ori) )?true:false; }

	//calculate the final score of single end results part se1 and se2
	int store_pair_end_score(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2){
		//detect whether reads is nearby
		int ISIZE =	read_is_proper_mated(se1, se2);
		//detect whether one of reads is new alignment
		int basic_score = ((se1 == NULL)?0:se1->align_score) + ((se2 == NULL)?0:se2->align_score); //basic scores
		int final_score = basic_score +
				+ ((ISIZE > 0)?(0):(-60)) //not nearby penalty
				+ (one_read_is_new(se1, se2)?(0):(+1)); // all reads are not new penalty

		//if(read_pair_full_match == false && final_score > FINAL_SCORE_CUTOFF && ISIZE > 0){ read_pair_full_match = true; }

		if(final_score >= max_score)
			set_score(final_score, se1, se2, ISIZE);
		return 0;
	}
};

#define K_T 22 		//kmer length used for building deBGA index

//data used for each ALIGN-PROCESSIONG thread
struct Classify_buff_pool
{
	//pair end handler:
	PE_score ps;
	//single end handler
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

	bam1_t 	*ori_b1;
	bam1_t 	*ori_b2;

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

	htsFile *output_file_ori = NULL;//results output file, used for storing reads that still not aligned well

};

struct deCOY_CLASSIFY_MAIN{
	CLASSIFY_SHARE_DATA *share = NULL;
	void init_run(int argc, char *argv[]);
private:
	static void *classify_pipeline(void *shared, int step, int tid, void *_data);									//pipeline
	static int load_reads(kstream_t *_fp1, kseq_t *_seqs1, kseq_t *_seqs2, int n_needed, MAP_PARA *o, Classify_buff_pool *cbp);	//function in step 1
	static void inline worker_for(void *_data, long data_index, int thread_index); 									//function in step 2
	static void output_results(int readNum, deBGA_INDEX * idx, htsFile *output_file, bam1_t *b1, bam1_t *b2, htsFile *output_file_ori, bam1_t *b1_ori, bam1_t *b2_ori);		//function in step 3
};

#endif /* READ_REALIGNMENT_HPP_ */
