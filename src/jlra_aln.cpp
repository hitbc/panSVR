/*
 * jlra_aln.cpp
 *
 *  Created on: 2021年3月10日
 *      Author: fenghe
 */

/*
 * classify_main.c
 *
 *  Created on: 2018-5-14
 *      Author: fenghe
 */

#include "cpp_lib/graph.hpp"
#include "jlra_aln.hpp"
#include <string.h>

#include <getopt.h>

#include <stdio.h>
#include <sys/time.h>

#include "zlib.h"

extern "C"{
#include "clib/kthread.h"
#include "clib/utils.h"
#include "clib/desc.h"
}

#include <math.h>

//**************************************class: CLASSIFY_MAIN********************************/
#define PIPELINE_T_NUM 3//one for reading; one for classifying; one for writing
#define STEP_NUM PIPELINE_T_NUM
#define N_NEEDED 2000000 //2M read pair per time

void CLASSIFY_MAIN::init_run(int argc, char *argv[]){

	//share data
	share = (CLASSIFY_SHARE_DATA*)xcalloc(1, sizeof(CLASSIFY_SHARE_DATA));

	//loading parameters
	share->o = (MAP_PARA *)xcalloc(1, sizeof(MAP_PARA));//mapping option
	if(share->o->get_option(argc, argv) != 0)
		return;
	//loading index
	share->idx = (deBGA_INDEX *)xcalloc(1, sizeof(deBGA_INDEX));
	//original header:
	fprintf(stderr, "Open original header file [%s]\n", share->o->ori_header_fn);
	htsFile *ori_header_file = hts_open(share->o->ori_header_fn, "r");//open output file
	bam_hdr_t *ori_header = sam_hdr_read(ori_header_file);
	hts_close(ori_header_file);
	share->idx->ori_header = ori_header;

	share->idx->load_index_file(share->o->indexDir);

	//get pipeline threads data
	share->data = (CLASSIFY_THREAD_DATA *)xcalloc(PIPELINE_T_NUM, sizeof(CLASSIFY_THREAD_DATA));
	for(int i = 0; i < PIPELINE_T_NUM; i++)
	{
		share->data[i].seqs1	= (kseq_t *)xcalloc(N_NEEDED, sizeof(kseq_t));
		share->data[i].seqs2	= (kseq_t *)xcalloc(N_NEEDED, sizeof(kseq_t));
		share->data[i].b1		= (bam1_t *)xcalloc(N_NEEDED, sizeof(bam1_t));
		share->data[i].b2		= (bam1_t *)xcalloc(N_NEEDED, sizeof(bam1_t));
		share->data[i].share_data_pointer = share;//anti-register for share data
	}

	//get alignment thread data
	share->buff = (Classify_buff_pool*)xcalloc(share->o->thread_n, sizeof(Classify_buff_pool));
	for(int i = 0; i < share->o->thread_n; i++){
		share->buff[i].SE_h[0].init(share->o, share->idx);
		share->buff[i].SE_h[1].init(share->o, share->idx);
	}

	//running alignment
	fprintf(stderr,"Start classify\n");
	double cpu_time = cputime();
	struct timeval start;
	gettimeofday(&start, NULL);
	//open SAM/BAM files
	//fprintf(stderr, "Open SAM/BAM file [%s]\n", share->o->read_fastq1);
	gzFile fp1 = xzopen(share->o->read_fastq1, "rb");
	//output bam file
	char bam_name[1024];
	if(share->o->sam_path == NULL){
		sprintf(bam_name, "%s.bam", share->o->read_fastq1);
		share->o->sam_path = bam_name;
	}
	if(share->o->output_sam)
		share->output_file =  hts_open(share->o->sam_path, "w");
	else
		share->output_file =  hts_open(share->o->sam_path, "wb");

	xassert(sam_hdr_write(share->output_file, share->idx->ori_header) >=0, "write header wrong!");//write header

	share->_fp1 = ks_init(fp1);
	fprintf(stderr, "Processing file: [%s].\n", share->o->read_fastq1);
	kt_pipeline(PIPELINE_T_NUM, classify_pipeline, share, STEP_NUM);
	//close FASTA files
	gzclose(fp1);

	hts_close(share->output_file);//close output file;
	fprintf(stderr, "Classify CPU: %.3f sec\n", cputime() - cpu_time);

}

#define MAX_read_size 100000000//100M

void *CLASSIFY_MAIN::classify_pipeline(void *shared, int step, int tid, void *_data) {
	CLASSIFY_SHARE_DATA * s = (CLASSIFY_SHARE_DATA*) shared;
	//step0: read read data from files; step1: process; step2: output result
	if 		(step == 0)	{ 	if((s->data[tid].readNum = load_reads(s->_fp1, s->data[tid].seqs1, s->data[tid].seqs2, N_NEEDED, s->o->max_use_read))) 	return (void *)1; }
	else if (step == 1)	{	kt_for(s->o->thread_n, worker_for, s->data + tid, s->data[tid].readNum);  							return (void *)1; }
	else if (step == 2)	{	output_results(s->data[tid].readNum, s->idx, s->output_file, s->data[tid].b1, s->data[tid].b2);    	return (void *)1; }
	return 0;
}

//function in step 1
static int load_read_number = 0;
int CLASSIFY_MAIN::load_reads(kstream_t *_fp1, kseq_t *_seqs1, kseq_t *_seqs2, int n_needed, int MAX_load_read_number)
{
	kseq_t *temp1 = _seqs1;
	kseq_t *temp2 = _seqs2;
	int i, rst1 = 0, rst2 = 0, total_length = 0;
	for( i = 0; i < n_needed &&  total_length < MAX_read_size && load_read_number < MAX_load_read_number; ++i, load_read_number++)
	{
		temp1[i].f = _fp1; //load two reads from same file
		temp2[i].f = _fp1;
		rst1 = kseq_read(temp1 + i);
		rst2 = kseq_read(temp2 + i);
		if(rst1 < 0 || rst2 < 0)
			break;
		total_length += (temp1[i].seq.l + temp2[i].seq.l);
	}
	return i;
}

extern void align_read_pair(kseq_t * read1,  kseq_t * read2, Classify_buff_pool * buff, bam1_t *b1, bam1_t *b2);
//function in step 2
void inline CLASSIFY_MAIN::worker_for(void *_data, long data_index, int thread_index)
{ // kt_for() callback
	CLASSIFY_THREAD_DATA *d = (CLASSIFY_THREAD_DATA *) _data;
	CLASSIFY_SHARE_DATA  *s = (CLASSIFY_SHARE_DATA  *) (d->share_data_pointer);
	align_read_pair(d->seqs1 + data_index, d->seqs2 + data_index, s->buff + thread_index, d->b1 + data_index, d->b2 + data_index);
}

//GLOBAL VALUABLE

static int read_block_count = 0;
void CLASSIFY_MAIN::output_results(int readNum, deBGA_INDEX * idx, htsFile *output_file, bam1_t *b1, bam1_t *b2){
	fprintf(stderr, "Processing %d reads, at block ID %d\n", readNum, read_block_count++);
	for(int i = 0; i < readNum; i++){
		if(b1[i].core.tid != -1)	xassert(sam_write1(output_file, idx->ori_header, b1 + i) >= 0, "");
		if(b2[i].core.tid != -1)	xassert(sam_write1(output_file, idx->ori_header, b2 + i) >= 0, "");
	}
}

//**************************************class: single_end_handler********************************/

uint8_t charToDna5n[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,//4 the last but one
    /*    A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//'Z'
    /*             T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*    a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*             t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

uint64_t inline getKmer(uint32_t read_off, uint64_t *read_bit){
	uint32_t index_word = read_off >> 5;
	uint32_t index_in_word = (read_off & 0X1f);
	uint64_t full_kmer = (((read_bit[index_word]) << ((index_in_word) << 1)) | ((index_in_word == 0)?0:(read_bit[index_word+ 1] >> ((32 - index_in_word) << 1))));
	return full_kmer >> (LEN_KMER_LEFT << 1);
	//return (((read_bit[index_word] & kmerMask[index_in_word]) << ((index_in_word - LEN_KMER_LEFT) << 1)) | (read_bit[index_word+ 1] >> ((64 - LEN_KMER - index_in_word) << 1)));
}

static int sort_output(Graph_handler &g, deBGA_INDEX * idx, MAX_IDX_OUTPUT &rst, int direction){
	//output the best
	if(g.vertexArr_size == 0) return 0;

	//get max score
	g.max_index = MAX_uint32_t; g.max_distance = 0;
	g.same_top_max_distance_id_list.clear();
	g.same_top_max_distance_id_list.emplace_back(g.max_index);
	int i = g.vertexArr_size - 1;
   	for(; i >= 0; i--){
   		if(g.dist_path[i].already_used)
   			continue;
   		float c_dist = g.dist_path[i].dist;
		if (g.max_distance < c_dist){
			g.max_distance = c_dist;
			g.max_index = i;
			g.same_top_max_distance_id_list.clear();
			g.same_top_max_distance_id_list.emplace_back(i);
		}
		else if(g.max_distance == c_dist){
			g.same_top_max_distance_id_list.emplace_back(i);
		}
   	}

//   	if(g.max_distance == 127){
//   		printf(" ");
//   	}
   	if(g.max_index == MAX_uint32_t)//without reult
   		return 0;

   	int node_num_already_used = 0;
   	int node_num_not_used = 0;
   	//get chain for that score
   	uint32_t same_size = g.same_top_max_distance_id_list.size();
   	if(same_size > 1){ //random select a new chain
   		g.max_index = g.same_top_max_distance_id_list[rand() % same_size];
   	}

	int first_node = g.max_index;
	int original_first_node = first_node;

	//debug code::
	//fprintf(stderr, "\n\ng.max_index: [%d], g.max_distance:[%f]\n", g.max_index, g.max_distance);
	for(; first_node != -1; ){
		//debug code::
		//g.vertexArr_[0][first_node].print();
		if(g.dist_path[first_node].already_used == true)
			node_num_already_used++;
		else
			node_num_not_used++;
		g.dist_path[first_node].already_used = true;
		int next_node = g.dist_path[first_node].pre_node;
		if(next_node == -1)
			break;
		first_node = next_node;
	}
	int original_final_node = first_node;

	//total node number within region is far bigger then all used node in a chain: treat as STR or VNTR
	if(original_first_node - original_final_node > ((node_num_not_used + node_num_already_used + 5) << 1)){
		for(int node_ID = original_final_node; node_ID < original_first_node; node_ID ++)
			g.dist_path[node_ID].already_used = true;
	}

	if(node_num_already_used >= node_num_not_used){
		return sort_output(g, idx, rst, direction);
	}

	//debug code::
	//fprintf(stderr, "read_begin: [%d], ref_begin:[%d]\n", g.vertexArr_[0][first_node].read_begin, g.vertexArr_[0][first_node].ref_begin);
	//set result
	int ref_begin = g.vertexArr_[0][first_node].ref_begin;
	int chr_ID = idx->get_chromosome_ID(ref_begin);

	rst.direction = direction;
	rst.max_index = g.max_index;
	rst.chain_score = g.max_distance;
   	rst.read_bg = g.vertexArr_[0][first_node].read_begin;
   	rst.chrID = chr_ID;
  	rst.ref_bg = ref_begin - idx->chr_end_n[chr_ID - 1];//end of last chr is the begin of this chr
	return 1;
}

static void binary_read_64_bit(int read_len, uint8_t * read_seq, uint64_t* read_bit64_FOR){
	unsigned long int read_bit_char = (((uint16_t )((read_len >> 5) + 1)) << 3);
	memset(read_bit64_FOR, 0, read_bit_char);
	for (int r_i = 0;r_i < read_len; r_i++)
		read_bit64_FOR[r_i >> 5] |= (((uint64_t )read_seq[r_i]) << ((31 - (r_i & 0X1f)) << 1));
}

#define MAX_BIN_READ_LEN_64 50 //1600/32 = 50
#define MIN_STR_REPEAT_COUNT 4 //when a kmer count >= MIN_STR_REPEAT_COUNT, that kmer will be treated as an STR kmer
#define MIN_STR_DETECT_LEN 15 //when a read has more than 15 STR-kmer, it will be treat as an STR read

//#define KMER_COUNT_ANALYSIS 1

int get_ksw_score(Graph_handler &g, int first_node, int read_l, KSW_ALN_handler &kswh){

	std::vector<PATH_t> &dist_path_p = g.dist_path;
	std::vector<UNI_SEED>& vertexArr_p = g.vertexArr_[0];

	//get ksw score
	int aln_read_begin = read_l;
	int aln_read_end = read_l;
	int aln_ref_begin = MAX_int32t;
	int aln_ref_end = MAX_int32t;
	if(kswh.read_score > 1000){	printf(" ");}

	int last_aln_begin = read_l;
	int last_ref_begin = MAX_int32t;
	int UNITIG_MIS = 0;
	for(; first_node != -1; ){
		//debug code::
		//vertexArr_p[first_node].print();
		int MEM_read_beg = vertexArr_p[first_node].read_begin;
		int MEM_read_end = vertexArr_p[first_node].read_end;
		int MEM_ref_beg = vertexArr_p[first_node].ref_begin;
		int MEM_ref_end = vertexArr_p[first_node].ref_end;

		aln_read_begin = MIN(aln_read_begin, MEM_read_end);
		aln_ref_begin  = MIN(aln_ref_begin,  MEM_ref_end);
		if(aln_read_begin <= aln_read_end){
			//fprintf(stderr, "##%d %d %d %d \n", aln_read_begin, aln_read_end, aln_ref_begin, aln_ref_end);
			//store cigar:
			if(aln_read_end < last_aln_begin){
				int MEM_LEN = last_aln_begin - aln_read_end;
				//kswh.alignment(aln_read_end, aln_read_end + MEM_LEN, last_ref_begin, last_ref_begin + MEM_LEN, KSW_ALN_end_to_end);
				UNITIG_MIS += kswh.get_misMatch(aln_read_end, aln_read_end + MEM_LEN, last_ref_begin, last_ref_begin + MEM_LEN);
				kswh.cigar_tmp.emplace_back('M', (uint16_t)(MEM_LEN));
			}
			last_aln_begin = aln_read_begin;
			if(aln_ref_end == MAX_int32t){
				aln_ref_end = aln_ref_begin + (aln_read_end - aln_read_begin) + 30;
				kswh.alignment(aln_read_begin, aln_read_end, aln_ref_begin, aln_ref_end, KSW_ALN_right_extend);
			}
			else
				kswh.alignment(aln_read_begin, aln_read_end, aln_ref_begin, aln_ref_end, KSW_ALN_end_to_end);
		}
		else{
			int distance_read = aln_read_end - aln_read_begin;
			int distance_ref =  aln_ref_end - aln_ref_begin;
			if(distance_read != distance_ref){
				int deletion_len = distance_ref - distance_read;
				//kswh.cigar_tmp.emplace_back('D', deletion_len);
				int ABS_deletion_len = ABS(deletion_len);
				int score_1 = kswh.gap_open_D + (ABS_deletion_len - 1) * kswh.gap_ex_D;
				int score_2 = kswh.gap_open2_D + (ABS_deletion_len - 1) * kswh.gap_ex2_D;
				int score = MIN(score_1, score_2);
				kswh.read_score -= score;
			}
		}
		//else fprintf(stderr, "@@%d %d \n", MEM_read_beg, aln_read_begin);
		aln_read_end = MEM_read_beg;
		last_ref_begin = MEM_ref_beg;
		aln_ref_end  = MEM_ref_beg;
		int next_node = dist_path_p[first_node].pre_node;
		if(next_node == -1)
			break;
		first_node = next_node;
	}

	if(aln_read_end < last_aln_begin){
		int MEM_LEN = last_aln_begin - aln_read_end;
		//kswh.alignment(last_aln_begin - last_aln_begin, last_aln_begin, last_ref_begin, last_ref_begin + MEM_LEN, KSW_ALN_end_to_end);
		UNITIG_MIS += kswh.get_misMatch(aln_read_end, aln_read_end + MEM_LEN, last_ref_begin, last_ref_begin + MEM_LEN);
		kswh.cigar_tmp.emplace_back('M', MEM_LEN);
	}
	//kswh.cigar_tmp.emplace_back('M', last_aln_begin - aln_read_end);

	aln_read_begin = 0;
	aln_ref_begin = 0;
	int read_begin_alignment = 0;
	if(aln_read_begin < aln_read_end){
		aln_ref_begin = aln_ref_end - (aln_read_end - aln_read_begin) - 30;
		aln_ref_begin = MAX(0, aln_ref_begin);
		//fprintf(stderr, "##%d %d %d %d \n", aln_read_begin, aln_read_end, aln_ref_begin, aln_ref_end);
		kswh.alignment(aln_read_begin, aln_read_end, aln_ref_begin, aln_ref_end, KSW_ALN_left_extend);
		//reset begin position:
		if(aln_ref_end > aln_ref_begin){
			if(kswh.is_simple_aln)
				read_begin_alignment = aln_ref_end - aln_ref_begin - 30;//only match
			else
				read_begin_alignment = aln_ref_end - aln_ref_begin;//with deletion or insertion
		}
	}
	kswh.read_score += (read_l - kswh.total_q_len)*kswh.match_D;
	kswh.read_score -= (UNITIG_MIS)*(kswh.match_D +kswh.mismatch_D);
	return read_begin_alignment;
}

#define MAX_CHAIN_SOCRE_DIFF 30 //when chain score results are less the [best chain score - MIN_CHAIN_SOCRE_DIFF], the search will be break;
#define MIN_CHAIN_SOCRE 30 //when chain score results are less the [best chain score - MIN_CHAIN_SOCRE_DIFF], the search will be break;
#define MIN_ALN_SOCRE 40

void single_end_handler::align(){

	//initiation alignment flags
	result_num = 0;
	primary_result = NULL;
	secondary_result = NULL;
	parse_ori_mapping_rst(c_read->comment.s, read_l);
	if(ori.chrID > 24) ORI_is_UNMAPPED = true;
	if(!ORI_is_UNMAPPED && ori.align_score == read_l * kswh.match_D) return; //refuse to align already full score read

	binary_read_2_bit();
	for(int ori_i = 0; ori_i < 2; ori_i++)
		chainning_one_read(ori_i);
	uint32_t max_chain_score = 0;
	for(int ori_i = 0; ori_i < 2; ori_i++){// for each direction
		int direction = (ori_i == 0)?FORWARD:REVERSE;
		for(int i = 0; i < MAX_OUTPUT_NUMBER; i++){ // for each results
			int rst = sort_output(g[ori_i], idx, result[result_num], direction);
			if(rst == 0) break; // stop search reason 1: no result
			else{
				uint32_t c_chain_score = result[result_num].chain_score;
				max_chain_score = MAX(c_chain_score, max_chain_score);
				//stop search reason 2: all other chain has low score
				if(c_chain_score + MAX_CHAIN_SOCRE_DIFF < max_chain_score || c_chain_score < MIN_CHAIN_SOCRE )	break;
				result_num ++;
			}
		}
	}
	//sort all result by score
	qsort(result, result_num, sizeof(MAX_IDX_OUTPUT), MAX_IDX_OUTPUT::cmp_chain_score);

	//uint32_t max_chain_score2 = (result_num[read_id] == 1)?0:result[read_id][1].chain_score;
	if(result_num == 0 || max_chain_score < MIN_CHAIN_SCORE)
		return;
	for(int rst_ID = 0; rst_ID < result_num; rst_ID++){
		MAX_IDX_OUTPUT & c_r = result[rst_ID];
		if(c_r.chain_score + MAX_CHAIN_SOCRE_DIFF < max_chain_score) { result_num = rst_ID; break; }
		int is_reverse = (c_r.direction == REVERSE)?true:false;
		kswh.setRead(bin_read[is_reverse]);
		int read_begin_alignment = get_ksw_score(g[is_reverse], c_r.max_index, read_l, kswh);
		//xassert(read_begin_alignment < (int)c_r.ref_bg, "");
		c_r.ref_bg -= read_begin_alignment;
		c_r.align_score = MAX(kswh.read_score, 0);
		if(c_r.reverseGIGAR(kswh.cigar_tmp, read_l) == false){
			fprintf(stderr, " %s\n", read_seq_char);
		}
	}
	qsort(result, result_num, sizeof(MAX_IDX_OUTPUT), MAX_IDX_OUTPUT::cmp_align_score);
	//store results
	uint32_t max_aln_score  = result[0].align_score;
	if(max_aln_score < MIN_ALN_SOCRE){	result_num = 0;	return;	}

	//set SV id and chr_id in original reference
	for(int i = 0; i < result_num; i++){
		uint32_t sv_ID = result[i].chrID;
		result[i].sv_info_p = &(idx->sv_info[sv_ID]);
		result[i].chrID = result[i].sv_info_p->chr_ID;
		result[i].ref_bg += result[i].sv_info_p->st_pos;
		if(result[i].ref_bg >= MAX_int32t)
			result[i].ref_bg = 5;
		result[i].is_ori = false;
		result[i].rst_idx = i;
		result[i].mapq = 0;
		result[i].has_mate = false;
	}
	if(result_num > 0){
		int32_t pri_minus_sec_score  = result[0].align_score - ((result_num > 1)?result[1].align_score:0);
		result[0].mapq = (pri_minus_sec_score > 40)?40:(pri_minus_sec_score);
	}
}

//set tid to be -1 to discard the output of BAM
void single_end_handler::output_BAM(bam1_t *b, bool is_first_read, bool read_proper_mate, char * comment_str)
{	//get primary alignment and mapQ
	//skip original alignment on decoy reference:
	//init:
	b->core.tid = -1;
	if(primary_result == NULL) return;//output nothing
	if(NOT_output_original_rst && primary_result->is_ori) return;

	int direction = primary_result->direction;

	uint8_t flag = ((is_first_read)?BAM_FIRST_READ:0) + ((direction == REVERSE)?BAM_STRAND:0) + ((primary_result->has_mate)?0:BAM_MATE_UNMAPPED);

	uint32_t read_l = c_read->seq.l;
	if(direction == REVERSE){
		getReverseStr_char(c_read->seq.s, read_l);
		getReverseStr_qual_char(c_read->qual.s, read_l);
	}
	s.l = 0;
	//BASIC PART:
	sprintf(s.s + s.l,	"%s\t%d\t%s\t%d\t%d\t", c_read->name.s, flag, idx->ori_header->target_name[primary_result->chrID], primary_result->ref_bg , primary_result->mapq); s.l += strlen(s.s + s.l);//read name//flag//chr_name //mapq, set to blank at first//ref begin
	//CIGAR
	std::vector<CIGAR_PATH> &cigar = primary_result->cigar;
	for(auto & c : cigar){ sprintf(s.s + s.l, "%d%c", c.size, OUT_BAM_CIGAR_STR[c.type]); s.l += strlen(s.s + s.l);	}
	sprintf(s.s + s.l, "\t"	); s.l ++;
	// mate pos and ISIZE: isize is simple 500 for mated read
	int isize = (read_proper_mate)?500:0;
	if(primary_result->has_mate){ sprintf(s.s + s.l, "%s\t%d\t%d\t", idx->ori_header->target_name[primary_result->mate_chrID], primary_result->mate_ref_bg, isize);} else sprintf(s.s + s.l,	"*\t0\t0\t"); s.l += strlen(s.s + s.l);
	//seq//qual //align score
	sprintf(s.s + s.l,  "%s\t%s\t" "AS:i:%d\t", c_read->seq.s, c_read->qual.s, primary_result->align_score); s.l += strlen(s.s + s.l);
	//ORIGINAL SCORE and original alignment
	//original score, alignment and position, "OA:Z:000782F,-60173,250M,8;"
	sprintf(s.s + s.l, "OS:i:%d\t" "OA:Z:%d,%d,%d,%d,%c;\t", ori.align_score, ori.chrID, ori.ref_bg, ori.read_bg, ori.mapq, ORI_is_UNMAPPED?'U':'M'); s.l += strlen(s.s + s.l);
	//chaining score ,when it is new alignment results
	if(!primary_result->is_ori){ sprintf(s.s + s.l, "CS:i:%d\t", primary_result->chain_score); s.l += strlen(s.s + s.l); }
	//SV_region_information, when it is not null
	if(primary_result->sv_info_p != NULL){ sprintf(s.s + s.l, "SV:Z:%s\t", primary_result->sv_info_p->vcf_print_string.c_str());s.l += strlen(s.s + s.l);}
	//mate SV info, when it is not null
	if(primary_result->mate_sv_info_p != NULL){	sprintf(s.s + s.l, "MV:Z:%s\t", primary_result->mate_sv_info_p->vcf_print_string.c_str()); s.l += strlen(s.s + s.l);	}
	//secondary alignment
	if(secondary_result != NULL){
		sprintf(s.s + s.l, "XA:Z:%d,%d,%d,"	"%d,%c,%s;\t",
				secondary_result->chrID, secondary_result->ref_bg, secondary_result->read_bg,
				secondary_result->align_score, (secondary_result->direction == FORWARD)?'F':'R', (secondary_result->sv_info_p == NULL)?"*":secondary_result->sv_info_p->vcf_id.c_str());
		s.l += strlen(s.s + s.l);
	}

	{sprintf(s.s + s.l, "RC:Z:%s\t", comment_str); s.l += strlen(s.s + s.l);}
	//output result:
	int sam_parse1_rst = sam_parse1(&s, idx->ori_header, b);
	if(sam_parse1_rst != 0)	fprintf(stderr, "@sam_parse1 ERROR\n");
//	b->core.tid = primary_result->chrID;
//	if(primary_result->has_mate)
//		b->core.mtid = primary_result->mate_chrID;
}

void single_end_handler::chainning_one_read(int read_is_reverse){
	//step 1: binary reads
	uint64_t read_bit_64[MAX_BIN_READ_LEN_64];
	//uint64_t read_bit_REV[MAX_BIN_READ_LEN_64];

	binary_read_64_bit(read_l, ((read_is_reverse)?(bin_read[1]):(bin_read[0])), read_bit_64);
	vertexm_v.clear();
	vertexu_v.clear();
	uniseed_v[read_is_reverse].clear();


	uint32_t kmer_number = read_l - LEN_KMER + 1;
	//test get repeat kmer, used for seed step 5.
	//seed list store whether a kmer can be selected as seed, 0:NO; >0: YES
	uint8_t * seed_list = repeat_seed_info;
	if(!read_is_reverse){// forward
		readIsSTR = false;
		std::map<uint64_t, int> &kmer_set = kmer_set_buff;
		kmer_set.clear();

		for(uint32_t read_off = 0; read_off < kmer_number; read_off ++) //every seed_l[tid] we choose a seed
		{
			uint64_t kmer_bit = getKmer(read_off, read_bit_64);
			std::map<uint64_t, int>::iterator it = kmer_set.find(kmer_bit);
			if(it!=kmer_set.end()){		it->second ++;			}
			else{						kmer_set[kmer_bit] = 1;	}
		}
		if(kmer_set.size() < kmer_number - MIN_STR_DETECT_LEN){//with at lest 15 repeat kmer
			readIsSTR = true;
			//mask unique kmer as 1
			for(uint32_t read_off = 0; read_off < read_l - LEN_KMER + 1; read_off ++) //every seed_l[tid] we choose a seed
			{
				uint64_t kmer_bit = getKmer(read_off, read_bit_64);
				std::map<uint64_t, int>::iterator it = kmer_set.find(kmer_bit);
				seed_list[read_off] = (it->second >= MIN_STR_REPEAT_COUNT)?0:1;
			}
			int STR_region_kmer_num_bg = 0;
			int STR_region_kmer_num_ed = 0;
			//set 5 kmer in the begin of read as 2 and 5 kmer at end of read as 4
			for(uint32_t read_off = 0; read_off < SEED_STEP; read_off ++) //every seed_l[tid] we choose a seed
			{
				STR_region_kmer_num_bg += (seed_list[read_off] == 0)?1:0;
				STR_region_kmer_num_ed += (seed_list[read_l - LEN_KMER - read_off] == 0)?1:0;
				seed_list[read_off] += 2;
				seed_list[read_l - LEN_KMER - read_off] += 4;
			}
			//if in the begin or end there is no STR kmer, select other 5 STR kmer to be used as seed
			if(STR_region_kmer_num_bg < SEED_STEP && STR_region_kmer_num_ed < SEED_STEP){
				int total_VNTR_KMER_numer = 0;
				for(uint32_t read_off = 0; total_VNTR_KMER_numer < SEED_STEP && read_off < read_l - LEN_KMER + 1; read_off ++) //every seed_l[tid] we choose a seed
				{
					if(seed_list[read_off] > 0) continue;
					seed_list[read_off] += 8;
					total_VNTR_KMER_numer++;
				}
			}
			//debug code:
			//for(uint32_t read_off = 0; read_off < read_l - LEN_KMER + 1; read_off ++)	{uint64_t kmer_bit = getKmer(read_off, read_bit_64); std::map<uint64_t, int>::iterator it = kmer_set.find(kmer_bit); fprintf(stderr, "[%d %d]\t", it->second, seed_list[read_off]);	}
			//fprintf(stdout, "%ld\t", kmer_set.size());	for(uint32_t i = 0; i < read_l; i++){	fprintf(stderr, "%c", "ACGT"[read_seq[i]]);	}fprintf(stderr, "\n");
		}
	}
	else{
		if(readIsSTR)
			getReverseStr_qual(seed_list, read_l - LEN_KMER + 1);//reverse the seed list
	}
	g[read_is_reverse].readIsSTR = readIsSTR;
#ifdef KMER_COUNT_ANALYSIS

if( readIsSTR && !read_is_reverse){
	fprintf(stdout, "%ld\t", kmer_set_buff.size());
	for(uint32_t i = 0; i < read_l; i++){	fprintf(stdout, "%c", read_seq_char);	}
	fprintf(stdout, "\n");
}
return;
#endif

	uint32_t max_search_right = 0;
	for(uint32_t read_off = 0; read_off < kmer_number; read_off += SEED_STEP) //every seed_l[tid] we choose a seed
	{
		if(read_off + LEN_KMER - 1 <= max_search_right) continue;
		if(readIsSTR && seed_list[read_off] == 0) continue;
		uint64_t kmer_bit = getKmer(read_off, read_bit_64);
		int64_t range[2];
		bool search_rst = idx->search_kmer(LEN_KMER, kmer_bit, range, SEED_OFFSET);//search kmer in the index
		if(!search_rst) continue;

		uint64_t hit_bg = range[0];
		uint64_t hit_ed = range[1];

		if ((hit_ed - hit_bg + 1) > UNI_POS_N_MAX)
			continue;
		//for each kmer in the index:
		uint32_t max_right_i = 1;
		for (uint64_t hit_i = hit_bg; hit_i <= hit_ed; ++hit_i)
			idx->UNITIG_MEM_search(hit_i, vertexm_v, read_bit_64, read_off, read_l, LEN_KMER, max_right_i);
		//skip seed when MEM covered
		max_search_right = read_off + LEN_KMER + max_right_i - 1;
	}

	//merge seed in the unipath:
	idx->merge_seed_in_unipath(vertexm_v, vertexu_v);
	//expand seed to reference
	idx->expand_seed(vertexu_v, uniseed_v[read_is_reverse], &rand_buff);
	//chaining
	//chaining step 1: create graph and SDP
	g[read_is_reverse].process(uniseed_v[read_is_reverse]);
}

void single_end_handler::binary_read_2_bit(){
	for (int r_i = 0; read_seq_char[r_i]; r_i++){
		char tmp_char = read_seq_char[r_i];
		if (tmp_char == 'N') tmp_char = "ACGT"[rand()%4]; //random
		uint8_t c_tmp = charToDna5n[(uint8_t)tmp_char];
		bin_read[0][r_i] = c_tmp;
		bin_read[1][(read_l - r_i - 1)] = (c_tmp ^ 0X3);
	}
}

struct PE_score{
	int max_same;
	int max_score;
	bool proper_mate;
	MAX_IDX_OUTPUT *max_1;
	MAX_IDX_OUTPUT *max_2;


	void init(){
		max_same = 1;
		max_score = 0;
		max_1 = NULL;
		max_2 = NULL;
		proper_mate = false;
	}

	void set_score(int new_score, MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2, bool proper_mate_){
		if(new_score > max_score){
			max_1 = se1; max_2 = se2; max_score = new_score; proper_mate = proper_mate_;
			max_same = 1;
		}
		else if(new_score == max_score){
			max_same ++;
			if(rand() % max_same == 0){
				max_1 = se1; max_2 = se2; max_score = new_score; proper_mate = proper_mate_;
			}
		}
	}

	int store_pair_end_score(MAX_IDX_OUTPUT *se1, MAX_IDX_OUTPUT *se2){
		//detect whether reads is nearby
		bool read_is_nearby = false;//when two reads are not nearby, score -= 60
		bool direction_right = false;//when two reads are not nearby or direction not right, score -= 20
		if(se1 != NULL && se2 != NULL && se1->chrID == se2->chrID){
			//when one read is ORI and its mapq is 0, read nearby is false
			if((se1->is_ori && se1->mapq == 0) || (se2->is_ori && se2->mapq == 0)){
				read_is_nearby = false;
			}
			else{
				//s1_p1 is the position of read 1 in first SV region, s1_p2 is the position in the second region
				int s1_p1 = se1->ref_bg; int s1_p2 = s1_p1 + ((se1->is_ori)?0:se1->sv_info_p->end_offset);
				int s2_p1 = se2->ref_bg; int s2_p2 = s2_p1 + ((se2->is_ori)?0:se2->sv_info_p->end_offset);
				//first round check, for distance
				if(		(ABS_U(s1_p1, s2_p1) < 2500) ||
						(ABS_U(s1_p1, s2_p2) < 2500) ||
						(ABS_U(s1_p2, s2_p1) < 2500) ||
						(ABS_U(s1_p2, s2_p2) < 2500)){
					read_is_nearby = true;
				}
				//second round check, for direction, the position of forward read must be smaller, the position of reverse must be bigger
				if(read_is_nearby && se1->direction != se2->direction){
					int min_forward = (se1->direction == FORWARD)?(MIN(s1_p1, s1_p2)):(MIN(s2_p1, s2_p2));
					int max_reverse = (se1->direction == REVERSE)?(MAX(s1_p1, s1_p2)):(MAX(s2_p1, s2_p2));
					if(min_forward <= max_reverse)
						direction_right = true;
				}
			}
		}
		//detect whether one of reads is new alignment
		bool one_read_is_new = false;
		if((se1 != NULL && !se1->is_ori) || (se2 != NULL && !se2->is_ori) )
			one_read_is_new = true;
		int score = ((se1 == NULL)?0:se1->align_score) + ((se2 == NULL)?0:se2->align_score) + (read_is_nearby?(0):(-60)) + (direction_right?(0):(-20)) + (one_read_is_new?(0):(-10));
		if(score >= max_score)
			set_score(score, se1, se2, read_is_nearby);
		return 0;
	}
};

void align_read_pair(kseq_t * read1,  kseq_t * read2, Classify_buff_pool * buff, bam1_t *b1, bam1_t *b2){
	single_end_handler *SE_h = &(buff->SE_h[0]);
	//step 1: align both read separately
	for(int read_id = 0; read_id < 2; read_id++){
		single_end_handler &c_SE_h = SE_h[read_id];
		SE_h[read_id].read_register((read_id == 0)?read1:read2);
		c_SE_h.align();
	}
	//step 2: re-select alignment results by pair end information
	//all results, including the original result will be used to get pair end score
	//all results stored in a list, first new results, then ori_results
	//PE score will be stored in a matrix
	PE_score ps;
	ps.init();

	//for both original
	int se0_rst_num = SE_h[0].result_num; int se1_rst_num = SE_h[1].result_num;
	if(!SE_h[0].ORI_is_UNMAPPED) se0_rst_num++;
	if(!SE_h[1].ORI_is_UNMAPPED) se1_rst_num++;

	//if(se0_rst_num == 0 )
	//when using only second end
	for(int i = 0; i < se0_rst_num; i++)
		ps.store_pair_end_score(&((i < SE_h[0].result_num)?SE_h[0].result[i]:SE_h[0].ori), NULL);
	//when using only first end
	for(int j = 0; j < se1_rst_num; j++)
		ps.store_pair_end_score(NULL, &((j < SE_h[1].result_num)?SE_h[1].result[j]:SE_h[1].ori));
	//when using both end
	for(int i = 0; i < se0_rst_num; i++){
		MAX_IDX_OUTPUT & se1 = (i < SE_h[0].result_num)?SE_h[0].result[i]:SE_h[0].ori;
		for(int j = 0; j < se1_rst_num; j++)
			ps.store_pair_end_score(&se1, &((j < SE_h[1].result_num)?SE_h[1].result[j]:SE_h[1].ori));
	}
	//discard result when both best single end results are original(or NULL) OR max score is 0
	if(ps.max_score <= 0 || ((ps.max_1 == NULL || ps.max_1->is_ori) && (ps.max_2 == NULL || ps.max_2->is_ori))){
		ps.max_1 = NULL;
		ps.max_2 = NULL;
	}
	for(int i = 0; i < 2; i++){
		MAX_IDX_OUTPUT *c_rst = (i == 0)?ps.max_1:ps.max_2;
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
		MAX_IDX_OUTPUT *mate_rst = (i == 0)?ps.max_2:ps.max_1;
		if(mate_rst != NULL){
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
			c_rst->mate_sv_info_p = NULL;
		}
	}

	//step3: output result by pre-defined 'PRIMARY' and 'SECONDARY' tag, stored in bam format
	for(int read_id = 0; read_id < 2; read_id++){
		single_end_handler &c_SE_h = SE_h[read_id];
		c_SE_h.output_BAM((read_id == 0)?b1:b2, 1 - read_id, ps.proper_mate, (read_id == 0)?read1->comment.s:read2->comment.s);
	}
}

void getReverseStr_uint8_t(uint8_t * q, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half; i++){
		int ri = len - 1 - i;
		uint8_t tmp = q[i ];
		q[i ] = q[ri];
		q[ri] = tmp;
	}
}

//**************************************class: KSW_ALN_handler********************************/

void KSW_ALN_handler::copy_option(MAP_PARA 	*o){
	match_D= o->match_D ;
	mismatch_D= o->mismatch_D ;
	gap_open_D= o->gap_open_D ;
	gap_ex_D= o->gap_ex_D ;
	gap_open2_D= o->gap_open2_D ;
	gap_ex2_D= o->gap_ex2_D ;
	zdrop_D= o->zdrop_D ; //for DNA zdrop = 400, 200 for RNA
	bandwith = 200;
	flag = 0;
}

void KSW_ALN_handler::ksw_gen_mat_D(){
	int8_t l,k,m;
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m) {
			mata_D[k] = l == m ? match_D : -(mismatch_D);	/* weight_match : -weight_mismatch */
			k++;
		}
		mata_D[k] = 0; // ambiguous base
		k++;
	}
	for (m = 0; m < 5; ++m) {
		mata_D[k] = 0;
		k++;
	}
}

void KSW_ALN_handler::init(MAP_PARA 	*o, deBGA_INDEX * idx_){
	tseq = (uint8_t *)xmalloc(1600);
	qseq_rev = (uint8_t *)xmalloc(1600);
	km = km_init();
	memset(&ez, 0, sizeof(ksw_extz_t));
	//mapping options
	copy_option(o);
	ksw_gen_mat_D();
	//index
	idx = idx_;
}

void KSW_ALN_handler::free_memory(){
	if(ez.cigar)   kfree(km, (void *)(ez.cigar));
	km_destroy(km);
	free(tseq);
	free(qseq_rev);
}

void KSW_ALN_handler::setRead(uint8_t *read_str_){
	cigar_tmp.clear();
	read_str = read_str_;
	//n_cigar = 0;
	read_score = 0;
	total_q_len = 0;
}

void KSW_ALN_handler::align_non_splice()
{
	if ((int64_t)tlen * qlen > 1000000)
	{
		ksw_reset_extz(&ez);
		if(ez.m_cigar < 2)
		{
			ez.n_cigar = 2;
			ez.m_cigar = (ez.n_cigar)<<2;
			ez.cigar = (uint32_t *)krealloc(km, (void *)ez.cigar, (ez.m_cigar)<<4);
		}
		ez.n_cigar = 2;
		ez.cigar[0] = qlen<<4 | 1;
		ez.cigar[1] = tlen<<4 | 3;
		ez.score = 0;
	}
	else{
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}
}

int KSW_ALN_handler::get_misMatch(int read_st_, int read_ed_, int ref_st_, int ref_ed_){
	qlen = read_ed_ - read_st_;
	qseq = read_str + read_st_;
	tlen = ref_ed_ - ref_st_;
	if(ref_ed_ < ref_st_){//when ref_len < 0
		tlen = 0;
		qlen += (ref_st_ - ref_ed_);
	}

	xassert(tlen < 1600, "");
	idx->get_refseq(tseq, tlen, ref_st_);
	int c_simple_NM = 0;
	for(uint32_t i = 0; i < qlen; c_simple_NM += (qseq[i] != tseq[i])?1:0, i++);
	if(c_simple_NM > 3)	c_simple_NM = 3;
	return c_simple_NM;
}

void KSW_ALN_handler::alignment(int read_st_, int read_ed_, int ref_st_, int ref_ed_, int type_){
	qlen = read_ed_ - read_st_;
	qseq = read_str + read_st_;
	tlen = ref_ed_ - ref_st_;
	if(ref_ed_ < ref_st_){//when ref_len < 0
		tlen = 0;
		qlen += (ref_st_ - ref_ed_);
	}

	xassert(tlen < 1600, "");
	idx->get_refseq(tseq, tlen, ref_st_);
	type = type_;

	if(type == KSW_ALN_left_extend){
		getReverseStr_uint8_t(tseq, tlen);
		memcpy(qseq_rev, qseq, qlen);
		getReverseStr_uint8_t(qseq_rev, qlen);
		qseq = qseq_rev;
	}

	total_q_len += qlen;

	//try simple compare
	is_simple_aln = false;
	simple_NM = 0;
	if(qlen == 0 || tlen == 0){
		is_simple_aln = true;
		simple_NM = qlen + tlen;
	}else if(qlen == tlen || type != KSW_ALN_end_to_end){//when [global search and length is equal] or [extension search]
		for(uint32_t i = 0; i < qlen && simple_NM < 6; simple_NM += (qseq[i] != tseq[i])?1:0, i++);
		if(simple_NM == 1 || (simple_NM < 6 && ((simple_NM << 3) < qlen)))
			is_simple_aln = true;
	}

	if(!is_simple_aln)
		 align_non_splice();

	//store score and cigar:

	//set read score
	if(is_simple_aln){
		if(qlen == 0 || tlen == 0){//deletion, insertion
			if(simple_NM != 0){
				int score_1 = gap_open_D + (simple_NM - 1) * gap_ex_D;
				int score_2 = gap_open2_D + (simple_NM - 1) * gap_ex2_D;
				int score = MIN(score_1, score_2);
				read_score -= score;
			}
		}
		else
			read_score += qlen * match_D - simple_NM * (match_D + mismatch_D);
		//cigar:
		if(qlen == 0)
			cigar_tmp.emplace_back('D', tlen);
		else if(tlen == 0)
			cigar_tmp.emplace_back('I', qlen);
		else
			cigar_tmp.emplace_back('M', qlen);
		if(ref_ed_ < ref_st_)//when ref_len < 0){ //add minus deletion, must follow the D/I/M
			cigar_tmp.emplace_back('D', ref_ed_ - ref_st_);
	}
	else{
		if(type == KSW_ALN_end_to_end){ //global
			read_score += ez.score;
			for(int i = ez.n_cigar - 1; i >= 0; i--)
				cigar_tmp.emplace_back(ez.cigar[i]);
		}
		else if(type == KSW_ALN_left_extend){ // left extension, soft clip at end
			read_score += (ez.mqe);
			for(int i = 0; i < ez.n_cigar; i++)
				cigar_tmp.emplace_back(ez.cigar[i]);
		}else{ //right extension, soft clip at begin
			read_score += (ez.mqe);
			for(int i = ez.n_cigar - 1; i >= 0; i--)
				cigar_tmp.emplace_back(ez.cigar[i]);
		}
	}
}

