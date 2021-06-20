/*
 * deBGA_index.hpp
 *
 *  Created on: 2021年3月15日
 *      Author: fenghe
 */

#ifndef DEBGA_INDEX_HPP_
#define DEBGA_INDEX_HPP_

#include <stdint.h>
#include <stdlib.h>
#include <string>
extern "C"
{
#include "clib/bam_file.h"
#include "clib/desc.h"
#include "clib/binarys_qsort.h"
}

#include <vector>
#include "cpp_lib/graph.hpp"

#define LEN_KMER 20 //kmer length used for search
#define ROUTE_LENGTH_MAX 1024
#define MAX_CHR_NAME_LENGTH 200
#define MAX_CHR_NUM 6000000
#define START_POS_REF 0

//MEM
struct vertex_MEM
{
	uint64_t uid;
	uint32_t seed_id;
	uint32_t read_pos;
	uint32_t uni_pos_off;
	uint32_t length;
	uint32_t pos_n;

	static int cmp(const void * a, const void * b)
	{
		vertex_MEM* sm1 = (vertex_MEM *)a;
		vertex_MEM* sm2 = (vertex_MEM *)b;

	    if(sm1->uid > sm2->uid)
	        return 1;
	    if(sm1->uid < sm2->uid)
	        return -1;
	    else
	    {
	        //according to read_pos and uni_pos_off detected if there is an inversion
	        if(sm1->read_pos > sm2->read_pos)
	            return 1;
	        if(sm1->read_pos < sm2->read_pos)
	            return -1;
	        else    return 0;
	    }
	}

	void print(){ fprintf(stderr, "vertex_MEM: uid:[%ld], seed_id:[%d], read_pos:[%d], uni_pos_off:[%d], length:[%d], pos_n:[%d] \n",
				uid, seed_id, read_pos, uni_pos_off, length, pos_n);}

};

//MEM joint region in uinpath
struct vertex_U
{
	uint64_t uid;
	uint32_t read_pos;
	uint32_t uni_pos_off;
	uint32_t length1; //length in read
	uint32_t length2;// length in reference
	uint32_t pos_n;
	uint32_t cov;

	void print(){
		fprintf(stderr, "vertex_U: uid:[%ld], read_pos:[%d], uni_pos_off:[%d], length1:[%d], length2:[%d], pos_n:[%d], cov:[%d] \n",
				uid, read_pos, uni_pos_off, length1, length2, pos_n, cov);}
};

struct SV_chr_info{
	SV_chr_info(char * chr_names, bam_hdr_t *ori_header){
		char target_name_tmp[1024];
		strcpy(target_name_tmp, chr_names);
		char * token = strtok(target_name_tmp, "_"); 	int id = atoi(token);
		token = strtok(NULL, "_"); 				int char_ID = bam_name2id(ori_header, token);		//char name
		token = strtok(NULL, "_"); 				uint32_t pos = atoi(token);

		token = strtok(NULL, "_"); 				int region_len = atoi(token);
		token = strtok(NULL, "_"); 				char * sv_type = token;

		token = strtok(NULL, "_"); 				int break_point1_pos =  atoi(token);//position of break point 1 in original reference
		token = strtok(NULL, "_"); 				int break_point2_pos =  atoi(token);//position of break point 2 in original reference
		token = strtok(NULL, "_"); 				int ed_pos =  atoi(token);//ending position of new reference in original reference
		token = strtok(NULL, "_"); 				char* vcf_id = token;//a name of vcf record in vcf file, like  "pbsv.INS.2"

		add_node(id, char_ID, pos ,
				region_len, sv_type,
				break_point1_pos ,break_point2_pos ,ed_pos , vcf_id);
	}

	SV_chr_info(uint32_t ID_, int32_t chr_ID_, uint64_t st_pos_, int region_len_, char * sv_type_,
			uint64_t break_point1_pos_, uint64_t break_point2_pos_, uint64_t ed_pos_, char * vcf_id_){
		add_node(ID_, chr_ID_, st_pos_,	region_len_,  sv_type_,	break_point1_pos_, break_point2_pos_, ed_pos_, vcf_id_);
	}

	void add_node(uint32_t ID_, int32_t chr_ID_, uint64_t st_pos_, int region_len_, char * sv_type_,
			uint64_t break_point1_pos_, uint64_t break_point2_pos_, uint64_t ed_pos_, char * vcf_id_){
		ID = ID_;
		chr_ID = chr_ID_;
		st_pos = st_pos_;

		region_len = region_len_;
		sv_type.append(sv_type_);

		break_point1_pos = break_point1_pos_;
		break_point2_pos = break_point2_pos_;
		ed_pos = ed_pos_;

		end_offset = ed_pos - st_pos - region_len;
		vcf_id.append(vcf_id_);

		//like: SV:Z:3_0_446371_1108_INS_pbsv.INS.2
		char vcf_print_string_buff[1000];
		sprintf(vcf_print_string_buff, "%d_%d_%ld_%d_%s_%s", ID, chr_ID, st_pos, region_len, sv_type_, vcf_id_);
		vcf_print_string.append(vcf_print_string_buff);

		//original_read_region_2_adjust_length
		if(sv_type.compare("INS") == 0)
			original_read_region_2_adjust_length = region_len - (break_point1_pos - st_pos) - (ed_pos - break_point2_pos);
		else if(sv_type.compare("DEL") == 0)
			original_read_region_2_adjust_length = break_point1_pos - break_point2_pos;
		else
			original_read_region_2_adjust_length = 0;
	}

	int getDeletionLen(){
		int del_len = break_point2_pos - break_point1_pos;
		if(del_len < 10)
			del_len = 0;
		return del_len;
	}

	uint32_t ID;
	uint32_t chr_ID;
	uint64_t st_pos;

	int region_len;
	std::string sv_type;

	//position for other important point
	uint64_t break_point1_pos;//position of break point 1 in original reference
	uint64_t break_point2_pos;//position of break point 2 in original reference
	uint64_t ed_pos;//ending position of new reference in original reference

	//original read region 2 adjust length
	int original_read_region_2_adjust_length;

	int end_offset;
	std::string vcf_id;//a name of vcf record in vcf file, like  "pbsv.INS.2"
	std::string vcf_print_string;
};

struct deBGA_INDEX{

	uint64_t* buffer_ref_seq = NULL;//original reference
	uint64_t* buffer_seq = NULL;// UNITIG sequence

	uint64_t* buffer_seqf = NULL;//start offset of [UINTIG] in [buffer_seq]
	uint64_t* buffer_off_g = NULL;//start offset of [KMER] in [buffer_seq]
	uint64_t* buffer_p = NULL;//start offset of [UINTIG] in [buffer_ref_seq]
	uint64_t* buffer_pp = NULL;//start index of [UINTIG] in [buffer_p], used to search offset of [UINTIG] in [buffer_ref_seq]
	uint64_t* buffer_hash_g = NULL;//index for the first 14 base pair, used to search kmer in the hash table
	uint32_t* buffer_kmer_g = NULL;//used to store other part of a kmer(except the first 14 byte)

	uint32_t chr_search_index_size = 0;
	uint32_t * chr_search_index = NULL;

	uint64_t reference_len = 0;
	uint32_t chr_end_n[MAX_CHR_NUM];
	char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
	char chr_line_content[MAX_CHR_NAME_LENGTH];
	int chr_file_n = 1;

	std::vector<SV_chr_info> sv_info;

	bam_hdr_t * header = NULL;
	bam_hdr_t * ori_header = NULL;

	uint64_t result_ref_seq = 0;
	uint64_t result_seq = 0;
	uint64_t result_seqf = 0;
	uint64_t result_p = 0;
	uint64_t result_pp = 0;
	uint64_t result_pu = 0;
	uint64_t result_hash_g = 0;
	uint64_t result_kmer_g = 0;
	uint64_t result_off_g = 0;
	uint64_t result_ref_g = 0;

	int load_index_file(char *index_dir);

	//search kmer in the index, the search result "index of kmer" stored in range, range[0] is the start index, and range[1] is the end index
	//return -1 when search failed, otherwise return 0
	bool search_kmer(int len_k, uint64_t kmer, int64_t range[2], int8_t seed_offset);
	//function: given an index of kmer, search the MEM within a UNITIG
	//return 0 when successfully stored data, -1 when failed
	int UNITIG_MEM_search(uint64_t kmer_index, std::vector<vertex_MEM>& vertexm_v, uint64_t *read_bit, uint32_t read_off, uint32_t read_length, int len_k, uint32_t & max_right_i);
	//function: merge all MEMs within a unipath, the result stored in vertexu_v
	void merge_seed_in_unipath(std::vector<vertex_MEM>& vertexm_v, std::vector<vertex_U>& vertexu_v);
	//store seed in uintig into seed in reference
	void expand_seed(std::vector<vertex_U>& vertexu_v, std::vector<UNI_SEED>& uniseed_v, random_data*  rand_num);

	void get_refseq(uint8_t *ref, uint32_t len, uint32_t start);
	void get_ref_jump(uint8_t *ref, uint32_t start1, uint32_t len1, uint32_t start2, uint32_t len2);
	void get_ref_onebyone(uint8_t *ref, uint32_t start, uint32_t len, uint32_t pre_pos);

	void building_chr_index();
	//get chr ID and begin
	int get_chromosome_ID(uint32_t position);

	void building_bam_header();

	void free_memory();

};
uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, bool useCalloc, int additional_byte);
int build_deBGA_index(int argc, char *argv[]);

#endif /* DEBGA_INDEX_HPP_ */
