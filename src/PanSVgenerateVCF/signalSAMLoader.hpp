/*
 * signalSAMLoader.hpp
 *
 *  Created on: 2021年5月31日
 *      Author: fenghe
 */

#ifndef PANSVGENERATEVCF_SIGNALSAMLOADER_HPP_
#define PANSVGENERATEVCF_SIGNALSAMLOADER_HPP_

#include "deBGA_index.hpp"
#include <vector>

#include "generateVCFoptions.hpp"
#include "samTag.hpp"
extern "C"
{
#include "../clib/bam_file.h"
}

struct assembly_sorter_item{
	int bam_id;
	int SV_ID;
	int st_pos_ref;

	void set(int bam_id_, int SV_ID_, int st_pos_ref_){
		bam_id = bam_id_;
		SV_ID = SV_ID_;
		st_pos_ref = st_pos_ref_;
	}

	static inline int cmp(const void*a_,const void*b_)
	{
		assembly_sorter_item *a = (assembly_sorter_item *)a_;
		assembly_sorter_item *b = (assembly_sorter_item *)b_;

		if(a->SV_ID != b->SV_ID)
			return a->SV_ID > b->SV_ID;
		if(a->st_pos_ref != b->st_pos_ref)
			return a->st_pos_ref > b->st_pos_ref;
		else
			return a->bam_id > b->bam_id;
	}
};

struct Remapped_read_loader{
public:

	//output
	bam1_t *read_load_buff;
	assembly_sorter_item * sort_l;

	//
	//return read number, when no region left, return -1;
	//the "st_idx" is the start index of read list
	bool get_next_sv_region(int * st_idx, int *read_cnt, int *SV_ID){
		if((uint32_t)current_sv_region + 1 >= sv_region_st.size()){
			load_block_sort();
			current_sv_region = 0;
			if(sv_region_st.size() == 0){
				return false;
			}
		}
		*st_idx = sv_region_st[current_sv_region];
		*read_cnt = sv_region_st[current_sv_region + 1] - sv_region_st[current_sv_region];
		*SV_ID = first_SV_ID + current_sv_region;
		current_sv_region++;
		return true;
	}

	bool get_sv_region_by_ID(int * st_idx, int *read_cnt, int SV_ID){
		int c_sv_region = SV_ID - first_SV_ID;
		if((uint32_t)c_sv_region + 1 >= sv_region_st.size()){ *st_idx = -1; *read_cnt = 0; return false; } //out of range, return false
		*st_idx = sv_region_st[c_sv_region];
		*read_cnt = sv_region_st[c_sv_region + 1] - sv_region_st[c_sv_region];
		return true;
	}

	int init(SignalAssemblyOption *o_, std::vector<SV_chr_info> *sv_info_list){
		o = o_;//set option
		all_finish = false;
		chr_finish = true;

		bam_file_open(o->input_bam_fn, NULL, NULL, &c_b);
		hdr = c_b._hdr;
		//read loading buff
		read_load_buff		=      (bam1_t *)xcalloc(o->max_load_read_per_time, sizeof(bam1_t));
		sort_l				= (assembly_sorter_item *)xcalloc(o->max_load_read_per_time, sizeof(assembly_sorter_item));

		//chr_ID check
		o->st_chr_ID = MAX(o->st_chr_ID, 0);
		o->ed_chr_ID = MIN(o->ed_chr_ID, hdr->n_targets);
		current_chr_ID = o->st_chr_ID - 1;

		current_sv_region = 0;
		total_load_number = 0;
		sv_info_list_p = (sv_info_list);
		return 0;
	}

	void destory(){
		bam_file_close(&c_b);
		for(int i = 0; i < o->max_load_read_per_time; i++){
			if(read_load_buff[i].data != NULL){
				 free(read_load_buff[i].data);
				 read_load_buff[i].data = NULL;
			}
		}
		free(read_load_buff);
	}

	std::vector<SV_chr_info> * sv_info_list_p;

private:

	//function in step 1, Bam_file *cf, bam1_t * b, int n_max_load can`t be global, for thread safety reason
	static int load_bam(Bam_file *cf, bam1_t * b, assembly_sorter_item * sort_l, int n_max_load, int MIN_score, std::vector<SV_chr_info> *sv_info_list)
	{
		int i, rst, score;
		for( i = 0; i < n_max_load; i++){
			rst = bam_next_dump(cf, b + i);
			if(rst <= 0) break;
			bam1_t *br = b + i;
			char* SV_tag_char = bam_get_string_tag(br, vcf_tag.SV);
			int chaininig_score = 0;
			bool has_SC = bam_get_num_tag(br, vcf_tag.CS, &chaininig_score);
			bam_get_num_tag(br, vcf_tag.AS, &score);
			if(score < MIN_score || SV_tag_char == NULL || (!has_SC && br->core.isize == 0)){
				i--;
				continue;
			}
			//CIGAR adjust:
			//adjust timing: before store sort item
			//adjust SAM: when M number less 2 and followed by long D(ignore I and other cigar type), like: 25I2M55D121M; 60D148M; etc.
			//adjust method: change small M into I, change length of D into 0, pos += (len_small_M + len_D)
			//original: 25I2M55D121M (pos: 1913086) ----> 25I2I0D121M (pos:1913143)
			uint32_t n_cigar = br->core.n_cigar;
			br->core.pos += cigar_adjust(&n_cigar, bam_get_cigar(br), true, 4);
			br->core.n_cigar = n_cigar;

			//get SV ID
			char * save_ptr;
			strtok_r(SV_tag_char, "_", &save_ptr);
			int SV_ID = atoi(SV_tag_char);
			save_ptr[-1] = '_';

			//position adjust for original mapping reads
			int chaining_score = -1;
			bam_get_num_tag(br, vcf_tag.CS, &chaining_score);
			if(chaining_score == -1 && br->core.pos > (int)sv_info_list[0][SV_ID].break_point2_pos && br->core.pos < (int)sv_info_list[0][SV_ID].ed_pos)
				br->core.pos += sv_info_list[0][SV_ID].original_read_region_2_adjust_length;

			//store sort item
			sort_l[i].set(i, SV_ID, br->core.pos);
		}
		return i;
	}

	bool load_block_sort(){
		sv_region_st.clear();
		cur_load_numer = 0;
		if(all_finish) return all_finish;
		if(chr_finish){// load new region
			current_chr_ID ++;
			if(current_chr_ID >  o->ed_chr_ID){	all_finish = true; return all_finish;	}
			chr_finish = false;
			int tar_len = hdr->target_len[current_chr_ID];
			R_region r;	r.chr_ID = current_chr_ID;
			r.st_pos = (current_chr_ID == o->st_chr_ID)?o->st_pos:0;
			r.ed_pos = (current_chr_ID == o->ed_chr_ID)?o->ed_pos:tar_len;//2954350
			r.st_pos = MAX(r.st_pos, 0);
			r.ed_pos = MIN(r.ed_pos, tar_len);
			resetRegion_ID(&c_b, &r);
		}

		//load reads
		int max_load = MIN(o->max_read - total_load_number, o->max_load_read_per_time);
		cur_load_numer = load_bam(&c_b, read_load_buff, sort_l, max_load, o->MIN_score, sv_info_list_p);
		//fprintf(stderr, "cur_load_numer: %d\n", cur_load_numer);
		if(cur_load_numer < 10)	{chr_finish = true;	if(!all_finish) return load_block_sort();}
		total_load_number += cur_load_numer;
		qsort(sort_l, cur_load_numer, sizeof(assembly_sorter_item), assembly_sorter_item::cmp);
		sv_region_st.emplace_back(0);
		first_SV_ID = sort_l[0].SV_ID;
		for(int read_idx = 0; read_idx < cur_load_numer - 1; read_idx++){
			if(sort_l[read_idx].SV_ID != sort_l[read_idx + 1].SV_ID){
				int store_number = sort_l[read_idx + 1].SV_ID - sort_l[read_idx].SV_ID;
				for(int c_sv = 0; c_sv < store_number; c_sv++) sv_region_st.emplace_back(read_idx + 1);
			}
		}
		sv_region_st.emplace_back(cur_load_numer);

		return all_finish;
	}

	//options:
	SignalAssemblyOption *o;
	int total_read;

	//states
	int cur_load_numer; //current read number in buff
	int current_sv_region;
	bool all_finish;
	bool chr_finish;
	int current_chr_ID;
	int total_load_number;
	int first_SV_ID;//the first sv id for "sv_region_st" list

	//bam file
	Bam_file c_b;
	bam_hdr_t* hdr;

	//buffs
	std::vector<int> sv_region_st;
};



#endif /* PANSVGENERATEVCF_SIGNALSAMLOADER_HPP_ */
