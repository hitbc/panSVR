/*
 * NovaSVRst.hpp
 *
 *  Created on: 2021年8月20日
 *      Author: fenghe
 */

#ifndef NOVASVGENERATEVCF_NOVASVRST_HPP_
#define NOVASVGENERATEVCF_NOVASVRST_HPP_

#include <string>
#include <vector>
#include <cstring>
#include "RefHandler.hpp"

extern "C"{
	#include "../clib/utils.h"
	#include "../clib/bam_file.h"
	#include "../kswlib/kalloc.h"
	#include "../kswlib/ksw2.h"
}

struct Genotyping_read_aligner{

public:
	//as input
	void init(){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option();
		ksw_gen_mat_D();
	}

	void setRef(uint8_t * ori_ref_string, int ori_ref_len_, const char * contig_string, int contig_len)
	{
		bin_contig.resize(contig_len);
		//change char contig into bin contig
		for (int i = 0; i < contig_len; ++i)
			bin_contig[i] = charToDna5n[(uint8_t)contig_string[i]];
		ori_ref = ori_ref_string;
		ori_ref_len = ori_ref_len_;
	}

	void destory(){
		if(ez.cigar != NULL)
			free(ez.cigar);
		km_destroy(km);
	}

	void align_non_splice(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		//simple check:
		if(0){
			//debug code:
			uint32_t debug_qlen = qlen;
			uint8_t* debug_qseq = qseq;
			for(uint32_t i = 0; i < debug_qlen; i++)
				fprintf(stderr, "%c", "ACGT"[ debug_qseq[i]]);
			fprintf(stderr, "\n");
			uint32_t debug_tlen = ref_end_pos - ref_st_pos;
			uint8_t* debug_tseq =  tseq + ref_st_pos;
			for(uint32_t i = 0; i < debug_tlen; i++)
				fprintf(stderr, "%c", "ACGT"[ debug_tseq[i]]);
			fprintf(stderr, "\n");
		}
		ksw_extd2_sse(km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

	int get_score_reach_end_of_read(){ return MAX(0, ez.mqe);}

	//when align_to_contig == true, aligned to contig, otherwise, aligned to refernece
	void align_genotyping(bool align_to_contig, uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		if(align_to_contig){
			tseq = &(bin_contig[0]);
			tlen = bin_contig.size();
		}else{
			tseq = ori_ref;
			tlen = ori_ref_len;
		}
		align_non_splice(qseq_, qlen_, ref_st_pos, ref_end_pos);
	}

	int gap_penalty(int gap_len){
		int penalty1 = gap_open_D + gap_len*gap_ex_D;
		int penalty2 = gap_open2_D + gap_len*gap_ex2_D;
		return MIN(penalty1, penalty2);
	}
	int getScoreByCigar_with_skip_region(bam1_t *br, int read_skip_left, int read_skip_right, uint8_t * qseq_buff, RefHandler *refHandler){

		int read_left_boundary = read_skip_left;
		int read_right_boundary = br->core.l_qseq - read_skip_right;

		int score = 0;
		uint32_t* bam_cigar = bam_get_cigar(br);
		int match_base, mis_match_base;
		int q_seq_idx = 0;
		int t_seq_idx = 0;
		uint8_t *qseq_str = qseq_buff;
		uint8_t *tseq_str = refHandler->getRefStr(br->core.pos);
		for (uint i = 0; i < br->core.n_cigar; ++i)
		{
			int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			switch (c_type){
			case CIGAR_MATCH:
			case CIGAR_SEQ_MATCH:
				match_base = 0;
				mis_match_base = 0;
				for(int i = 0; i < c_size; i++, q_seq_idx++, t_seq_idx++){
					if(q_seq_idx < read_left_boundary || q_seq_idx >= read_right_boundary)
						continue;
					if(qseq_str[q_seq_idx] != tseq_str[t_seq_idx])
						mis_match_base++;
					else
						match_base++;
				}
				score += ((match_base * match_D) - (mis_match_base * mismatch_D)); break;
			case CIGAR_INSERT:	case CIGAR_SOFT_CLIP:	case CIGAR_HARD_CLIP:
				mis_match_base = 0;
				for(int i = 0; i < c_size; i++, q_seq_idx++){
					if(q_seq_idx < read_left_boundary || q_seq_idx >= read_right_boundary)
						continue;
					mis_match_base++;
				}
				score -= gap_penalty(mis_match_base); break;
			case CIGAR_DELETE:
				t_seq_idx += c_size;
				if(q_seq_idx < read_left_boundary || q_seq_idx >= read_right_boundary)
					break;
				score -= gap_penalty(c_size); break;
			default:
				break;
			}
		}
		score = MAX(0, score);
		return score;
	}

	int getScoreByMismatch(int search_length, int mismatch_num){
		return ((search_length - mismatch_num) * match_D) - (mismatch_num * mismatch_D);
	}

	int adjustCIGAR(){
		uint32_t n_cigar = ez.n_cigar;
		int adj_size = cigar_adjust(&n_cigar, ez.cigar, false, 15);
		ez.n_cigar = n_cigar;
		return adj_size;
	}

	void printf_alignment_detail(FILE * output, int suggest_st_pos){
		log_output = output;
		int cigar_len = ez.n_cigar;
		fprintf(output," CIGAR number: %d ", cigar_len);
		if(cigar_len == 0) return;
		fprintf(output," \tCIGAR: ", cigar_len);
		uint32_t* bam_cigar =  ez.cigar;
		 for (int i = 0; i < cigar_len; ++i){
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			char c_type = segment_type_to_cigar_code(type);
			int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			fprintf(output,"%d%c",length,  c_type);
		}
		fprintf(output," ");
		if(0){
			print_X_E_sequence(suggest_st_pos); //print X/= sequence
			fprintf(output, "\n");
		}
	}

	inline uint8_t getTseq(int i){ return tseq[i];}

	void setZdrop(uint16_t zdrop_D_, int bandwith_){
		zdrop_D = zdrop_D_;
		bandwith = bandwith_;
	}

	uint8_t * get_bin_contig(){ return &(bin_contig[0]); }
	uint8_t * get_ori_ref(){ return ori_ref; }

	int getGenotypingMinScore(int used_read_length){
		int min_score = (used_read_length -80) * match_D;
		int global_min_match_score = 50 * match_D;
		return MAX(global_min_match_score, min_score);
	}

private:
	//part2: reference and target
	uint8_t *tseq; int tlen;
	uint8_t *qseq; uint32_t qlen; //query temp: used for analysis
	//part3: alignment result
	ksw_extz_t ez;
	//part4: buffs
	void *km;
	//part5: options
	int8_t mata_D[25]; int8_t match_D; int8_t mismatch_D;
	int8_t gap_open_D; int8_t gap_ex_D; int8_t gap_open2_D; int8_t gap_ex2_D;
	uint16_t zdrop_D; int bandwith;
	int flag;
	//part6: logs
	FILE * log_output;
	std::vector<uint8_t> bin_contig;

	uint8_t* ori_ref; int ori_ref_len;

	void copy_option(){
			match_D = 2;
			mismatch_D= 6;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 1;
			zdrop_D= gap_open2_D + 30;
			bandwith = 30;//zdrop_D;
			flag = 0;
	}

	void ksw_gen_mat_D(){
		int8_t l,k,m;
		for (l = k = 0; l < 4; ++l) {
			for (m = 0; m < 4; ++m) { mata_D[k] = l == m ? match_D : -(mismatch_D);	/* weight_match : -weight_mismatch */ k++; }
			mata_D[k] = 0; // ambiguous base
			k++;
		}
		for (m = 0; m < 5; ++m) { mata_D[k] = 0; k++; }
	}

	int print_X_E_sequence(int suggest_st_pos){
		int contig_coverage = 0;
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:
				for(int i = 0; i < cigar_len; i++, seq_i++)
				{
					if(qseq[seq_i] == tseq[output_index]){
						contig_coverage++;
						fprintf(log_output, "=");//1 == M; 2 == X; 0 == -;
					}
					else{
						fprintf(log_output, "X");//1 == M; 2 == X; 0 == -;
					}
					output_index++;
				}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\n");
		return contig_coverage;
	}
};

bool bam_aligned_analysis(bam1_t *b, int *clip_left, int *clip_right, int *gap_mismatch_inside);

struct NOVA_SV_FINAL_RST_item{

public:
	//new a blank nodes
	NOVA_SV_FINAL_RST_item(){
		chr_ID = 0; st_pos = 0; final_genotype = 0;
		suggenst_SV_length = 0;	SV_length = 0;	endPos = 0;
		will_be_output_to_vcf = false;
		dis_SV_len_suggset_length = 0;
		cigar_idx_bg = 0;
		cigar_idx_ed = 0;
		contig_st_pos_in_ref = 0;
		region_ref_global_position = 0;
		region_is_overlap = false;
		deletion_pass_DR_filter = true;
		is_duplication_sv = false;
		is_BND = false;
	}

	static inline int cmp_by_position(const NOVA_SV_FINAL_RST_item &a, const NOVA_SV_FINAL_RST_item &b){
		if(a.chr_ID == b.chr_ID)	return a.st_pos < b.st_pos;
		else							return a.chr_ID < b.chr_ID;
	}

	NOVA_SV_FINAL_RST_item(int chr_ID_, int st_pos_, const char * SV_type_,const char * ref_,const char * alt_,
			int suggest_sv_len, uint8_t * contig_, int contig_len, uint32_t *cigar_, int cigar_num, int cigar_idx_bg_, int cigar_idx_ed_, int contig_st_pos_in_ref_, int region_ref_global_position_){
		//basic informations
		chr_ID = chr_ID_; st_pos = st_pos_; final_genotype = 0;
		ref.clear(); alt.clear(); SV_type.clear();
		SV_type.append(SV_type_);
		SV_name.append("TODO_NO_NAME");
		ref.append(ref_);
		alt.append(alt_);
		SV_length = alt.size() - ref.size() + 1;
		endPos = st_pos + ref.size();

		//additional informations
		suggenst_SV_length = suggest_sv_len;
		will_be_output_to_vcf = false;
		dis_SV_len_suggset_length = ABS_U(SV_length, suggenst_SV_length);
		contig.resize(contig_len);
		for(int i = 0; i < contig_len; i++) contig[i] = "ACGT"[contig_[i]];
		cigar.clear();
		for(int i = 0; i < cigar_num; i++)
			cigar.emplace_back(cigar_[i]);
		cigar_idx_bg = cigar_idx_bg_;
		cigar_idx_ed = cigar_idx_ed_;
		contig_st_pos_in_ref = contig_st_pos_in_ref_;
		region_ref_global_position = region_ref_global_position_;
		region_is_overlap = false;
		deletion_pass_DR_filter = true;
		is_duplication_sv = false;
		region_support_numnber[2][0] = 0;
		region_support_numnber[2][1] = 0;
		region_support_numnber[2][2] = 0;
		is_BND = false;
	}

	static void add_to_vector(std::vector<NOVA_SV_FINAL_RST_item> & v, int chr_ID_, int st_pos_, const char * SV_type_, const char * ref_,const char * alt_,
			int suggest_sv_len, uint8_t * contig_, int contig_len, uint32_t *cigar, int cigar_num, int cigar_idx_bg_, int cigar_idx_ed_, int contig_st_pos_in_ref, int region_ref_global_position_){
		v.emplace_back(chr_ID_, st_pos_, SV_type_, ref_, alt_,
				suggest_sv_len,  contig_, contig_len, cigar, cigar_num, cigar_idx_bg_, cigar_idx_ed_, contig_st_pos_in_ref, region_ref_global_position_);
	}

	//store BND mode:
	NOVA_SV_FINAL_RST_item(int chr_ID_, int st_pos_, const char * ref_,const char * alt_, const char * BND_type,
			int region_support_numnber_read, int region_support_numnber_ref, int region_support_numnber_unknown,
			int genotype){
		//basic informations
		chr_ID = chr_ID_; st_pos = st_pos_;
		final_genotype = genotype;
		ref.clear(); alt.clear(); SV_type.clear();
		SV_type.append("BND");
		SV_name.append("TODO_NO_NAME");
		ref.append(ref_);
		alt.append(alt_);
		bnd_type.append(BND_type);
		SV_length = 99999;
		endPos = 0;

		//additional informations
		suggenst_SV_length = 0;
		will_be_output_to_vcf = true;
		dis_SV_len_suggset_length = 0;
		cigar_idx_bg = 0;
		cigar_idx_ed = 0;
		contig_st_pos_in_ref = 0;
		region_ref_global_position = 0;
		region_is_overlap = false;
		deletion_pass_DR_filter = true;
		is_duplication_sv = false;
		region_support_numnber[2][0] = region_support_numnber_read;
		region_support_numnber[2][1] = region_support_numnber_ref;
		region_support_numnber[2][2] = region_support_numnber_unknown;
		is_BND = true;
	}

	static void store_BND(std::vector<NOVA_SV_FINAL_RST_item> & v, int chr_ID_, int st_pos_, const char * ref_,const char * alt_, const char * BND_type,
			int region_support_numnber_read, int region_support_numnber_ref, int region_support_numnber_unknown,
			int genotype){
		v.emplace_back(chr_ID_, st_pos_, ref_, alt_, BND_type,
				region_support_numnber_read,  region_support_numnber_ref, region_support_numnber_unknown, genotype);
	}

	void writeVCF(kstring_t *s, bam_hdr_t * header){
		s->l = 0;
		//show results
		sprintf(s->s + s->l, "%s\t%d\t%s\t", header->target_name[chr_ID], st_pos, SV_name.c_str()); s->l += strlen(s->s + s->l);
		sprintf(s->s + s->l, "%s\t%s\t", 	 ref.c_str(), alt.c_str()); s->l += strlen(s->s + s->l);

		//filters:
		if(SV_length < 50 && SV_length > -50)			sprintf(s->s + s->l, ".\t%s\t", "lt50bp");
		else if(is_duplication_sv)						sprintf(s->s + s->l, ".\t%s\t", "Duplication_sv");
		else if(final_genotype == 0)					sprintf(s->s + s->l, ".\t%s\t", "LOW_DEPTH");
		else if(!deletion_pass_DR_filter)				sprintf(s->s + s->l, ".\t%s\t", "DR_not_support");
		else											sprintf(s->s + s->l, ".\t%s\t", "PASS");
		s->l += strlen(s->s + s->l);
		if(is_BND){
			sprintf(s->s + s->l, "SVTYPE=%s;BND_TYPE=%s\t" , SV_type.c_str(), bnd_type.c_str()); s->l += strlen(s->s + s->l);
		}else{
			sprintf(s->s + s->l, "SVTYPE=%s;END=%d;SVLEN=%d\t" , SV_type.c_str(), endPos, SV_length); s->l += strlen(s->s + s->l);
		}
		sprintf(s->s + s->l, "GT:SR:SG\t" "%s:", (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1")); s->l += strlen(s->s + s->l);
		sprintf(s->s + s->l, "%d,%d,%d:%d", region_support_numnber[2][1],  region_support_numnber[2][0], region_support_numnber[2][2], suggenst_SV_length); s->l += strlen(s->s + s->l);
		sprintf(s->s + s->l, "\n"); s->l += strlen(s->s + s->l);
	}

	void printContigSeq(int suggest_st_pos, FILE * log){
		uint32_t* bam_cigar = &(cigar[0]);
		int cigar_len = cigar.size();
		int output_index = 0;
		const char * qseq = contig.c_str();
		int seq_i = 0;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "%c", qseq[seq_i]); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "-"); output_index++;}	break;//S, print -
			default: fprintf(log, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log, "\n");
	}

	void print_X_E_sequence(int suggest_st_pos, uint8_t * tseq, FILE * log){
		int output_index = 0;
		int seq_i = 0;
		const char * qseq = contig.c_str();
		uint32_t* bam_cigar = &(cigar[0]);
		int cigar_len = cigar.size();
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:
				for(int i = 0; i < cigar_len; i++, seq_i++){ fprintf(log, "%c", (qseq[seq_i] == "ACGT"[tseq[output_index]])?'=':'X'); output_index++;}
				break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "-"); output_index++;}	break;//S, print -
			default: fprintf(log, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log, "\n");
	}

	void print(FILE * log,  bam_hdr_t * header, RefHandler *refHandler){
		fprintf(log, "%s\t%d\t%s\t", header->target_name[chr_ID], st_pos, SV_name.c_str());
		fprintf(log, "%s\t%s\t", 	 ref.c_str(), alt.c_str());
		fprintf(log, ".\t%s\tSVTYPE=%s;END=%d;SVLEN=%d\t" "GT:DP:SG\t" "%s:", (final_genotype == 0)?"LOW_DEPTH":"PASS" ,SV_type.c_str(), endPos, SV_length , (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1"));

		fprintf(log, "%d,%d:%d", region_support_numnber[2][0],  region_support_numnber[2][1], suggenst_SV_length);
		fprintf(log, "\nother info: contig: %s; \ncigar_idx: [%d, %d]; contig in ref: %d, ref in global %d \t", contig.c_str(), cigar_idx_bg, cigar_idx_ed, contig_st_pos_in_ref, region_ref_global_position);

		for(uint32_t c: cigar)
			fprintf(log, "%d%c", (c >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(c & BAM_CIGAR_MASK)]);
		fprintf(log, "\n");

		printContigSeq(contig_st_pos_in_ref, log);
		print_X_E_sequence(contig_st_pos_in_ref, refHandler->getRefStr(region_ref_global_position), log);

		fprintf(log, "\n");
	}

	static void resultFilter(std::vector<NOVA_SV_FINAL_RST_item> & rst_l, int MIN_sv_len){
		int smallest_SV_dis = MAX_int32t;
		int rst_num = rst_l.size();
		for(int i = 0; i < rst_num; i++)
			if(rst_l[i].suggenst_SV_length <= -MIN_sv_len || rst_l[i].suggenst_SV_length >= MIN_sv_len)
				{ smallest_SV_dis = MIN(smallest_SV_dis, rst_l[i].dis_SV_len_suggset_length); }
		if(smallest_SV_dis == MAX_int32t) return;
		std::map<int, int> duplication_count;
		for(int i = 0; i < rst_num; i++)
			if((rst_l[i].suggenst_SV_length <= -40 || rst_l[i].suggenst_SV_length >= 40) && rst_l[i].dis_SV_len_suggset_length <= smallest_SV_dis){
				rst_l[i].will_be_output_to_vcf = true;
				std::map<int, int>::iterator find_rst = duplication_count.find(rst_l[i].suggenst_SV_length);
				if(find_rst == duplication_count.end())
					duplication_count[rst_l[i].suggenst_SV_length] = 1;
				else
					find_rst->second++;
			}
		//remove duplications
		for(auto &dup: duplication_count){
			if(dup.second > 1){
				int c_suggenst_SV_length = dup.first;
				int c_index = 0;
				for(int i = 0; i < rst_num; i++){
					if(rst_l[i].will_be_output_to_vcf && rst_l[i].suggenst_SV_length == c_suggenst_SV_length){
						if(c_index > 0)
							rst_l[i].will_be_output_to_vcf = false;
						c_index++;
					}
				}
			}
		}
		//remove other wrong results
		for(int i = 0; i < rst_num; i++)
			if(rst_l[i].will_be_output_to_vcf && ((rst_l[i].SV_length < 0 && rst_l[i].suggenst_SV_length > 0) || (rst_l[i].SV_length > 0 && rst_l[i].suggenst_SV_length < 0)))
				rst_l[i].will_be_output_to_vcf = false;
	}

	//return:
	//0: not overlap with both
	//1: overlap with bp1
	//2: overlap with bp2
	//3: overlap with both

	int read_overlap_breakpoint(bam1_t *br, int region_ID, bool with_supp, int normal_read_length, int * left_clip_output){

		int breakpoint1 = st_pos;
		int breakpoint2 = endPos;
		*left_clip_output = 0;
		//check whether reads overlap with breakpoint
		int read_st_pos = -1;
		if(br->core.n_cigar <= 1)
			read_st_pos = br->core.pos;
		else{
			int clip_left, clip_right, gap_mismatch_inside; //the results of analysis alignment
			bam_aligned_analysis(br, &clip_left, &clip_right, &gap_mismatch_inside);
			read_st_pos = br->core.pos - clip_left;
			*left_clip_output = clip_left;
		}
		bool read_overlap_with_breakpoint1 = (read_st_pos <= breakpoint1 && read_st_pos + normal_read_length > breakpoint1);
		bool read_overlap_with_breakpoint2 = (read_st_pos <= breakpoint2 && read_st_pos + normal_read_length > breakpoint2);
		if(0){
			if(read_overlap_with_breakpoint1)	fprintf(stderr, "overlap bp1\t"); else	fprintf(stderr, "NOTover bp1\t");
			if(read_overlap_with_breakpoint2)	fprintf(stderr, "overlap bp2\t"); else  fprintf(stderr, "NOTover bp2\t");
		}
		int overlap_mode = 0;

		if(with_supp){
			if(region_ID == 0 && read_overlap_with_breakpoint1)
				overlap_mode = read_overlap_with_breakpoint1;
			if(region_ID == 1 && read_overlap_with_breakpoint2 )
				overlap_mode = read_overlap_with_breakpoint2 * 2;
		}else
			overlap_mode = read_overlap_with_breakpoint1 + read_overlap_with_breakpoint2 * 2;

		return overlap_mode;
	}

	int get_contig_golbal_position_core(bool is_bp1){
		int used_cigar_num = cigar_idx_bg;
		if(is_bp1 == false)
			used_cigar_num = cigar_idx_ed + 1;
		//suggest contig_region_st for BP1
		int contig_region_st = 0;
		int small_indel_count = 0;

		for(int c_idx = 0; c_idx < used_cigar_num; c_idx++){
			int c_size = cigar[c_idx] >> BAM_CIGAR_SHIFT;
			int c_type = cigar[c_idx] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:
				//contig_region_st += c_size;
				if(c_size < 20){ small_indel_count -= c_size; } else small_indel_count = 0;
				break;//M
			case 1:
				contig_region_st += c_size;
				if(c_size < 20){ small_indel_count -= c_size; } else small_indel_count = 0;
				break;//I
			case 2:
				contig_region_st -= c_size;
				if(c_size < 20){ small_indel_count += c_size; } else small_indel_count = 0;
				break;//M or D
			default: break;
			}
		}
		contig_region_st -= small_indel_count;

		return region_ref_global_position + contig_st_pos_in_ref - contig_region_st;
	}

	void get_contig_golbal_position(int *contig_pos_bp1, int *contig_pos_bp2){
		*contig_pos_bp1 = get_contig_golbal_position_core(true);
		*contig_pos_bp2 = get_contig_golbal_position_core(false);
		fprintf(stderr, "SUG:contig:bp1 %d; SUG:contig:bp2 %d\n", *contig_pos_bp1, *contig_pos_bp2);
	}

	int get_ori_alignment_score(Genotyping_read_aligner * gra, bam1_t *br, uint8_t * qseq_buff, int read_skip_left, int read_skip_right, RefHandler *refHandler){
		int score = 0;
		if(0){//running re-alignment
			//re-alignment and get score for both ref-region
			int read_in_ref_st_pos = br->core.pos - region_ref_global_position;
			int read_in_ref_ed_pos = read_in_ref_st_pos + br->core.l_qseq + 50;
			gra->align_genotyping(false, qseq_buff + read_skip_left, br->core.l_qseq - read_skip_left - read_skip_right, read_in_ref_st_pos, read_in_ref_ed_pos);
			//gra->adjustCIGAR();
			gra->printf_alignment_detail(stderr, read_in_ref_st_pos);
			return 0;
		}else{//direct get score from cigar:
			score =  gra->getScoreByCigar_with_skip_region(br, read_skip_left, read_skip_right, qseq_buff, refHandler);
		}
		return score;
	}

	int get_contig_alignment_score_core(Genotyping_read_aligner * gra, bam1_t *br, uint8_t * qseq_buff, int read_in_contig_st_pos, int *read_skip_left_, int *read_skip_right_){
		//realigned for contig regions
		int read_skip_left = 0;
		int read_skip_right = 0;
		if(read_in_contig_st_pos < 0){
			read_skip_left -= read_in_contig_st_pos;
			read_in_contig_st_pos = 0;
		}
		int read_in_contig_ed_pos = read_in_contig_st_pos + br->core.l_qseq;
		if(read_in_contig_ed_pos > (int)contig.size()){
			read_skip_right = read_in_contig_st_pos + br->core.l_qseq - contig.size();
			read_in_contig_ed_pos = contig.size();
		}

		*read_skip_left_ = read_skip_left;
		*read_skip_right_ = read_skip_right;

		int qlen = br->core.l_qseq - read_skip_left;
		uint8_t * qseq = qseq_buff + read_skip_left;
		int tar_len = read_in_contig_ed_pos - read_in_contig_st_pos;
		uint8_t* tar_str = gra->get_bin_contig() + read_in_contig_st_pos;
		int search_len = MIN(qlen, tar_len);
		int wrong_base = 0;
		for(int i = 0; i < search_len && wrong_base < 6; i++)
			if(tar_str[i] != qseq[i]) wrong_base ++;
		if(wrong_base < 6){
			int score = gra->getScoreByMismatch(search_len, wrong_base);
			if(false) fprintf(stderr,"( Simple search used, read region: [%d-%d, len: %d], mismatch %d, score %d) \t", read_skip_left, search_len - read_skip_left, search_len, wrong_base, score);
			return score;
		}
		gra->align_genotyping(true, qseq, qlen, read_in_contig_st_pos, read_in_contig_ed_pos);
		//gra->adjustCIGAR();
		if(false){
			fprintf(stderr,"( Detail search used:");
			gra->printf_alignment_detail(stderr, read_in_contig_st_pos);
			fprintf(stderr,")\t");
		}
		return gra->get_score_reach_end_of_read();
	}

	int get_contig_alignment_score(Genotyping_read_aligner * gra, bam1_t *br, uint8_t * qseq_buff, int read_in_contig_st_pos1, int read_in_contig_st_pos2, int *read_skip_left_, int *read_skip_right_){
		int score = 0;
		if(read_in_contig_st_pos1 == read_in_contig_st_pos2){
			score = get_contig_alignment_score_core(gra, br, qseq_buff, read_in_contig_st_pos1, read_skip_left_, read_skip_right_);
		}else{
			int read_contig_score1 = get_contig_alignment_score_core(gra, br, qseq_buff, read_in_contig_st_pos1, read_skip_left_, read_skip_right_);
			int read_skip_left_bp2, read_skip_right_bp2;
			int read_contig_score2 = get_contig_alignment_score_core(gra, br, qseq_buff, read_in_contig_st_pos2, &read_skip_left_bp2, &read_skip_right_bp2);
			if(read_contig_score1 >= read_contig_score2){
				score = read_contig_score1;
			}else{
				*read_skip_left_ = read_skip_left_bp2; *read_skip_right_ = read_skip_right_bp2;
				score = read_contig_score2;
			}
		}
		if(false)
			fprintf(stderr, "\tread_skip_left %d, read_skip_right % d\t",  *read_skip_left_, *read_skip_right_);
		return score;
	}

	//get the read position of the first base and the last base in the reference
	void get_true_read_pos(bam1_t *br, int * true_read_pos_bg, int *true_read_pos_ed){
		int begin_pos = br->core.pos;
		int end_pos = br->core.pos;

		uint32_t* bam_cigar = bam_get_cigar(br);
		uint32_t n_cigar = br->core.n_cigar;
		for (uint i = 0; i < n_cigar; ++i)
		{
			int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			switch (c_type){
			case CIGAR_MATCH: case CIGAR_SEQ_MATCH: end_pos += c_size; break;
			case CIGAR_INSERT:	break; //do nothing
			case CIGAR_DELETE:	end_pos += c_size; break;
			case CIGAR_SOFT_CLIP:	case CIGAR_HARD_CLIP:
				if(i == 0) 	begin_pos -= c_size; else end_pos += c_size;
				break;
			default:	break;
			}
		}
		*true_read_pos_bg = begin_pos;
		*true_read_pos_ed = end_pos;

		if(false){
			fprintf(stderr,"\t Ori Cigar: ");
			for (unsigned int i = 0; i < n_cigar; ++i){
				int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
				char c_type = segment_type_to_cigar_code(type);
				int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				fprintf(stderr,"%d%c",length,  c_type);
			}
		}
	}

	//this step running before genotyping
	//return false only when distance of two SVs is over 300
	bool duplication_SV_filter(NOVA_SV_FINAL_RST_item & to_compare){
		//position is nearby:
		int ABS_POS = ABS_U(st_pos, to_compare.st_pos);
		if(chr_ID != to_compare.chr_ID || ABS_POS > 300) return false;
		//filter:
		if(is_duplication_sv || to_compare.is_duplication_sv) return true;
		//SV type is same and SV length is similar
		float length_rate = (float)SV_length/to_compare.SV_length;
		if(length_rate < 0.85 || length_rate > 1.15) return true;

		//skip short SVs
		//if(SV_length < 50 && SV_length > -50)		return true;
		//if(to_compare.SV_length < 50 && to_compare.SV_length > -50)		return true;
		//skip BND SVs
		if(is_BND || to_compare.is_BND) return true;

		if(SV_length < 0){//for deletions:
			if(SV_length < to_compare.SV_length)	to_compare.is_duplication_sv = true;
			else									is_duplication_sv = true;
		}else{//insertions
			if(SV_length > to_compare.SV_length)	to_compare.is_duplication_sv = true;
			else									is_duplication_sv = true;
		}
		return true;
	}
	bool SV_is_duplicated(){ return is_duplication_sv; }

	int long_deletion_DR_filter(Bam_file *bam_f, int normal_read_length, int max_isize){

		RefRegion withInR1(chr_ID, st_pos - (max_isize - normal_read_length), st_pos);
		RefRegion withInR2(chr_ID, endPos, endPos + (max_isize - normal_read_length - normal_read_length));
		RefRegion withInDel(chr_ID, st_pos, st_pos + (max_isize - normal_read_length - normal_read_length));

		//only consider reads in the first region:
		R_region bam_load_region;
		bam_load_region.chr_ID = chr_ID;
		bam_load_region.st_pos = st_pos - max_isize;
		bam_load_region.ed_pos = st_pos + normal_read_length;

		//Analyzer
		int read_num_support_deletion = 0;
		int64_t total_isize_deletion = 0;
		int read_num_support_reference = 0;
		int64_t total_isize_reference = 0;


		resetRegion_ID(bam_f, &bam_load_region);
		while (bam_next(bam_f)) {
			bam1_t *br = &(bam_f->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;

			//get iSIZE
			int isize = br->core.isize;
			if(isize <= 0 || isize > 1000000) continue;

			int pos = br->core.pos + br->core.l_qseq;
			int mpos = br->core.mpos;

			if(	withInR1.pos_within_region(pos) && withInR2.pos_within_region(mpos)){//for read pairs that support the deletion
				read_num_support_deletion++;
				total_isize_deletion += isize;
			}

			if(withInR1.pos_within_region(pos) && withInDel.pos_within_region(mpos)){//for read pairs that support the reference
				read_num_support_reference++;
				total_isize_reference += isize;
			}
		}

		fprintf(stderr, "Read_num_support_deletion %d %f\n",
				read_num_support_deletion, (float)total_isize_deletion/read_num_support_deletion);
		fprintf(stderr, "Read_num_support_reference %d %f\n",
				read_num_support_reference, (float)total_isize_reference/read_num_support_reference);

		if((read_num_support_deletion >= read_num_support_reference * 2))			return 2;
		else if((read_num_support_deletion * 3 >= read_num_support_reference))		return 1;
		else																		return 0;
	}

	void genotyping(int normal_read_length, int max_isize, Bam_file *bam_f, Genotyping_read_aligner * gra, RefHandler *refHandler){
		if(is_BND) return;

		//step1: show break point region:
		//load reference:
		region_is_overlap = false;
		//load reference in break point 1:
		int edge_len = normal_read_length;
		RefRegion main(chr_ID, st_pos - 10, st_pos + edge_len);
		RefRegion supp(chr_ID, endPos - 10, endPos + edge_len);
		if(main.region_overlap(supp)){ main.Combine(supp, true); region_is_overlap = true; }
		std::cerr << "Region 1: " << main;
		if(!region_is_overlap) std::cerr << "Region 2: " << supp;
		std::cerr << std::endl;

		uint8_t qseq_buff[1024]; int left_clip = 0;

		fprintf(stderr, "SV breakpoint1 %d SV; breakpoint2 %d\n", st_pos, endPos);

		//suggest contig_region_st for BP1/2
		int contig_pos_bp1; int contig_pos_bp2;
		get_contig_golbal_position(&contig_pos_bp1, &contig_pos_bp2);

		for(int region_ID = 0; region_ID < 3; region_ID++)
			for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)
				region_support_numnber[region_ID][read_genotype_type] = 0;

		for(int region_ID = 0; region_ID < 2; region_ID++){ //mode == 0 for region1, mode == 1 for region2
			if(region_is_overlap == true && region_ID == 1) continue;

			R_region region;
			region.chr_ID = main.chr_ID;
			region.st_pos = ((region_ID == 0)?main.st_pos:supp.st_pos) + 1;
			region.ed_pos = ((region_ID == 0)?main.ed_pos:supp.ed_pos) + 1;

			gra->setRef(refHandler->getRefStr(region_ref_global_position), 100000, contig.c_str(), contig.size());

			fprintf(stderr, "Loading reads for region %d\n", region_ID + 1);

			resetRegion_ID(bam_f, &region);	//reset region
			//reference check:
			while (bam_next(bam_f)) {
				bam1_t *br = &(bam_f->_brec);
				if(bam_is_secondary(br))		continue;
				if(bam_is_supplementary(br))	continue;
//				if(248617642 == br->core.pos){
//					fprintf(stderr, " 0");
//				}

				int overlap_mode = read_overlap_breakpoint(br, region_ID, !region_is_overlap, normal_read_length, &left_clip);
				if(overlap_mode == 0) {
					if(0) fprintf(stderr, "@Read name: %s Overlap_mode %d; SKIP\n", bam_qname(br), overlap_mode);
					continue;
				}
				if(false) fprintf(stderr, "@Read name: %s POS: %d Overlap_mode %d;\t", bam_qname(br), br->core.pos, overlap_mode);

				get_bam_seq_bin(0, br->core.l_qseq, qseq_buff, br);
				if(false){
					for(int i = 0; i < br->core.l_qseq; i++)
						fprintf(stderr, "%c", "ACGT"[qseq_buff[i]]);
				}
				//realigned for contig regions
				int read_contig_score = 0;

				int true_read_pos_bg; int true_read_pos_ed;
				get_true_read_pos(br, &true_read_pos_bg, &true_read_pos_ed);

				if(false) fprintf(stderr, "\ttrue_read_pos: [%d, %d]\t", true_read_pos_bg, true_read_pos_ed);

				int read_in_contig_st_pos_bp1_rbg = (true_read_pos_bg) - contig_pos_bp1;
				int read_in_contig_st_pos_bp2_rbg = (true_read_pos_bg) - contig_pos_bp2;

				int read_in_contig_st_pos_bp1_red = (true_read_pos_ed - br->core.l_qseq) - contig_pos_bp1;
				int read_in_contig_st_pos_bp2_red = (true_read_pos_ed - br->core.l_qseq) - contig_pos_bp2;

				int read_skip_left, read_skip_right;

				if(overlap_mode == 1){
					read_contig_score = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp1_rbg, read_in_contig_st_pos_bp1_red, &read_skip_left, &read_skip_right);
				}else if(overlap_mode == 2){
					read_contig_score = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp2_rbg, read_in_contig_st_pos_bp2_red, &read_skip_left, &read_skip_right);
				}else if(overlap_mode == 3){
					int read_contig_score1 = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp1_rbg, read_in_contig_st_pos_bp1_red, &read_skip_left, &read_skip_right);
					int read_skip_left_bp2, read_skip_right_bp2;
					int read_contig_score2 = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp2_rbg, read_in_contig_st_pos_bp2_red, &read_skip_left_bp2, &read_skip_right_bp2);
					if(read_contig_score1 >= read_contig_score2){
						read_contig_score = read_contig_score1;
					}else{
						read_skip_left = read_skip_left_bp2; read_skip_right = read_skip_right_bp2;
						read_contig_score = read_contig_score2;
					}
				}

				int read_ori_score = get_ori_alignment_score(gra, br, qseq_buff, read_skip_left, read_skip_right, refHandler);

				fprintf(stderr,"Read_contig_score %d, Read_ori_score is %d ", read_contig_score, read_ori_score);
				if(0){
					fprintf(stderr, "@Read string: ");
					for(int i = 0; i < br->core.l_qseq; i++)
						fprintf(stderr, "%c", "ACGT"[qseq_buff[i]]);
				}
				fprintf(stderr, "\n");

				//analysis:
				int min_score = gra->getGenotypingMinScore(br->core.l_qseq - read_skip_left - read_skip_right);
				if(read_contig_score > read_ori_score + 4 && read_contig_score > min_score)		region_support_numnber[region_ID][0]++;
				else if(read_contig_score + 4 < read_ori_score && read_ori_score> min_score)	region_support_numnber[region_ID][1]++;
				else																		region_support_numnber[region_ID][2]++;
			}
		}

		fprintf(stderr,"Region analysis: for region1: ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_numnber[0][0], region_support_numnber[0][1], region_support_numnber[0][2]);
		fprintf(stderr,"Region analysis: for region2: ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_numnber[1][0], region_support_numnber[1][1], region_support_numnber[1][2]);

		for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)
			region_support_numnber[2][read_genotype_type] = region_support_numnber[0][read_genotype_type] + region_support_numnber[1][read_genotype_type];

		//signal read number adjust for insertions and deletions
		if(SV_length < 0)
			region_support_numnber[2][1] /= 1.5;
		if(SV_length > 0)
			region_support_numnber[2][0] /= 1.5;

		if(region_support_numnber[2][0] > region_support_numnber[2][1] * 3){
			final_genotype = 2;//1/1
		}else if(region_support_numnber[2][0] * 3 < region_support_numnber[2][1]){
			final_genotype = 0;//0/0
		}else{
			final_genotype = 1;//0/1
		}

		if((region_support_numnber[2][0] + region_support_numnber[2][1]) * 3 < region_support_numnber[2][2]){
			final_genotype = 0;//0/0
		}

		fprintf(stderr,"Region analysis: for region(1+2): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_numnber[2][0], region_support_numnber[2][1], region_support_numnber[2][2]);
		fprintf(stderr, "final_genotype : %s", (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1"));

		//DR filter:
		deletion_pass_DR_filter = true;
		if(final_genotype > 0 && SV_length < -200){
			fprintf(stderr, "Running long deletion DR filter\n");
			bool pass_dr_filter = long_deletion_DR_filter(bam_f, normal_read_length, max_isize);
			if( 0 == pass_dr_filter){
				deletion_pass_DR_filter = false;
				fprintf(stderr, "Deletion failed the DR filter\n");
			}
		}
		else if(final_genotype == 0 && SV_length < -700){
			fprintf(stderr, "Running long deletion DR filter\n");
			bool pass_dr_filter = long_deletion_DR_filter(bam_f, normal_read_length, max_isize);
			if(pass_dr_filter > 0){
				final_genotype = pass_dr_filter;
				fprintf(stderr, "Long deletion genotype reset using DR signals\n");
			}
		}
	}

	static void write_to_vcf_header(FILE *vcf_output, bam_hdr_t * bam_header){
		//BASIC PART
		time_t c_time; struct tm*p;
		time(&c_time);
		p = gmtime(&c_time);
		//date
		fprintf(vcf_output, "##fileformat=VCFv4.1\n");
		fprintf(vcf_output, "##fileDate=%d%d%d\n",1990+p->tm_year, 1 + p->tm_mon, p->tm_mday);
		fprintf(vcf_output, "##source=%s%s\n",PACKAGE_NAME, PACKAGE_VERSION);//software
		///INFO PART
		fprintf(vcf_output, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"sample_id from dbVar submission; every call must have SAMPLE\">\n");
		fprintf(vcf_output, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
		fprintf(vcf_output, "##INFO=<ID=BND_TYPE,Number=1,Type=String,Description=\"Type of BND variant\">\n");
		fprintf(vcf_output, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
		fprintf(vcf_output, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
		///FORMAT PART
		fprintf(vcf_output, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		fprintf(vcf_output, "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Reads number that supported for the reference, allele and unknown when genotyping\">\n");
		fprintf(vcf_output, "##FORMAT=<ID=SG,Number=.,Type=Integer,Description=\"Suggested SV length by the assembler.\">\n");
		///FILTER PART
		fprintf(vcf_output, "##FILTER=<ID=LOW_DEPTH,Description=\"Variant with genotype [0/0]\">\n");
		fprintf(vcf_output, "##FILTER=<ID=lt50bp,Description=\"Supported variant but smaller than 50bp\">\n");
		fprintf(vcf_output, "##FILTER=<ID=DR_not_support,Description=\"Long deletion (>300) without enough Discordant read pair signal support\">\n");
		fprintf(vcf_output, "##FILTER=<ID=Duplication_sv,Description=\"When it is true, this SV is similar with others, and will be removed\">\n");

		//ALT PART
		fprintf(vcf_output, "##ALT=<ID=DEL,Description=\"Deletion\">\n");
		fprintf(vcf_output, "##ALT=<ID=INS,Description=\"Insertion\">\n");
		fprintf(vcf_output, "##ALT=<ID=TRA,Description=\"Trans-location\">\n");
		fprintf(vcf_output, "##ALT=<ID=INV_1,Description=\"Inversion type1\">\n");
		fprintf(vcf_output, "##ALT=<ID=INV_2,Description=\"Inversion type2\">\n");

		//CONTIG PART
		for(int i = 0; i < bam_header->n_targets; i++)
			fprintf(vcf_output, "##contig=<ID=%s,length=%d>\n", bam_header->target_name[i], bam_header->target_len[i]);
		fprintf(vcf_output, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample\n");
	}

	bool will_be_output_to_vcf;

private:
	//basic
	int chr_ID;
	int st_pos;
	std::string SV_type;
	std::string SV_name;

	//variation information
	int SV_length;
	int endPos;
	std::string ref;
	std::string alt;
	std::string bnd_type;

	///supplementary informations///

	//suggest SV length, those data are generated from the analysis of positions of assembled reads
	int suggenst_SV_length;//the suggested SV length given by the assembly problems
	int dis_SV_len_suggset_length;//distance between SV length and suggest SV length

	//assembly results
	std::string contig;//the assembly results this is derived from
	std::vector<uint32_t> cigar; //the CIGAR of contig aligned to the reference
	int cigar_idx_bg;//the cigar index of this SV
	int cigar_idx_ed;//the cigar index of this SV
	int region_ref_global_position;//the global position of reference region, a reference region is the region collected signals and calling SVs in the SV calling steps
	int contig_st_pos_in_ref;//the position of contig in the reference region

	//genotyping results
	bool region_is_overlap;//region around bp1 is overlapped with region around bp2
	//region include: [region1, region2 and region(1+2)]; read type include: [support ALT, support REF, and (not) support both]
	int region_support_numnber[3][3]; //[region number] * [support read type number]
	//final_genotype = 0: 0/0; final_genotype = 1: 0/1; final_genotype = 2: 1/1
	int final_genotype;
	bool deletion_pass_DR_filter;
	bool is_duplication_sv;//when it is true, this SV is similar with others, and will be removed

	//
	bool is_BND;
};

#endif /* NOVASVGENERATEVCF_NOVASVRST_HPP_ */
