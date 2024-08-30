/*
 * SignalAssembly.hpp
 *
 *  Created on: 2021年4月20日
 *      Author: fenghe
 */

#ifndef SIGNALASSEMBLY_HPP_
#define SIGNALASSEMBLY_HPP_

#include <vector>
#include <algorithm>
extern "C"
{

#include "../kswlib/kalloc.h"
#include "../kswlib/ksw2.h"

#include "../clib/utils.h"
#include "../clib/bam_file.h"
#include "../clib/vcf_file.h"
#include "../clib/desc.h"
}

#include "../cpp_lib/Assembler/mantaAssembler.hpp"
#include "../cpp_lib/get_option_cpp.hpp"
#include "deBGA_index.hpp"

#include "../PanSVgenerateVCF/generateVCFoptions.hpp"
#include "../PanSVgenerateVCF/signalSAMLoader.hpp"
#include "../PanSVgenerateVCF/SV_ref_sequence.hpp"

struct GlobalDepthItem{
	int16_t ACGTD_num[6];//4 for a deletion ,5 for insertion
	int16_t ACGTD_num_tmp[6];//4 for a deletion ,5 for insertion at next base
	int16_t total_depth;
	int16_t cur_ab_block;
	uint8_t ref_base;
	uint8_t max_base;
	void set_base(uint8_t base, int ab_block, int depth){
		if(cur_ab_block == ab_block){
			ACGTD_num_tmp[base] = MAX(ACGTD_num_tmp[base], depth);
		}else{
			cur_ab_block = ab_block;
			ACGTD_num[base] += ACGTD_num_tmp[base];
			ACGTD_num_tmp[base] = depth;
		}
	}

	void add_tmp_data(){
		int max_base_count = 0;
		for(int base = 0; base < 6; base++){
			ACGTD_num[base] += ACGTD_num_tmp[base];
			total_depth += ACGTD_num[base];
			if(ACGTD_num[base] > max_base_count){
				max_base = base; max_base_count = ACGTD_num[base];
			}
		}
		//ref_base = ref_base_;
	}

	int event_info(){
		if(total_depth == 0)
			return 1;
		else if(max_base != ref_base)//3~8
			return 3 + max_base;
		else if(ACGTD_num[max_base] != total_depth){
			return 2;
		}
		return 0;
	}

	void print(FILE * output){
		fprintf(output, "[ Global depth: [ref %c, MAX %c]"
				" [depth: ",
				"ACGT-I"[ref_base],
				(ref_base == max_base)?'=':"ACGT-I"[max_base]);
		for(uint8_t base = 0; base < 6; base++){
			int base_depth = ACGTD_num[base];
			if(base_depth != 0) fprintf(output, "%c:%d ", "ACGT-I"[base], base_depth);
		}
		fprintf(output, " ]]\t");
	}

	static void print_full_data(FILE * output, std::vector<GlobalDepthItem> & di, int sv_length, int break_point1, int break_point2){
		//detail around the break point
		print_data(output, di, 0, sv_length);
		fprintf(output, "\n");
		bool has_ins_part = (break_point2 > break_point1 + 10);
		int bp_chech_region_len;
		if(has_ins_part) bp_chech_region_len = 10; else bp_chech_region_len = 20; //10 for insertion; 20 for deletion or other type
		print_data(output, di, break_point1 - bp_chech_region_len , break_point1 + bp_chech_region_len);
		fprintf(output, "Above is break point 1, [ %d %d ].\n", break_point1 - bp_chech_region_len , break_point1 + bp_chech_region_len);

		print_data(output, di, break_point2 - bp_chech_region_len , break_point2 + bp_chech_region_len);
		fprintf(output, "Above is break point 2, [ %d %d ].\n", break_point2 - bp_chech_region_len , break_point2 + bp_chech_region_len);

		if(break_point2 > break_point1 + 10){
			print_data(output, di, break_point1, break_point2 );
			fprintf(output, "Above is INS part, [ %d %d ] .\n", break_point1 , break_point2);
		}
		fprintf(output, "\n");
	}

private:
	static void print_data(FILE* output, std::vector<GlobalDepthItem> & di, int st_pos, int ed_pos){
			//for reference string
			for(int i = st_pos; i < ed_pos; i++)
				fprintf(output, "%c", "ACGT"[di[i].ref_base]);
			fprintf(output, "\n");
			//for action analysis string
			for(int i = st_pos; i < ed_pos; i++){
				int ei = di[i].event_info();
				if(ei == 0) fprintf(output, " ");
				else		fprintf(output, "%d", ei);
			}
			fprintf(output, "\n");
			//for base information string
			for(uint8_t base = 0; base < 6; base++){
				for(int i = st_pos; i < ed_pos; i++){
					int base_depth = di[i].ACGTD_num[base];
					if(base_depth == 0) fprintf(output, " ");
					else				fprintf(output, "%c", '#' + base_depth);
				}
				fprintf(output, "\n");
			}
		}
};

//'ALL' 'DEL' 'INS' 'DUP' 'TRA' 'INV' 'BND'");
struct VAR_TYPE{
	static const int SNP = 0;
	static const int DEL = 1;
	static const int INS = 2;
	static const int UNKNOWN = 4;
};

struct VariationInfo{
	VariationInfo(const VariationInfo & cp){
		ref = cp.ref;
		alt = cp.alt;
		set_basic(cp.ref_position, cp.alt_position, cp.depth, cp.assembly_part, cp.contigID, cp.var_type);
	}

	VariationInfo(	std::string ref_, std::string alt_, int ref_position_, int alt_position_, int depth_, int assembly_part_, int contigID_, int var_type_){
		ref = ref_;
		alt = alt_;
		set_basic(ref_position_, alt_position_, depth_, assembly_part_, contigID_, var_type_);
	}

	VariationInfo(int ref_position_, int alt_position_, int depth_, int assembly_part_, int contigID_, int var_type_){
		set_basic(ref_position_, alt_position_, depth_, assembly_part_, contigID_, var_type_);
	}

	void set_basic(int ref_position_, int alt_position_, int depth_, int assembly_part_, int contigID_, int var_type_){
		ref_position = ref_position_;
		alt_position = alt_position_;
		depth = depth_;
		assembly_part = assembly_part_;
		contigID = contigID_;
		var_type = var_type_;
		pass_filter = false;
	}

	bool try_merge(const VariationInfo &try_item){
		if(same_var(try_item)){
			if(assembly_part != try_item.assembly_part){
				depth += try_item.depth;
				assembly_part = try_item.assembly_part;
			}else
				depth = MAX(depth, try_item.depth);
			return true;
		}
		else
			return false;
	}

	bool same_var(const VariationInfo &try_item){
		if(ref_position == try_item.ref_position && var_type == try_item.var_type && ref == try_item.ref && alt == try_item.alt)
			return true;
		return false;
	}

	static inline int cmp_by_pos(const VariationInfo &a, const VariationInfo &b){
		//var basic
		if(a.ref_position != b.ref_position) 	return a.ref_position < b.ref_position;
		if(a.var_type != b.var_type)			return a.var_type < b.var_type;
		if(a.ref != b.ref)						return a.ref < b.ref;
		if(a.alt != b.alt)						return a.alt < b.alt;

		//contig info
		if(a.assembly_part != b.assembly_part)	return a.assembly_part < b.assembly_part;
		if(a.contigID != b.contigID)			return a.contigID < b.contigID;

		return a.depth > b.depth;
	}

	void print(FILE * output,  std::vector<GlobalDepthItem> &global_depth){
		fprintf(output, "[ref: %s @ %d; alt: %s @ %d; depth: %d, "
						"[%d %d %d]] %s\t",
						ref.c_str(), ref_position, alt.c_str(), alt_position, depth,
						assembly_part, contigID, var_type, pass_filter?"PASS":"FAIL");
		global_depth[ref_position].print(output);
		fprintf(output, "\n");
	}

	std::string ref;
	std::string alt;
	int ref_position;
	int alt_position;
	int depth;
	int assembly_part;//from which assembly part
	int contigID;//from which contig
	int var_type;
	int pass_filter;
	//cigar
};

const char* ACGT = "ACGT";

struct VI_list{
	std::vector<VariationInfo> vi;
	std::vector<VariationInfo> vi_tmp;

	void clear(){vi.clear();}
	void add_data(char * ref, int ref_len, char * alt, int alt_len, int ref_position_, int alt_position_, int depth_, int assembly_part_, int contigID_, int var_type_){
		vi.emplace_back(ref_position_, alt_position_, depth_, assembly_part_, contigID_, var_type_);
		vi.back().ref.append(ref, ref_len);
		vi.back().alt.append(alt, alt_len);
	}
	void add_data(uint8_t * ref, int ref_len, uint8_t * alt, int alt_len, int ref_position_, int alt_position_, int depth_, int assembly_part_, int contigID_, int var_type_){
		vi.emplace_back(ref_position_, alt_position_, depth_, assembly_part_, contigID_, var_type_);
		std::string &c_ref = vi.back().ref; for(int i = 0; i < ref_len; i++) c_ref += ACGT[ref[i]];
		std::string &c_alt = vi.back().alt; for(int i = 0; i < alt_len; i++) c_alt += ACGT[alt[i]];
	}
	void sort_merge(FILE * output, std::vector<GlobalDepthItem> &global_depth){
		std::sort(vi.begin(), vi.end(), VariationInfo::cmp_by_pos);

		if(false) for(auto & c_vi: vi)	c_vi.print(output, global_depth);

		//merge and run filter
		vi_tmp.clear();
		auto vi_ed = vi.end();
		for(auto c_vi = vi.begin(); c_vi < vi_ed;)
			for(auto try_item = c_vi + 1; try_item <= vi_ed; try_item++)
				if(try_item == vi_ed || (c_vi->try_merge(try_item[0]) == false)){ 	vi_tmp.emplace_back(c_vi[0]); 	c_vi = try_item; 	break;	}
		std::swap(vi_tmp, vi);

		if(false){	fprintf(output, "After merge: \n");	for(auto & c_vi: vi) c_vi.print(output, global_depth);	}
		//simple filter
		for(auto & c_vi: vi)
			if(c_vi.depth * 4 >= global_depth[c_vi.ref_position].total_depth && c_vi.depth > 2) c_vi.pass_filter = true;

		fprintf(output, "After simple filter: \n");
		for(auto & c_vi: vi)	if(c_vi.pass_filter) c_vi.print(output, global_depth);
	}
};

struct Contig_aligner{

//#define KSW_ALN_left_extend_forward_cigar 3
	SV_ref_sequence *sf;
	uint8_t *tseq;
	int tlen;

	//query temp: used for analysis
	uint8_t *qseq;
	uint32_t qlen;

	//read info
	int32_t read_score;
	uint32_t total_q_len;

	//sequence buff
	uint8_t type;

	//result
	bool is_simple_aln;
	uint32_t simple_NM;
	ksw_extz_t ez;

	//kmalloc
	void *km;

	std::vector<uint8_t> ref_buff;

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
	deBGA_INDEX * idx;

	//as input
	void printf_alignment_detail(FILE * output, int suggest_st_pos, uint16_t * contig_depth, int ref_pos_enough_match_base){

		fprintf(output, "\n");
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0){
			for(int i = 0; i < (int)qlen; i++) {fprintf(output, "%c", "ACGT"[qseq[i]]);}
			fprintf(output, "\n ???NO alignment result\n");
			return;
		}
		//print contig sequence
		if(true){
			int output_index = 0;
			int seq_i = 0;
			for(int i = 0; i < suggest_st_pos; i++) {fprintf(output, "-");  output_index++;}
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
				case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "%c", "ACGT"[qseq[seq_i]]); output_index++;}	break;//M
				case 1:	seq_i += cigar_len;	break;//I, print nothing
				case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(output, "-");  output_index++;}	break;//D, print '-'
				case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "N"); output_index++;}	break;//N, print N
				case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "-"); output_index++;}	break;//S, print -
				default: fprintf(output, "ERROR CIGAR  %d %d ", type, cigar_len);
				}
			}
			fprintf(output, "\n");
		}

		//info:
		int contig_coverage = 0;
		int number_mismatch = 0;
		int coverage_increase = 0;

		//print X/= sequence
		if(true){
			int output_index = 0;
			int seq_i = 0;
			for(int i = 0; i < suggest_st_pos; i++) {fprintf(output, "-");  output_index++;}
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
				case 0:
					for(int i = 0; i < cigar_len; i++, seq_i++)
					{
						if(qseq[seq_i] == tseq[output_index]){
							contig_coverage++;
							if(output_index < ref_pos_enough_match_base) fprintf(output, "M");//1 == M; 2 == X; 0 == -;
							else fprintf(output, "=");//1 == M; 2 == X; 0 == -;
						}
						else{
							number_mismatch++;
							fprintf(output, "X");//1 == M; 2 == X; 0 == -;
						}
						output_index++;
					}	break;//M
				case 1:	seq_i += cigar_len;	break;//I, print nothing
				case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(output, "-");  output_index++;}	break;//D, print '-'
				case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "N"); output_index++;}	break;//N, print N
				case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "-"); output_index++;}	break;//S, print -
				default: fprintf(output, "ERROR CIGAR  %d %d ", type, cigar_len);
				}
			}
			fprintf(output, "\n");
		}

		//print coverage sequence
		if(true){
			int output_index = 0;
			int seq_i = 0;
			for(int i = 0; i < suggest_st_pos; i++) {fprintf(output, "-");  output_index++;}
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
				case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "%c", contig_depth[seq_i] + '#'); output_index++;}	break;//M
				case 1:	seq_i += cigar_len;	break;//I, print nothing
				case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(output, "-");  output_index++;}	break;//D, print '-'
				case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "0"); output_index++;}	break;//N, print N
				case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(output, "-"); output_index++;}	break;//S, print -
				default: fprintf(output, "ERROR CIGAR  %d %d ", type, cigar_len);
				}
			}
			fprintf(output, "\t");
		}

		if(true){
			fprintf(output, "\nCigar sequence: ");
			for(int i = 0; i < cigar_len; i++){
				fprintf(output, "%d%c", (bam_cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(bam_cigar[i] & BAM_CIGAR_MASK)]);
			}
			fprintf(output, "contig_coverage: [%d %f]\tnumber_mismatch[%d]\t", contig_coverage, (double)contig_coverage/tlen, number_mismatch);
			fprintf(output, "coverage_increase: [%d %f]\t", coverage_increase, (double)coverage_increase/tlen);

			fprintf(output, "\t");
		}
		fprintf(output, "\n");
	}

	void printf_ref(FILE * output){
		fprintf(output, "Reference sequence: len %d\n", tlen);
		for(int i = 0; i < tlen; i++){
			fprintf(output, "%c", "ACGT"[tseq[i]]);
		}
		fprintf(output, "\n");
	}

	void copy_option(){
			match_D = 2;
			mismatch_D= 10;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 1;
			zdrop_D= gap_open2_D + 100; //for DNA zdrop = 400, 200 for RNA
			bandwith = zdrop_D;
			flag = 0;
	}

	void init(SV_ref_sequence *sf_){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option();
		ksw_gen_mat_D();
		sf = sf_;
		//index
		//idx = idx_;
	}

	void setRef(int sv_id){
		tseq = sf->get_sv_str_by_id(sv_id, tlen);
	}

	void ksw_gen_mat_D(){
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

	void free_memory(){
		if(ez.cigar)   {kfree(km, (void *)(ez.cigar)); }
		km_destroy(km);
	}

	void align_non_splice(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		ksw_extd2_sse(km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

};

struct ASS_add_reads{
	ASS_add_reads(bam1_t  * br_, int read_in_ref_offset_, bool is_from_main_SV_){
		br = br_; read_in_ref_offset = read_in_ref_offset_; is_from_main_SV = is_from_main_SV_;
		int new_score = 0; int ori_score = 0;
		bam_get_num_tag(br, vcf_tag.AS, &new_score);
		bam_get_num_tag(br, vcf_tag.OS, &ori_score);
		char* OA_tag_char = bam_get_string_tag(br, vcf_tag.OA);
		char UM = OA_tag_char[strlen(OA_tag_char) - 2];
		xassert(UM == 'U' || UM == 'M', "");
		new_score_is_bigger = (UM == 'U' || new_score > ori_score);
		new_score_is_far_bigger= (UM == 'U' || new_score > ori_score + 40 || new_score > ori_score * 2);
	}
	bam1_t  * br;
	int read_in_ref_offset;
	bool is_from_main_SV;
	bool new_score_is_bigger;
	bool new_score_is_far_bigger;
};

struct AssemblyBlock{
	//input read info
	std::vector<ASS_add_reads> ass_read_list;
	//input read data
	std::vector<std::string> reads;
	//output: read result
	std::vector<AssembledContig> contigs;

	void init(){
		ass_read_list.clear();
		reads.clear();
	}

	void add_read(bam1_t  * br, char * read_str, int sv_region_st, bool is_from_main_SV){
		ass_read_list.emplace_back(br, br->core.pos - sv_region_st, is_from_main_SV);
		reads.emplace_back(read_str);
	}

	void run_assembly(AssemblyManager *am){
		std::swap(am->reads, reads);
		am->assembley();
		std::swap(am->contigs, contigs);
		std::swap(am->reads, reads);//swap back the read list
	}
};

struct READ_ST_IDX{
	READ_ST_IDX(int sv_ID_, int load_read_number_, assembly_sorter_item * index_st_, bam1_t  * read_buff_, bool is_main_){
		sv_ID = sv_ID_;
		load_read_number= load_read_number_;
		index_st = index_st_;
		read_buff = read_buff_;
		is_main = is_main_;
	}

	int sv_ID;
	int load_read_number;
	assembly_sorter_item * index_st;
	bam1_t  * read_buff;
	bool is_main;
};

struct AssemblyHandler{
	bool output_original_read;
	bool output_assembly_contig;

	void init(SignalAssemblyOption *o_, SV_ref_sequence *sf_){
		o = o_;//get option
		FILE * output_file = NULL;
		if(strcmp(o->output_vcf, "stdout") == 0) output_file = stdout;
		else output_file = xopen(o->output_vcf, "w");
		summary_output = output_file;
		if(!o->not_output_vcf_header)
			sf_->output_VCF_header(summary_output);

		detail_output = stderr;
		output_original_read = false;
		output_assembly_contig = false;

		if(o->print_detail_information){
			output_original_read = true;
			output_assembly_contig = true;
		}

		depth_counter = (uint32_t *)xmalloc(100000 * sizeof(uint32_t));
		read_str = (char *)xmalloc(10000);
		am = new AssemblyManager[1];
		ca.init(sf_);

		normal_read_len = o->normal_read_length;
		ave_read_depth = o->ave_read_depth;
		max_read_depth = o->ave_read_depth * 2;
		rsf_block_size = 32;
		rsf_shift_offset = 5;
		ab_block_size_option = 300;
		max_read_in_block = max_read_depth*rsf_block_size / normal_read_len;
	}

	void destroy(){
		if(summary_output != NULL && summary_output != stdout){ fclose(summary_output); }
		free(read_str);
		delete[] am;
		free(depth_counter);
	}

	//int read_register(Remapped_read_loader &rl);
	int read_register(Remapped_read_loader &rl, SV_ref_sequence &sf);
	void run_assembly_analysis(){
		if(SV_already_used) return;
		addReads2assembler();
//		fprintf(detail_output, "detail of read depth");
//		for(int i = 0; i < sv_length; i++){
//			fprintf(detail_output, "index: %d read depth: %d\n", i, depth_counter[i]);
//		}
		if(o->show_coverage_depth_analysis) coverage_analysis();
		getContig();
	}

	std::vector<uint8_t> bin_contig;

	//	Contig_aligner
	Contig_aligner ca;

private:

	//basic I/O and options
	FILE * summary_output;
	FILE * detail_output;
	SignalAssemblyOption *o;
	int normal_read_len;

	//*********************Basic SV information********************************//
	//sv info
	int sv_id;
	std::vector<READ_ST_IDX> cluster_used_sv_id;//store the SV id list of used SVs in this cluster
	int sv_length;
	int sv_region_st;//sv region start position in reference
	char* SV_tag_char;
	SV_chr_info * sv_info;
	bool SV_already_used;
	//original position filter
	int ori_accept_region1_bg;
	int ori_accept_region1_ed;
	int ori_accept_region2_bg;
	int ori_accept_region2_ed;

	//*********************Read loader ********************************//
	//read info
	int start_index;
	int load_read_number;
	assembly_sorter_item * index_st;
	bam1_t  * read_buff;
	char * read_str;//read data buff

	//**************Depth counter and score filter*****************//
	//global options
	double ave_read_depth;
	double max_read_depth; // == 1.5 * ave_read_depth
	int max_read_in_block;
	int rsf_block_size;
	int rsf_shift_offset;
	//depth analysis and score filter
	//coverage analysis Variable
	uint32_t * depth_counter;
	//std::vector<uint32_t> max_score;
	std::vector<uint8_t> read_filter;
	std::map<int, int> score_counter;//key is score, second is the count of score
	struct read_score_filter{
		int read_idx_bg;
		int read_idx_ed;
		int filter_score;
		int passed_read_num;
	};
	std::vector<read_score_filter> rsf_buff;
	int total_filter_block_num;
	read_score_filter * get_rsf_by_pos(int pos);

	//**************Assembly blocks*****************//
	//assembly manager
	AssemblyManager *am;
	//contig information
	std::vector<uint16_t> contig_depth;
	std::set<int> remove_read_set;
	int ab_block_size_option;
	int ab_block_size;
	std::vector<AssemblyBlock> ab_list;
	int ab_size;

	std::map<int, int > suggest_st_pos_map;//the suggested alignment positions for contig to be aligned at original reference and the read number that support that position
	std::vector<int> suggent_pos_list;//the suggested alignment positions for contig to be aligned at original reference

	//coverage analysis
	//**************Variation informations*****************//
	VI_list vi_l;
	std::vector<GlobalDepthItem> global_depth;//depth used for all sv

	bool get_read_filter1(bam1_t  * br);// a position filter
	bool read_depth_filter(bam1_t  * br, bool read_from_main_SV, int read_idx);
	void addRead2depthCounter(bam1_t  * br);
	void output_reads(bam1_t  * br, int SCORE_FILTER);
	void addReads2assembler();

	void show_assembly_block_info(int ab_idx ,size_t contig_size);
	bool get_suggention_alignment_position_list(AssembledContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_add_reads> &ass_read_list );
	void contig_alignment_and_get_var(int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len);
	void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);
	void getContig();
	void print_seq(bam1_t  * br);
	void print_info(bam1_t  * br );
	float averageDepth(int begin_pos, int end_pos);
	void coverage_analysis();
};

#endif /* SIGNALASSEMBLY_HPP_ */
