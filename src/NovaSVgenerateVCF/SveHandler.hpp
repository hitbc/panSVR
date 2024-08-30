/*
 * SveHandler.hpp
 *
 *  Created on: 2020年4月29日
 *      Author: fenghe
 */

#ifndef SRC_SIGNAL_SVEHANDLER_HPP_
#define SRC_SIGNAL_SVEHANDLER_HPP_

#include "../cpp_lib/statistics/StatsManager.hpp"
#include "../cpp_lib/Assembler/assembler.hpp"
#include "../NovaSVgenerateVCF/ReadHandler.hpp"
#include "../NovaSVgenerateVCF/RefHandler.hpp"
#include "NovaSVRst.hpp"

extern "C"{
	extern int vcf_write_line(htsFile *fp, kstring_t *line);
	#include "../kswlib/kalloc.h"
	#include "../kswlib/ksw2.h"
	#include "../clib/desc.h"
}

#include <array>
#include <map>

struct BAM_STATUS_NOVA_SV{
	//global analysis variables
	//analysis results before getting all signals by sampling
	uint32_t minInsertLen;
	uint32_t middleInsertLen;
	uint32_t maxInsertLen;
	std::vector<float> isize_distribution;//from minInsertLen to maxInsertLen
	int analysis_read_length;
	double ave_read_len;
	//analysis results after all data is processed
	double ave_read_depth;

	void load(char * read_status_fn){
		FILE * status_file = xopen(read_status_fn, "r");
		uint32_t minInsertLen_l2;
		uint32_t maxInsertLen_l2;
		int fs_rst = fscanf(status_file, "%lf_%d_%d_%d_%d_%d\n", &ave_read_depth, &analysis_read_length, &minInsertLen_l2, &maxInsertLen_l2, &minInsertLen, &maxInsertLen);
		if(fs_rst == 0) fprintf(stderr, "fs_rst error");
		for(uint32_t i = minInsertLen; i < maxInsertLen; i++){
			float c_isize;
			fs_rst = fscanf(status_file, "%f\n", &c_isize);
			if(fs_rst == 0) fprintf(stderr, "fs_rst error");
			isize_distribution.emplace_back(c_isize);
		}
		fprintf(stderr, "BAM/CRAM status: read length: [Normal: %d] ISIZE: [MIN: %d MAX: %d]\n", analysis_read_length, minInsertLen, maxInsertLen);

		int isize_idx = minInsertLen;
		fprintf(stderr, "ISIZE distribution\n");
		for(float & isize_rate: isize_distribution)
			fprintf(stderr, "ISIZE %d, rate %f%%\n", isize_idx++, isize_rate*100);

		fclose(status_file);
	}

	void generate_read_length(const char * referenceFilename, const char * bamFile){
		//S1 :
		int MAX_ANA_READ_LEN = 1000; //1000
		Bam_file c_b;
		memset(&c_b, 0, sizeof(Bam_file));
		bam_file_open(bamFile, referenceFilename, NULL, &c_b);
		bam_hdr_t* hdr = c_b._hdr;
		int total_read_number = 0;
		bam1_t b1 = {0};//BAM record for the first read in a pair
		int sam_rst1 = 0;
		uint64_t *read_length_analysis = (uint64_t *)xcalloc(MAX_ANA_READ_LEN, sizeof(uint64_t));
		while (1){
			//load SAM 1 & 2
			do{	sam_rst1 = sam_read1(c_b._hfp, hdr, &b1); } while( (sam_rst1 >= 0) && (bam_is_secondary(&b1) || bam_is_supplementary(&b1)));
			if(sam_rst1 < 0)		break;
			total_read_number ++;
			if(total_read_number == 100000)
				break;
			//global analysis:
			if(b1.core.l_qseq < MAX_ANA_READ_LEN)
				read_length_analysis[b1.core.l_qseq]++;
		}
		bam_file_close(&c_b);

		//output analysis results
		//part1: get normal read length
		analysis_read_length = -1;
		double total_read_len = 0;
		for(int i = 0; i < MAX_ANA_READ_LEN; i++){
			total_read_len += i*read_length_analysis[i];
			if(read_length_analysis[i] > 0.6*total_read_number){
				analysis_read_length = i; break;
			}
		}
		ave_read_len = total_read_len/total_read_number;
		if(analysis_read_length == -1)
			analysis_read_length = ave_read_len;
		free(read_length_analysis);
	}

	void generate_isize_distribution(const char * referenceFilename, const char * bamFile){
		StatsManager rstats(referenceFilename, "");
		rstats.handleBamCramStats(bamFile, &ave_read_depth);
		minInsertLen = rstats.getInsertLen(bamFile, 0.01f);
		middleInsertLen = rstats.getInsertLen(bamFile, 0.5f);
		maxInsertLen = rstats.getInsertLen(bamFile, 0.99f);

		unsigned idx = rstats.getGroupIndex(StatLabel(bamFile, ""));
		int totalPairedReadCount = rstats.getStats(idx).readCounter.totalHighConfidenceReadPairCount() + 1;
		const SizeDistribution &fragStats = rstats.getStats(idx).fragStats;
		for(uint32_t i = minInsertLen; i < maxInsertLen; i++){
			int index_count = fragStats.getSizeCount(i);
			isize_distribution.emplace_back((float)index_count/totalPairedReadCount);
		}
	}

	void generate_state(const char * referenceFilename, const char * bamFile){
		fprintf(stderr, "Generating BAM status...\n");
		generate_read_length(referenceFilename, bamFile);
		//get stats: M2
		generate_isize_distribution(referenceFilename, bamFile);
		if(true)
			fprintf(stderr, "BAM/CRAM status: read length: [Normal: %d, AVE: %f] ISIZE: [MIN: %d MIDDLE:%d MAX: %d] ave_read_depth [%f]\n",
				analysis_read_length, ave_read_len,
				minInsertLen, middleInsertLen, maxInsertLen, ave_read_depth);
		if(true){
			int isize_idx = minInsertLen;
			fprintf(stderr, "ISIZE distribution\n");
			for(float & isize_rate: isize_distribution)
				fprintf(stderr, "ISIZE %d, rate %f%%\n", isize_idx++, isize_rate*100);
		}
	}

	void getBreakPoint_Distribution( std::vector<float> &DR_bp_distribution, std::vector<float> &SH_bp_distribution,
			std::vector<float> &UM_stPos_distribution, int &START_OFFSET_UM){
		//for discordant /hard clip signals break point distribution
		int min_probability_size = minInsertLen -  2*analysis_read_length; if(min_probability_size < 1) min_probability_size = 1;
		int max_probability_size = maxInsertLen -  2*analysis_read_length;
		DR_bp_distribution.clear();
		DR_bp_distribution.resize(max_probability_size, 0);

		//for DR signals
		for(int i = min_probability_size; i < max_probability_size; i++){
			float PI_per_position = ( isize_distribution[i + 2*analysis_read_length - minInsertLen])/i;
			for(int j = 0; j < i; j++ ) DR_bp_distribution[j] += PI_per_position;
		}
		float sum_p = 0; for(unsigned int i = 0; i < DR_bp_distribution.size(); i ++) sum_p += DR_bp_distribution[i];
		float levelup_rate = 1/sum_p; for(unsigned int i = 0; i < DR_bp_distribution.size(); i ++) DR_bp_distribution[i] *= levelup_rate; //sum will be 1

		//for SR signals
		//for Soft/hard clip signals break point distribution
		int max_probability_size_SH = 10;
		SH_bp_distribution.resize(max_probability_size_SH, 0.0);
		float SH_bp_distribution_sum = 0;
		for(int i = 0; i < max_probability_size_SH; i++){
			SH_bp_distribution[i] = (max_probability_size_SH - i) * (max_probability_size_SH - i);
			SH_bp_distribution_sum += SH_bp_distribution[i];
		}
		for(int i = 0; i < max_probability_size_SH; i++) SH_bp_distribution[i] /= SH_bp_distribution_sum;
		//for UM start position distribution:
		START_OFFSET_UM = minInsertLen -  analysis_read_length;
		UM_stPos_distribution.resize(maxInsertLen - minInsertLen, 0);
		for(uint32_t i = 0; i < maxInsertLen - minInsertLen; i++)
			UM_stPos_distribution[i] = isize_distribution[i];

		//debug show distribution:
		if(0){
			fprintf(stderr, "DR_bp_distribution\n");
			for(float d: DR_bp_distribution)
				fprintf(stderr, "%f\n", d);

			fprintf(stderr, "SH_bp_distribution\n");
			for(float d: SH_bp_distribution)
				fprintf(stderr, "%f\n", d);

			fprintf(stderr, "UM_stPos_distribution\n");
			for(float d: UM_stPos_distribution)
				fprintf(stderr, "%f\n", d);
		}
	}
};

struct SIGNAL_PARAMETER{
	//for DR signal
	int insert_size_min; //0.01%
	int insert_size_max;//0.99%
	int insert_size_combine;    //0.80%
	int insert_region_len;
	int MaxReadLen;
	double read_depth;

	int max_del_dup_length;//the max length of deletion or duplication, when deletion or duplication longer than this value, variations will be treated as BND

	int SVE_MIN_SOLID_SCORE;
	int SVE_MIN_READ_NUM;

	int SVE_combine_min_score_step1;

	void init( BAM_STATUS_NOVA_SV *bs){
		insert_size_max = bs->maxInsertLen;
		insert_size_min = bs->minInsertLen; if(insert_size_min < 0) insert_size_min = 0;
		insert_size_combine = bs->maxInsertLen;
		MaxReadLen = bs->analysis_read_length;
		insert_region_len = insert_size_max - MaxReadLen*2;
		read_depth = bs->ave_read_depth;

		SVE_MIN_SOLID_SCORE = bs->ave_read_depth * 0.08;// 40X ----> 3; 60X -----> 4; 29X -------> 2;
		SVE_MIN_READ_NUM = bs->ave_read_depth * 0.1;//40X ---> 4; 60X ------> 6; 29X ------> 2;
		SVE_MIN_SOLID_SCORE = MAX(2, SVE_MIN_SOLID_SCORE);//at lease 2
		SVE_MIN_READ_NUM = MAX(3, SVE_MIN_READ_NUM);//at least 3

		SVE_combine_min_score_step1 = bs->ave_read_depth * 0.4; //66-> 26; 44 -> 17; 30 -> 12
		SVE_combine_min_score_step1 = MIN(SVE_combine_min_score_step1, 29);
		SVE_combine_min_score_step1 = MAX(SVE_combine_min_score_step1, 16);

		max_del_dup_length = 50000;
	}
};

struct ORI_REF_Depth_Item{
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

};

struct SUM_variations{
	int ref_idx;
	int read_idx;

	int total_var_number;
	int total_var_size;
	char ref_alt_char;//

	int cigar_index_bg;
	int cigar_index_ed;

	void clear(){
		memset(this, 0, sizeof(SUM_variations));
		ref_idx = -1;
	}
	void store_signal(int ref_idx_, int read_idx_, int c_length, char ref_alt_char_, int cigar_idx){
		if(ref_idx == -1) {
			ref_idx = ref_idx_;
			read_idx = read_idx_;
			ref_alt_char = ref_alt_char_;
			cigar_index_bg = cigar_idx;
		}
		total_var_number ++;
		total_var_size += c_length;
		cigar_index_ed = cigar_idx;
	}

	bool minSV_filter(int minSVlen){
		if(total_var_number > 1 && ref_idx != -1 && total_var_size >= minSVlen)
			return true;
		return false;
	}

};

struct Contig_String_aligner{

public:
	//as input
	void init(){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option();
		ksw_gen_mat_D();
	}

	void setRef(uint8_t * ref_string, int ref_len, 	int chr_ID_, int global_ref_pos_){ tseq = ref_string; tlen = ref_len; chr_ID = chr_ID_;global_ref_pos = global_ref_pos_;}

	void destory(){
		if(ez.cigar != NULL)
			free(ez.cigar);
		km_destroy(km);
	}

	void align_non_splice(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		last_ref_st_pos = ref_st_pos;
		last_ref_end_pos = ref_end_pos;

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

	int adjustCIGAR(){
		uint32_t n_cigar = ez.n_cigar;
		int adj_size = cigar_adjust(&n_cigar, ez.cigar, false, 15);
		ez.n_cigar = n_cigar;
		return adj_size;
	}

	void printf_alignment_detail(FILE * output, int suggest_st_pos, uint16_t * contig_depth, int ref_pos_enough_match_base){
		log_output = output;
		printCIGAR();
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0) return;
		if(true) printContigSeq(suggest_st_pos);//print contig sequence
		int contig_coverage = 0;
		if(true){ contig_coverage = print_X_E_sequence(suggest_st_pos, ref_pos_enough_match_base); } //print X/= sequence
		if(true){ print_coverage(suggest_st_pos, contig_coverage, contig_depth); } 		//print coverage sequence
		print_SV_canditate(suggest_st_pos);
		fprintf(output, "\n");
	}

	int select_suggest_sv_length(int SV_length, std::vector<int> &suggest_SV_length){
		int best_suggest = 0; int min_dis = 9999999;
		for(int sug: suggest_SV_length){
			float diff_ratio = (float)SV_length/sug - 1;
			diff_ratio = ABS(diff_ratio);
			int dis = ABS_U(SV_length, sug);

			if(diff_ratio > 0.12 && dis > 8)
				continue;

			if(dis < min_dis){
				min_dis = dis;
				best_suggest = sug;
			}
		}
		return best_suggest;
	}

	void get_canditate_SVs(std::vector<NOVA_SV_FINAL_RST_item> &SVs, int minSVlen, int ref_st_pos, std::vector<int> &suggest_SV_length, int region_ref_global_position){
		uint32_t* bam_cigar = ez.cigar;
		int n_cigar = ez.n_cigar;
		int ref_ed_pos = last_ref_end_pos;
		int ref_index = ref_st_pos;
		int read_index = 0;
		uint8_t * ref = NULL;
		uint8_t * alt = NULL;
		int ABS_ref_distance_to_region_end;
		//for SUM deletions
		SUM_variations del_sum; del_sum.clear();
		SUM_variations ins_sum; ins_sum.clear();
		//for SUM insertions
		int suggest_sv_len = 0;
		for(int cigar_ID = 0;cigar_ID < n_cigar; cigar_ID++){
			int c_length =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:	ref_index += c_length; read_index += c_length; break;//M
			case 1:	//insertion
				if(c_length < minSVlen) break;
				ABS_ref_distance_to_region_end = ABS_U(ref_index, ref_ed_pos);
				if(ABS_ref_distance_to_region_end < 20) break;//
				ref = tseq + ref_index; ref_buff.resize(2);//ref
				ref_buff[0] = "ACGT"[ref[0]]; ref_buff[1] = 0;
				alt = qseq + read_index; alt_buff.resize(c_length + 1);//read
				for(int i = 0; i < c_length; i++) alt_buff[i] = "ACGT"[alt[i]];
				alt_buff[c_length] = 0;
				//sum insertions
				ins_sum.store_signal(ref_index + global_ref_pos, read_index, c_length, "ACGT"[qseq[read_index]], cigar_ID);

				//suggest SV length:
				suggest_sv_len = select_suggest_sv_length(c_length, suggest_SV_length);
				if(suggest_sv_len != 0)
					NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, ref_index + global_ref_pos, "INS", &(ref_buff[0]), &(alt_buff[0]),
							suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, cigar_ID, cigar_ID, ref_st_pos, region_ref_global_position);
				read_index += c_length;
				break;//I, int chr_ID; int
			case 2:	 //Deletions
				if(c_length < minSVlen) break;
				ref = tseq + ref_index; ref_buff.resize(c_length + 1);//ref
				for(int i = 0; i < c_length; i++) ref_buff[i] = "ACGT"[ref[i]];
				ref_buff[c_length] = 0;
				alt = qseq + read_index; alt_buff.resize(2);//read
				alt_buff[0] = "ACGT"[alt[0]]; alt_buff[1] = 0;
				//sum deletion:
				del_sum.store_signal(ref_index + global_ref_pos, read_index, c_length, "ACGT"[qseq[read_index]], cigar_ID);

				suggest_sv_len = select_suggest_sv_length(-c_length, suggest_SV_length);
				if(suggest_sv_len != 0)
					NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, ref_index + global_ref_pos, "DEL", &(ref_buff[0]), &(alt_buff[0]),
							suggest_sv_len,  qseq, qlen, bam_cigar, n_cigar, cigar_ID, cigar_ID, ref_st_pos, region_ref_global_position);
				ref_index += c_length;
				break;
			case 3:	ref_index += c_length; read_index += c_length; break;//N, print N
			case 4:	ref_index += c_length; break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}

		//sum deletions:
		SUM_variations &c_sum = del_sum;
		if(c_sum.minSV_filter(minSVlen)){
			ref = tseq + c_sum.ref_idx - global_ref_pos; ref_buff.resize(c_sum.total_var_size + 1);//ref
			for(int i = 0; i < c_sum.total_var_size; i++) ref_buff[i] = "ACGT"[ref[i]];
			ref_buff[c_sum.total_var_size] = 0;
			alt_buff.resize(2);	alt_buff[0] = c_sum.ref_alt_char; alt_buff[1] = 0; //read

			int suggest_sv_len = select_suggest_sv_length(-c_sum.total_var_size, suggest_SV_length);

			if(suggest_sv_len != 0)
				NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, c_sum.ref_idx, "DEL", &(ref_buff[0]), &(alt_buff[0]),
						suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, c_sum.cigar_index_bg, c_sum.cigar_index_ed, ref_st_pos, region_ref_global_position);
		}

		//sum insertion:
		c_sum = ins_sum;
		if(c_sum.minSV_filter(minSVlen)){
			ref_buff.resize(2); ref_buff[0] = c_sum.ref_alt_char; ref_buff[1] = 0;
			alt = qseq + c_sum.read_idx; alt_buff.resize(c_sum.total_var_size + 1);//ref
			for(int i = 0; i < c_sum.total_var_size; i++) alt_buff[i] = "ACGT"[alt[i]];
			alt_buff[c_sum.total_var_size] = 0;

			int suggest_sv_len = select_suggest_sv_length(c_sum.total_var_size, suggest_SV_length);
			if(suggest_sv_len != 0)
				NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, c_sum.ref_idx, "INS", &(ref_buff[0]), &(alt_buff[0]),
						suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, c_sum.cigar_index_bg, c_sum.cigar_index_ed, ref_st_pos, region_ref_global_position);
		}
	}

	inline uint8_t getTseq(int i){ return tseq[i];}

	void setZdrop(uint16_t zdrop_D_, int bandwith_){
		zdrop_D = zdrop_D_;
		bandwith = bandwith_;
	}

private:
	//part1 : basic informations
	int chr_ID;
	int global_ref_pos;
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
	uint16_t zdrop_D; int bandwith; //for DNA zdrop = 400, 200 for RNA
	int flag;
	//part6: logs
	FILE * log_output;
	//part7: candidate SVs
	int last_ref_st_pos;  //the reference start and end position of last alignment
	int last_ref_end_pos;
	std::vector<char> ref_buff;
	std::vector<char> alt_buff;

	void copy_option(){
			match_D = 2;
			mismatch_D= 10;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 0;
			zdrop_D= gap_open2_D + 600; //for DNA zdrop = 400, 200 for RNA
			bandwith = 600;//zdrop_D;
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

	void printCIGAR(){
		fprintf(log_output, "\n");
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0){
			for(int i = 0; i < (int)qlen; i++) {fprintf(log_output, "%c", "ACGT"[qseq[i]]);}
			fprintf(log_output, "\n ???NO alignment result\n");
		}
	}

	void printContigSeq(int suggest_st_pos){
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		int output_index = 0;
		int seq_i = 0;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "%c", "ACGT"[qseq[seq_i]]); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\n");
	}

	int print_X_E_sequence(int suggest_st_pos, int ref_pos_enough_match_base){
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
						if(output_index < ref_pos_enough_match_base) fprintf(log_output, "M");//1 == M; 2 == X; 0 == -;
						else fprintf(log_output, "=");//1 == M; 2 == X; 0 == -;
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

	void print_SV_canditate(int suggest_st_pos){
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		bool have_long_SV = false;
		int minSVlen = 50;

		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	break;//M
			case 1:	case 2: case 3:	case 4:
				if(cigar_len >= minSVlen) have_long_SV = true;
				break;//I, print nothing
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}

		int ref_index = suggest_st_pos;
		if(have_long_SV){
			fprintf(log_output, "Long SV detected:");
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
				case 0:	ref_index += cigar_len; break;//M
				case 1:	if(cigar_len >= minSVlen) {fprintf(log_output, "Insertion: "); fprintf(log_output, "Start %d:%d; length %d: ", chr_ID, ref_index + global_ref_pos, cigar_len);  } seq_i += cigar_len; break;//I, int chr_ID; int
				case 2:	if(cigar_len >= minSVlen) {fprintf(log_output, "Deletion: ");  fprintf(log_output, "Start %d:%d; length %d: ", chr_ID, ref_index + global_ref_pos, cigar_len);  } ref_index += cigar_len; break;//D, print '-'
				case 3:	ref_index += cigar_len; break;//N, print N
				case 4:	ref_index += cigar_len; break;//S, print -
				default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
				}
			}
		}
		fprintf(log_output, "\n");
	}

	void print_coverage(int suggest_st_pos, int contig_coverage, uint16_t * contig_depth){
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "%c", contig_depth[seq_i] + '#'); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "0"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\t");
		fprintf(log_output, "\nCigar sequence: ");
		for(int i = 0; i < cigar_len; i++){
			fprintf(log_output, "%d%c", (bam_cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(bam_cigar[i] & BAM_CIGAR_MASK)]);
		}
		fprintf(log_output, "contig_coverage: [%d bp]\t", contig_coverage);

		fprintf(log_output, "\t");
		fprintf(log_output, "\n");
	}
};

struct read_depth_counter{

private:
	struct read_depth_counter_item{
		int read_idx_bg;
		int read_idx_ed;
		int passed_read_num;
	};

	std::vector<read_depth_counter_item> depth_counter;
	int dc_region_st;
	int dc_shift_offset;//5
	int dc_total_filter_block_num;
	int max_high_read_counter_in_block;

public:
	void init(int region_length, int average_depth, int normal_read_length){
		dc_shift_offset = 5;
		dc_total_filter_block_num = (region_length >> dc_shift_offset);
		max_high_read_counter_in_block = average_depth* (0x1 << dc_shift_offset) / normal_read_length;
		depth_counter.resize(dc_total_filter_block_num);
	}
	read_depth_counter_item * get_rsf_by_pos(int pos){
		int st_idx = pos - dc_region_st;
		//left and right has additional blocks to store additional read
		int sf_idx = (st_idx >> dc_shift_offset);
		if(sf_idx < 0 || sf_idx >= dc_total_filter_block_num){
			return NULL;
		}
		return &(depth_counter[sf_idx]);
	}
	void add_read_depth_item(int read_idx, int pos){
		read_depth_counter_item * rdc = get_rsf_by_pos(pos);
		if(rdc == NULL)	 return;
		rdc->read_idx_bg = MIN(read_idx, rdc->read_idx_bg);
		rdc->read_idx_ed = MAX(read_idx, rdc->read_idx_ed);
	}
	void clear_read_depth_list(){
		for(read_depth_counter_item & rdc : depth_counter){
			rdc.read_idx_ed = 0;
			rdc.read_idx_bg = MAX_int32t;
			rdc.passed_read_num = 0;
		}
	}
	void set_dc_st(int dc_region_st_){
		dc_region_st = dc_region_st_;
	}

	float get_ave_depth(int pos_st, int pos_ed, float *max_depth){
		read_depth_counter_item * rdc_st = get_rsf_by_pos(pos_st);
		read_depth_counter_item * rdc_ed = get_rsf_by_pos(pos_ed);
		if(rdc_st == NULL && rdc_ed == NULL) return 0;
		if(rdc_st == NULL) rdc_st = rdc_ed;
		if(rdc_ed == NULL) rdc_ed = rdc_st;
		int total_depth_read_num = 0;
		int max_depth_read_number = 0;
		for(read_depth_counter_item * c_rdc = rdc_st; c_rdc <= rdc_ed; c_rdc ++){
			if(c_rdc->read_idx_ed - c_rdc->read_idx_bg < 0)
				continue;
			int current_block_read_number = (c_rdc->read_idx_ed - c_rdc->read_idx_bg + 1);
			total_depth_read_num += current_block_read_number;
			max_depth_read_number = MAX(max_depth_read_number, current_block_read_number);
		}
		*max_depth = (float)max_depth_read_number/max_high_read_counter_in_block;
		return (float)total_depth_read_num/(rdc_ed - rdc_st + 1)/max_high_read_counter_in_block;
	}

	bool is_high_coverage(int pos){
		read_depth_counter_item * rdc = get_rsf_by_pos(pos);
		if(rdc == NULL)	 return false;
		return (rdc->read_idx_ed - rdc->read_idx_bg + 1) > (max_high_read_counter_in_block * 2);
	}
};
//Suggest
struct SUGGEST_POS_LIST_ITEM{
	SUGGEST_POS_LIST_ITEM(int suggest_pos_){
		suggest_pos = suggest_pos_;
		read_count = 0;
		low_wrong_base_read_number = 0;//
		high_wrong_base_read_number = 0;
		read_start_pos_sum = 0;
		ave_read_start = 0;
	}

	void add_read_start_pos(int read_start_pos, int wrong_base){
		if(wrong_base < 2)
			low_wrong_base_read_number ++;
		else if(wrong_base < 200)
			high_wrong_base_read_number ++;
		else return;

		read_count++;
		read_start_pos_sum += read_start_pos;
	}

	void count_ave_read_start(){ ave_read_start = (float)read_start_pos_sum/((float)(read_count) + 0.1); }

	void printf(FILE * log){
		if(ave_read_start < -5000 || ave_read_start > 5000)
			fprintf(log, "Fatal ERROR, ave_read_start wrong");
		fprintf(log, "[SG: %d, read count: %d, ave_read_start: %f]\t", suggest_pos, read_count, ave_read_start);
	}

	static inline int cmp_by_read_count(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.read_count > b.read_count; }
	static inline int cmp_by_read_position(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.ave_read_start < b.ave_read_start; }

	int suggest_pos;
	int read_count;
	int low_wrong_base_read_number;
	int high_wrong_base_read_number;
	int read_start_pos_sum;
	float ave_read_start;

};

struct SveHandler{
private:
	//part 1: reference, it is only a reference for the refHandler
	RefHandler *ref;
	//part 2: bam files, the SVE handler keep and maintained the BAM files
	char * bamFileName;
	BAM_handler read;
	int read_counter;
	//buffs for unmapped signals handler
	uint8_t *UMQueryBuff;//buff used in handleUMSignal
	//buff for depth counter
	read_depth_counter rdc;
	std::vector<int> mismatch_position;
	//part 3: BUFFs for local alignment:
	//part 4: buff for signals combination
	float *possibility_r;
	std::vector<SVE> cmb_store_tmp;
	std::vector<int> cmb_try_list;
	//part 5: parameters
	int MIN_sv_len;
	SIGNAL_PARAMETER sig_para;
	//part 6: data for break point distributions, those data used for calculating break point distribution for each signals
	std::vector<float> DR_bp_distribution;
	float MIN_ACCEPT_POSSIBILITY_DR = 0;
	std::vector<float> SH_bp_distribution;
	float MIN_ACCEPT_POSSIBILITY_SH = 0;
	std::vector<float> UM_stPos_distribution;
	float MIN_ACCEPT_POSSIBILITY_UM = 0;
	int START_OFFSET_UM = 0;
	//part 7.0: buffs for SVE signals
	SVE_L sve[SIG::LEN][SV::LEN]; //SVE signals list in the first combining loop
	//SVE_L final; //SVE signals list in the second combining loop

	//part 7: assembly manager
	Ass_Block ass_block_r1;
	Ass_Block ass_block_r2;
	std::vector<uint8_t> bin_contig;
	MainAssemblyHandler *am;
	Contig_String_aligner ca;
	//for contig suggestion alignment positions
	std::map<int, int > suggest_st_pos_map;
	std::set<int> remove_read_set;
	std::vector<SUGGEST_POS_LIST_ITEM> SUGGEST_pos_list;
	std::vector<int> suggest_SV_length;
	int region_ref_global_position;
	//contig realignment
	int region_addition_load;

	std::vector<uint16_t> contig_depth;
	//AssemblyManager *am;
	//part 10: buffs used to VCF records
	FILE* vcf_w;
	kstring_t vcfBuffStr;
	std::vector<NOVA_SV_FINAL_RST_item> SVE_SVs;
	std::vector<NOVA_SV_FINAL_RST_item> region_SVs;

	//genotyping buffs
	Genotyping_read_aligner ga;

public:
	/*************************************PUBLIC FUNCTIONs**********************************************/
	void set_min_accpet_possibility(){
		float MAX_POSSIBILITY_DR = 0;
		for(float poss:DR_bp_distribution)
			MAX_POSSIBILITY_DR = MAX(MAX_POSSIBILITY_DR, poss);
		MIN_ACCEPT_POSSIBILITY_DR = 2*MAX_POSSIBILITY_DR;

		float MAX_POSSIBILITY_SH = 0;
		for(float poss:SH_bp_distribution)
			MAX_POSSIBILITY_SH = MAX(MAX_POSSIBILITY_SH, poss);
		MIN_ACCEPT_POSSIBILITY_SH = 2*MAX_POSSIBILITY_SH;

		float MAX_POSSIBILITY_UM = 0;
		for(float poss:UM_stPos_distribution)
			MAX_POSSIBILITY_UM = MAX(MAX_POSSIBILITY_UM, poss);
		MIN_ACCEPT_POSSIBILITY_UM = 3*MAX_POSSIBILITY_UM;
	}

	void init(RefHandler *ref_, char * BamFileName_, BAM_STATUS_NOVA_SV *bs, char *outputFile, bool is_compression, int MIN_sv_len_, bool output_vcf_header){
		ref = ref_;
		bamFileName = BamFileName_;
		//S1: basic parameters
		sig_para.init(bs);
		//S2: open BAM file
		read.init(bamFileName, ref->get_refFileName());
		//anti registration for BAM header of REF handler
		ref->registHeader(read.file._hdr);
		rdc.init(SEGMENT_LEN, bs->ave_read_depth, bs->analysis_read_length);
		//S3: get distribution
		bs->getBreakPoint_Distribution(DR_bp_distribution, SH_bp_distribution, UM_stPos_distribution, START_OFFSET_UM);
		set_min_accpet_possibility();
		//S4: allocated buffs for signaling
		UMQueryBuff = (uint8_t *)xcalloc(5000, sizeof(uint8_t));
		possibility_r = (float *)xcalloc(5000, sizeof(float));
		//S5: init assembler
		MIN_sv_len = MIN_sv_len_;
		am = new MainAssemblyHandler[1];
		ca.init();
		ga.init();
		//S6 open result vcf files
		vcf_w = xopen(outputFile, "w");
		vcfBuffStr.s = (char *)xcalloc(1000000, sizeof(char));
		vcfBuffStr.l = 0; vcfBuffStr.m = 1000000;
		if(output_vcf_header)
			NOVA_SV_FINAL_RST_item::write_to_vcf_header(vcf_w, read.file._hdr);
	}

	void distory(){
		read.destroy();
		free(UMQueryBuff);
		free(possibility_r);
		delete []am;
		ca.destory();
		ga.destory();
		fclose(vcf_w);
		free(vcfBuffStr.s);
	}

	void print_signal_list(SIG::T sigT, SV::T svt){
		SVE_L * cur_var_signals = &(sve[sigT][svt]);
		fprintf(stderr, "Print %s::%s signals, size %ld\n", SIG::STR[sigT].c_str(), SV::STR[svt].c_str(), cur_var_signals->size());
		for(unsigned int i = 0; i < cur_var_signals->size(); i++)
			std::cerr << cur_var_signals[0][i];
	}

	void SVE_handle_region(){
		clear_sve_signals_in_a_region();//S1
		get_original_signals_from_reads();//S2
		cluster_and_combine_original_signals();	//step3: JC for single sv, multi-sample and multi-signal type
		int signal_list_size;
		SVE_L * cur_var_signals = NULL;

		region_SVs.clear();

		if(true){
			//main type to call INS and DEL
			print_signal_list(SIG::DR, SV::DEL);
			print_signal_list(SIG::SH, SV::INS);
			//other type to call INS
			print_signal_list(SIG::DR, SV::INS);
			//types to call INV and BND
			print_signal_list(SIG::DR, SV::INV_1);
			print_signal_list(SIG::DR, SV::INV_2);
			print_signal_list(SIG::DR, SV::TRA);
			print_signal_list(SIG::DR, SV::TRA_INV);
		}

		//part 0: BND calling process
		am->setRepeatMode();
		fprintf(stderr, "BND(INV) calling process\n" );
		//INV
		bool main_read_before_BP; bool main_read_forward;
		bool supp_read_before_BP; bool supp_read_forward;

		for(int mode = 0; mode < 4; mode++){
			if(mode == 0)			{ cur_var_signals = &(sve[SIG::DR][SV::INV_1]); 	main_read_before_BP = true;   main_read_forward = true; supp_read_before_BP = true; supp_read_forward = true;}//
			else if(mode == 1)		{ cur_var_signals = &(sve[SIG::DR][SV::INV_2]); 	main_read_before_BP = false;  main_read_forward = false; supp_read_before_BP = false; supp_read_forward = false;}
			else if(mode == 2)		{ cur_var_signals = &(sve[SIG::DR][SV::TRA]); 		main_read_before_BP = true;   main_read_forward = true; supp_read_before_BP = false; supp_read_forward = false;}//
			else if(mode == 3)		{continue; /*Not handle TRA_INV*/cur_var_signals = &(sve[SIG::DR][SV::TRA_INV]); 	main_read_before_BP = false;  main_read_forward = false; supp_read_before_BP = true; supp_read_forward = true; }
			if(!cur_var_signals->empty()){
				for(unsigned int i = 0; i < cur_var_signals->size(); i++)
					std::cerr << cur_var_signals[0][i];
				signal_list_size = cur_var_signals->size();
				for(int i = 0; i < signal_list_size; i++){
					if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
					assembly_variations_BND(cur_var_signals[0][i], main_read_before_BP, main_read_forward, supp_read_before_BP, supp_read_forward, region_SVs, read.file._hdr);//todo::
				}
			}
		}

		am->setNormalMode();
		//part 1 large deletions: (using DR/DEL signals)
		fprintf(stderr, "Deletions calling process\n" );
		cur_var_signals = &(sve[SIG::DR][SV::DEL]);
		if(!cur_var_signals->empty()){
			for(unsigned int i = 0; i < cur_var_signals->size(); i++)
				std::cerr << cur_var_signals[0][i];
			signal_list_size = cur_var_signals->size();
			for(int i = 0; i < signal_list_size; i++){
				if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
				set_region_addition_load_long();
				int total_region_length = cur_var_signals[0][i].r2.ed_pos - cur_var_signals[0][i].r1.st_pos;
				if(total_region_length > 1800){
					if(total_region_length < 5000)
						set_region_addition_load_extramly_long();
					else if(total_region_length < 50000)
						set_region_addition_load_super_super_long();
					else
						continue;
				}
				assembly_variations(cur_var_signals[0][i]);
			}
		}

		//part 2: small variations: (using SH/INS signals) (< 300 bp ins/del/dup)
		am->setRepeatMode();
		fprintf(stderr, "Small variations calling process\n" );
		cur_var_signals = &(sve[SIG::SH][SV::INS]);
		for(unsigned int i = 0; i < cur_var_signals->size(); i++)
			std::cerr << cur_var_signals[0][i];
		//std::sort(cur_var_signals.begin(), cur_var_signals.end(), SVE::cmp_by_position);//needed to be sorted?
		if(!cur_var_signals->empty()){
			signal_list_size = cur_var_signals->size();
			for(int i = 0; i < signal_list_size; i++){
				fprintf(stderr, "\nCurrent SVE ID: [%d]", i);
				if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
				set_region_addition_load_short();
				assembly_variations(cur_var_signals[0][i]);
			}
		}

		//part 3: BND calling process
		fprintf(stderr, "True SVs output\n" );
		std::sort(region_SVs.begin(), region_SVs.end(), NOVA_SV_FINAL_RST_item::cmp_by_position);

		//remove duplication
		int SV_number = region_SVs.size();
		for(int c_SV_idx = 0; c_SV_idx < SV_number - 1; c_SV_idx++)
			for(int cmp_sv_idx = c_SV_idx + 1;
					cmp_sv_idx < SV_number && region_SVs[c_SV_idx].duplication_SV_filter(region_SVs[cmp_sv_idx]);
					cmp_sv_idx++);

		for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: region_SVs ){
			//running genotyping process:
			//if(sv.SV_is_duplicated()) continue;

			fprintf(stderr, "Genotyping for: \n ");
			sv.print(stderr, read.file._hdr, ref);
			sv.genotyping(sig_para.MaxReadLen, sig_para.insert_size_max, &read.file, &ga, ref);

			sv.writeVCF(&vcfBuffStr, read.file._hdr);
			fprintf(vcf_w, "%s", vcfBuffStr.s);
		}
	}

private:
	/*************************************BASIC FUNCTIONs**********************************************/
	void clear_sve_signals_in_a_region(){//used only within function "process"
		for(int i = 0; i < SIG::LEN; i++) for(int j = 0; j < SV::LEN; j++) sve[i][j].clear();
		read.clear();
	}

	SveHandler (SveHandler &B);//SveHandler can`t be copied
	/*************************************SIGNALS FUNCTIONs**********************************************/
	//S1: get signals from all reads, in this step, a read pair may be used as DR or SA or UM signals, for each read signals, a SVE will be stored //signals functions
	bool is_mate_fwd(uint16_t flag) { return (!((flag & BAM_MATE_STRAND) != 0));}
	void get_original_signals_from_reads();
	void handleDRSignal(bam1_core_t *core, int middle_size);
	void storeMismatchSignals(bam1_t *br, READ_record &c_r);
	void storeClipSignals(bool isClipAtRight, uint32_t ori_pos, uint8_t read_mapq);
	void handleSASignal(READ_record& c_r, bam1_t *br);//handle a SA signal and store results
	void handleUMSignal(bam1_t *br);

	/*************************************CLUSTERING FUNCTIONs**********************************************/
	//S2: combine signal from multiple reads. In this steps, signals from multiple reads will be clustered and joint together. At last, signals combined from many reads(SVEs) will form.
	void cluster_and_combine_original_signals();
	int  getTopPossibilityIdx(int r_min, int r_max, SVE_L & l, bool isR1, bool is_forward, float &max_poss, std::vector<float> &dis_bp_percent);//used in sve_combine_STEP1_DR
	void single_type_sve_combine(SVE_L & l, int min_score, SIG::T sigT, SV::T svt);
	void combine_duplication(SVE_L & l);
	void DR_SH_signal_combining(SVE_L &SH, SVE_L &DR);//S3: Joint Calling for single sv, multi-sample and multi-signal type //r1 and r2 of SH must be within r1 and r2 of DR

	/*************************************ASSEMBLY FUNCTIONs**********************************************/
	//S4: read assembly
	bool SV_region_depth_filter(SVE & c_sve); //return whether pass assembly filter
	bool assembly_load_read(SVE &sv, RefRegion &main, RefRegion &supp);
	bool assembly_variations(SVE &sv);
	bool assembly_variations_BND(SVE &sv, bool main_read_before_BP, bool main_read_forward, bool supp_read_before_BP, bool supp_read_forward, std::vector<NOVA_SV_FINAL_RST_item> &region_SVs, bam_hdr_t * header);//buff to store contig

	//S4.1: assembling contigs re-alignment
	void get_suggention_alignment_position_list(AssemblyContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_reads_info> &ass_read_list, FILE* log_f );
	void getSuggestSVlength();
	void alignment_and_get_var(int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len, int ref_region_length);

	void set_region_addition_load_short(){
		region_addition_load = 500;
		ca.setZdrop(1200, 1200);
	}

	void set_region_addition_load_long(){
		region_addition_load = 2000;
		ca.setZdrop(2000, 2000);
	}

	void set_region_addition_load_extramly_long(){
		region_addition_load = 5000;
		ca.setZdrop(5000, 5000);
	}

	void set_region_addition_load_super_super_long(){
		region_addition_load = 50000;
		ca.setZdrop(50000, 50000);
	}

	//abandoned functions
	void jointAssemblingUmReads();	//abandoned
	void umAssembly(int beginIdx);	//abandoned
	int getTopPossibilityIdx_UM(int r_min, int r_max, float &max_poss);	//abandoned
	/*************************************VCF writing FUNCTIONs**********************************************/
};

bam_hdr_t* checkBamHeadersAreSAME(std::vector<SveHandler>& h);
#endif /* SRC_SIGNAL_SVEHANDLER_HPP_ */
