#ifndef GETSIGNALREAD_HPP_
#define GETSIGNALREAD_HPP_

#include <iostream>

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "../cpp_lib/get_option_cpp.hpp"
#include "../cpp_lib/statistics/StatsManager.hpp"
extern "C"
{
#include "../clib/utils.h"
#include "../clib/bam_file.h"
#include "../clib/desc.h"
//#include "clib/vcf_file.h"
}

#define GAP_OPEN 16
#define GAP_EXT 1
#define GAP_OPEN2 32
#define GAP_EXT2 0
#define MATCH_SCORE 2
#define MISMATCH_SCORE 12

#define MAX_TID 24

#define MAX_ISIZE 100000 //0.1M
#define MAX_ANA_READ_LEN 1000 //1000
#define HUMAN_GENOME_SIZE 3100000000 //3.1G

struct BAM_STAT{

	//global analysis variables
	uint32_t reasonFlagCounter[1024];
	uint64_t *isize_analysis;
	uint64_t *read_length_analysis;
	double read_normal_percent;
	uint64_t total_read_number;
	//analysis results before getting all signals by sampling
	uint32_t minInsertLen;
	uint32_t middleInsertLen;
	uint32_t maxInsertLen;
	std::vector<float> isize_distribution;
	int analysis_read_length;
	double ave_read_len;
	//analysis results after all data is processed
	uint32_t minInsertLen_l2;
	uint32_t maxInsertLen_l2;
	double ave_read_depth_sample;//get ave read depth by sampling
	double ave_read_depth;

	bool already_init = false;

	void init(){
		//memory init
		isize_analysis = (uint64_t *)xmalloc(MAX_ISIZE * sizeof(uint64_t));
		read_length_analysis = (uint64_t *)xmalloc(MAX_ANA_READ_LEN * sizeof(uint64_t));
		reset();
	}

	void reset(){
		total_read_number = 0;
		memset(isize_analysis, 0, MAX_ISIZE* sizeof(uint64_t));
		memset(read_length_analysis, 0, MAX_ANA_READ_LEN* sizeof(uint64_t));
		memset(reasonFlagCounter, 0, 1024* sizeof(uint32_t));
	}

	void destory(){
		free(isize_analysis);
		free(read_length_analysis);
	}

	void collect_signal(bam1_t &c_read){
		//global analysis:
		//get ISIZE:
		int isize = ABS(c_read.core.isize);
		if(isize > 0 && isize < MAX_ISIZE){
			isize_analysis[isize]++;
		}
		//get read_length
		if(c_read.core.l_qseq < MAX_ANA_READ_LEN)
			read_length_analysis[c_read.core.l_qseq]++;
	}

	void global_analysis_stat(){
		//output analysis results
		//part1: get normal read length
		analysis_read_length = -1;
		double total_read_len = 0;
		for(int i = 0; i < MAX_ANA_READ_LEN; i++){
			total_read_len += i*read_length_analysis[i];
			if(read_length_analysis[i] > 0.6*total_read_number){
				analysis_read_length = i; read_normal_percent = (double)read_length_analysis[i]/total_read_number; break;
			}
		}
		ave_read_len = total_read_len/total_read_number;
		if(analysis_read_length == -1)
			analysis_read_length = ave_read_len;
		//part2:  get read depth:
		ave_read_depth = (double)analysis_read_length * total_read_number /HUMAN_GENOME_SIZE;
		//part3: max and min ISIZE
		minInsertLen_l2 = 0; maxInsertLen_l2 = 0;
		float max_percent_abnormal_isize = 0.01;
		int max_abnormal_read_number_one_side = max_percent_abnormal_isize * total_read_number;
		int read_number_sum = 0;
		for(int i = 0; i < MAX_ISIZE; i++){
			read_number_sum += isize_analysis[i];
			if(read_number_sum > max_abnormal_read_number_one_side){
				minInsertLen_l2 = i; break;
			}
		}
		read_number_sum = 0;
		for(int i = MAX_ISIZE - 1; i > 0; i--){
			read_number_sum += isize_analysis[i];
			if(read_number_sum > max_abnormal_read_number_one_side){
				maxInsertLen_l2 = i; break;
			}
		}
	}

	void sampling_analysis_stat(const char * ref_fn, char * input_bam_fn, bool bam_sort_by_name){
		//reset
		if(already_init){ reset(); } else{ init(); already_init = true; }
		//get read length and simple ISIZE
		//it is recommended to use sorted bam/cram files
		Bam_file c_b;
		memset(&c_b, 0, sizeof(Bam_file));
		bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
		bam_hdr_t* hdr = c_b._hdr;
		total_read_number = 0;
		bam1_t b1 = {0};//BAM record for the first read in a pair
		int sam_rst1 = 0;
		while (1){
			//load SAM 1 & 2
			do{	sam_rst1 = sam_read1(c_b._hfp, hdr, &b1); } while( (sam_rst1 >= 0) && (bam_is_secondary(&b1) || bam_is_supplementary(&b1)));
			if(sam_rst1 < 0)		break;
			total_read_number ++;
			if(total_read_number == 100000)
				break;
			collect_signal(b1);
		}
		bam_file_close(&c_b);
		global_analysis_stat();

		if(!bam_sort_by_name){
			//get ISZIE by sampling (from MANTA)
			StatsManager rstats(ref_fn, "");
			rstats.handleBamCramStats(input_bam_fn, &ave_read_depth_sample);
			minInsertLen = rstats.getInsertLen(input_bam_fn, 0.01f);
			middleInsertLen = rstats.getInsertLen(input_bam_fn, 0.5f);
			maxInsertLen = rstats.getInsertLen(input_bam_fn, 0.99f);

			unsigned idx = rstats.getGroupIndex(StatLabel(input_bam_fn, ""));
			int totalPairedReadCount = rstats.getStats(idx).readCounter.totalHighConfidenceReadPairCount() + 1;
			const SizeDistribution &fragStats = rstats.getStats(idx).fragStats;
			for(uint32_t i = minInsertLen; i < maxInsertLen; i++){
				int index_count = fragStats.getSizeCount(i);
				isize_distribution.emplace_back((float)index_count/totalPairedReadCount);
			}
		}else{
			//it is recommended to use sorted bam/cram files
			minInsertLen = minInsertLen_l2;
			middleInsertLen = (minInsertLen_l2 + maxInsertLen_l2)/2;
			maxInsertLen = maxInsertLen_l2;
			for(uint32_t i = minInsertLen; i < maxInsertLen; i++){
				int index_count = isize_analysis[i];
				isize_distribution.emplace_back((float)index_count/(total_read_number + 1));
			}
		}
		ave_read_depth = ave_read_depth_sample;
		fprintf(stderr, "BAM/CRAM status: read length: [Normal: %d @ %f%%, AVE: %f] ISIZE: [MIN: %d MIDDLE:%d MAX: %d] ave_read_depth [%f]\n", analysis_read_length, read_normal_percent*100, ave_read_len, minInsertLen, middleInsertLen, maxInsertLen, ave_read_depth_sample);
	}

	void output_full_stat(FILE * output){
		for(int i = 0; i < 128; i++)
			fprintf(output, "flag:\t%d\tcount:\t%d\n", i, reasonFlagCounter[i]);
	}

	void output_final_stat(FILE * output){
		fprintf(output, "%f_%d_%d_%d_%d_%d\n", ave_read_depth, analysis_read_length, minInsertLen_l2, maxInsertLen_l2, minInsertLen, maxInsertLen);
		for(uint32_t i = minInsertLen; i < maxInsertLen; i++){
			fprintf(output, "%f\n", isize_distribution[i - minInsertLen]);
		}
	}
};

struct READ_SIGNAL_HANDLER{
	//output options
	FILE * output_file1 = NULL;
	FILE * output_file2 = NULL;

	//score options
	int gap_open;
	int gap_ex;
	int gap_open2;
	int gap_ex2;
	int match;
	int mismatch;

	//other options
	int maxTid;
	bool bam_sort_by_name;
	bool not_filter_low_quality;
	char *ref_fn;
	const char *out_header_fn;
	const char *out_status_fn;
	const char *tmp_file_pairing;
	bool NOT_USING_FILTER;
	bool discard_both_full_match;
	double sample_rate;
	int sample_max_number_int;

	//buffs
	Bam_file c_b;
	bam_hdr_t* hdr = NULL;

	//global analysis variables
	BAM_STAT bs;
	int isize_max;
	int isize_min;

	bool already_output_isize = false;

	//functions
	uint32_t getScoreByCigar(bam1_t* b);
	//return 0 when no XA, 1~5 when has <= 5 XA , 6 when has more than 6 XA
	//when mapq == 0, BWA-MEM return at most 5 alt alignment when score are same, other wise XA is blank
	void all_signal_records_read_pair(bam1_t &read1, bam1_t &read2, bool will_Be_used);
	void stat_analysis(bam1_t &c_read);
	bool signal_be_used();
	void SURVIVOR_SV_region_get_all_signal_records_SORT_BY_pos(char *input_bam_fn);
	void SURVIVOR_SV_region_get_all_signal_records_SORT_BY_NAME();

	int init_run(int argc, char *argv[]){
		//option list
		options_list l;
		l.show_command(stderr, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("Usage:     ");  l.add_title_string(PACKAGE_NAME);
	    l.add_title_string("  fc_signal  [Options] <BAM/CRAM file> \n");

	    l.add_title_string("Basic: [BAM/CRAM file]  FILES    input sam/bam/cram file, only one file can be accept\n");
	    l.add_title_string("   For a cram file, a reference file is needed;\n");
	    l.add_title_string("   Input sam/bam/cram should be sorted by alignment position, when input sorted by name, ");
	    l.add_title_string("using [--sort-by-name] option, \n");
	    l.add_title_string("   The output file will be output into stdout as fastq format, ");
	    l.add_title_string("using [ | pigz -p 8 -- > [out_fn.fq.gz]] to compact it\n");

		//score
		l.add_option("gap-open1", 		'O', "Gap open penalty 1", true, GAP_OPEN); l.set_arg_pointer_back((void *)&gap_open);
		l.add_help_msg_back("a k-long gap costs min{O+k*E,P+k*F}..");
		l.add_option("gap-open2", 		'P', "Gap open penalty 2.", true, GAP_OPEN2); l.set_arg_pointer_back((void *)&gap_open2);
		l.add_option("gap-extension1", 	'E', "Gap extension penalty 1.", true, GAP_EXT); l.set_arg_pointer_back((void *)&gap_ex);
		l.add_option("gap-extension2", 	'F', "Gap extension penalty 2.", true, GAP_EXT2); l.set_arg_pointer_back((void *)&gap_ex2);
		l.add_option("match-score", 	'M', "Match score for SW-alignment.", true, MATCH_SCORE); l.set_arg_pointer_back((void *)&match);
		l.add_option("mis-score", 		'm', "Mismatch score for SW-alignment.", true, MISMATCH_SCORE); l.set_arg_pointer_back((void *)&mismatch);

		//other part
		l.add_option("max-tid-filter",  'I', "filter, the read will be treated as signal when tid > [max-tid-filter]", true, MAX_TID); l.set_arg_pointer_back((void *)&maxTid);
		l.add_option("sort-by-name",  	'N', "the input file sorted by name"); l.set_arg_pointer_back((void *)&bam_sort_by_name);
		l.add_option("not-ignore-low-q", 'L', "do not ignore the NM or clip filter in low quality read segment"); l.set_arg_pointer_back((void *)&not_filter_low_quality);
		l.add_option("reference",  		'r', "the reference file used for CRAM file", false, ""); l.set_arg_pointer_back((void *)&ref_fn);
		l.add_option("header-file",     'H', "output BAM/CRAM header file of input file", true, "./header.sam"); l.set_arg_pointer_back((void *)&out_header_fn);
		l.add_help_msg_back("this file will be used in [aln] command");
		l.add_option("status-file",     'S', "output BAM/CRAM status file of input file", true, "./status.sam"); l.set_arg_pointer_back((void *)&out_status_fn);
		l.add_help_msg_back("Include informations about depth, AVERAGE ISZIE and normal read length, etc. This may be used in assembly command");
		l.add_option("tmp_file_pairing", 't', "A tmp file used to store unpaired SAM records", true, "./tmp_file_pairing"); l.set_arg_pointer_back((void *)&tmp_file_pairing);


		l.add_option("not-use-filter",  'D', "do not using signal filter but output all reads as signal"); l.set_arg_pointer_back((void *)&NOT_USING_FILTER);
		l.add_option("discard-full-match",  'U', "-U has higher priority level than -D, discard read pair when both reads are aligned without any clip or error"); l.set_arg_pointer_back((void *)&discard_both_full_match);
		double default_sample_rate = 1;
		l.add_option("sample-rate",  	'R', "Only use part of the reads, when set to be 1, all reads are used, when set to be 0.5, only random selected half of reads", true, default_sample_rate); l.set_arg_pointer_back((void *)&sample_rate);
		fprintf(stderr, "V1.23\n");
		if(l.default_option_handler(argc, argv)){
			l.show_c_value(stderr);
			return 1;
		}
		l.show_c_value(stderr);

		if (argc - optind < 1)	return l.output_usage();
		char *input_bam_fn = argv[optind];

		bs.sampling_analysis_stat(ref_fn, input_bam_fn, bam_sort_by_name);
		bs.output_final_stat(stderr);
		//output analysis result
		FILE * status_output = xopen(out_status_fn, "w");
		bs.ave_read_depth *= sample_rate;
		bs.output_final_stat(status_output);
		fclose(status_output);


		bs.reset();
		isize_max = bs.maxInsertLen + 150;
		isize_min = bs.minInsertLen - 150; if(isize_min < 1) isize_min = 1;
		//debug::test
//		FILE * status_output = xopen(out_status_fn, "w");
//		bs.ave_read_depth *= sample_rate;
//		bs.output_final_stat(status_output);
//		fclose(status_output);
//		if(1){ exit(0); }

		output_file1 = stdout;//xzopen(out_fn_read1, "wb");
		output_file2 = stdout;//xzopen(out_fn_read2, "wb");
		if(sample_rate < 0.9999){
			sample_max_number_int = sample_rate * RAND_MAX;
			fprintf(stderr, "Sample_rate: [%f] sample_max_number_int: [%d]\n", sample_rate, sample_max_number_int);
		}
		//read bam/cram file:
		memset(&c_b, 0, sizeof(Bam_file));
		bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
		hdr = c_b._hdr;
		//print header
		if(out_header_fn != NULL){
			htsFile *out_header_file = hts_open(out_header_fn, "w");//open output file
			xassert(sam_hdr_write(out_header_file, hdr) >=0, "write header wrong!");//write header
			hts_close(out_header_file);
		}

		if(bam_sort_by_name)	SURVIVOR_SV_region_get_all_signal_records_SORT_BY_NAME();
		else					SURVIVOR_SV_region_get_all_signal_records_SORT_BY_pos(input_bam_fn);

		bam_file_close(&c_b);

		bs.global_analysis_stat();
		fprintf(stderr, "BAM/CRAM status: ave_read_depth: [%f] read length: [Normal: %d @ %f%%, AVE: %f] ISIZE: [MIN: %d MAX: %d]\n", bs.ave_read_depth,
				bs.analysis_read_length, bs.read_normal_percent*100, bs.ave_read_len, bs.minInsertLen_l2, bs.maxInsertLen_l2);
		return 0;
	}

};

#endif /* GETSIGNALREAD_HPP_ */
