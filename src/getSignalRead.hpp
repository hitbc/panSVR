/*
 * getSignalRead.hpp
 *
 *  Created on: 2021年3月16日
 *      Author: fenghe
 */

#ifndef GETSIGNALREAD_HPP_
#define GETSIGNALREAD_HPP_

#include <iostream>

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#include "cpp_lib/get_option_cpp.hpp"

extern "C"
{
#include "clib/utils.h"
#include "clib/bam_file.h"
#include "clib/desc.h"
//#include "clib/vcf_file.h"
}

#define GAP_OPEN 16
#define GAP_EXT 1
#define GAP_OPEN2 32
#define GAP_EXT2 0
#define MATCH_SCORE 2
#define MISMATCH_SCORE 12

#define MAX_TID 24

struct READ_SIGNAL_HANDLER{
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
	bool NOT_USING_FILTER;
	bool discard_both_full_match;
	double sample_rate;
	int sample_max_number_int;

	Bam_file c_b;
	bam_hdr_t* hdr = NULL;
	uint32_t reasonFlagCounter[1024];
	int penalty1, penalty2, penalty;

	uint32_t getScoreByCigar(bam1_t* b);
	//return 0 when no XA, 1~5 when has <= 5 XA , 6 when has more than 6 XA
	//when mapq == 0, BWA-MEM return at most 5 alt alignment when score are same, other wise XA is blank
	void all_signal_records_read_pair(bam1_t &read1, bam1_t &read2, bool will_Be_used);
	bool signal_be_used();
	void SURVIVOR_SV_region_get_all_signal_records_SORT_BY_pos(char *input_bam_fn);
	void SURVIVOR_SV_region_get_all_signal_records_SORT_BY_NAME();

	int init_run(int argc, char *argv[]){

		//option list
		options_list l;

	    l.add_title_string("\n");
	    l.add_title_string("Usage:     ");  l.add_title_string(PACKAGE_NAME);
	    l.add_title_string("  signal  [Options] <BAM/CRAM file> \n");

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
		l.add_option("not-use-filter",  'D', "do not using signal filter but output all reads as signal"); l.set_arg_pointer_back((void *)&NOT_USING_FILTER);
		l.add_option("discard-full-match",  'U', "-U has higher priority level than -D, discard read pair when both reads are aligned without any clip or error"); l.set_arg_pointer_back((void *)&discard_both_full_match);
		double default_sample_rate = 1;
		l.add_option("sample-rate",  	'R', "Only use part of the reads, when set to be 1, all reads are used, when set to be 0.5, only random selected half of reads", true, default_sample_rate); l.set_arg_pointer_back((void *)&sample_rate);

		fprintf(stderr, "V1.23\n");
		if(l.default_option_handler(argc, argv)) return 1;

		if (argc - optind < 1)	return l.output_usage();
		char *input_bam_fn = argv[optind];

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

		memset(reasonFlagCounter, 0, 1024* sizeof(uint32_t));
		if(bam_sort_by_name)	SURVIVOR_SV_region_get_all_signal_records_SORT_BY_NAME();
		else					SURVIVOR_SV_region_get_all_signal_records_SORT_BY_pos(input_bam_fn);

		bam_file_close(&c_b);
		for(int i = 0; i < 128; i++)
			fprintf(stderr, "flag:\t%d\tcount:\t%d\n", i, reasonFlagCounter[i]);
		return 0;
	}

};

#endif /* GETSIGNALREAD_HPP_ */
