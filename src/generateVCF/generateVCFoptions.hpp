/*
 * GenerateVCFoptions.hpp
 *
 *  Created on: 2021年5月31日
 *      Author: fenghe
 */

#ifndef GENERATEVCF_GENERATEVCFOPTIONS_HPP_
#define GENERATEVCF_GENERATEVCFOPTIONS_HPP_

extern "C"
{
#include "../clib/utils.h"
#include "../clib/desc.h"
}

struct SignalAssemblyOption{
	//option for read loader
	int st_chr_ID; int st_pos;
	int ed_chr_ID; int ed_pos;
	int max_read;
	int max_load_read_per_time;
	int MIN_score;

	//option for SV assembler
	int edge_len;
	char *input_bam_fn;

	//for reference
	char *index_dir;
	char * ori_header_fn;
	char * ori_reference_fn;

	//output
	char *output_vcf;

	//bool print detail information
	bool print_detail_information;
	bool show_coverage_depth_analysis;

	int max_cluster_distance;

	int get_option(int argc, char *argv[]){
		fprintf(stderr, "\n\n V1.015\n\n");

		options_list l;

		l.add_title_string("\n");
		l.add_title_string("Usage:     ");  l.add_title_string(PACKAGE_NAME);
		l.add_title_string("  assembly  [Options] <IndexDir> [BAM file] [ori_header_fn.sam] [reference.fa]\n");

	    l.add_title_string("Basic: n");
	    l.add_title_string("<IndexDir>      FOLDER   the directory contains deBGA index\n");
		l.add_title_string("[BAM file]  FILES    input sam/bam file, only one file can be accept, only sam/bam can be accept\n");
		l.add_title_string("   Input sam/bam should be sorted by alignment position");
		l.add_title_string("   The coverage reads will be output into stdout\n");
		l.add_title_string("   The analysis result will be output into stderr\n");
		l.add_title_string("[ori_header.sam]  FILES  Header file of original BAM/CRAM file, using [signal] command or [samtools view -H] to generate it\n");
		l.add_title_string("[reference.fa]  FILES  Original reference file to generate SV database, like hs37d5 or GRCH38.\n");

		l.add_option("st_chr_ID", 		'S', "The start chr_ID for analysis", true, 0); l.set_arg_pointer_back((void *)&st_chr_ID);
		l.add_option("st_pos", 			's', "The start position for analysis", true, 0); l.set_arg_pointer_back((void *)&st_pos);
		l.add_option("ed_chr_ID", 		'E', "The end chr_ID for analysis", true, 10000); l.set_arg_pointer_back((void *)&ed_chr_ID);
		l.add_option("ed_pos", 			'F', "The end position for analysis(500M)", true, 500000000); l.set_arg_pointer_back((void *)&ed_pos);
		l.add_option("max_read", 		'N', "MAX number of read loaded in whole analysis", true, MAX_int32t); l.set_arg_pointer_back((void *)&max_read);
		l.add_option("load-size", 		'L', "MAX number of read loaded into memory per time(10M)", true, 10000000); l.set_arg_pointer_back((void *)&max_load_read_per_time);
		l.add_option("MIN_score", 		'M', "MIN alignment score of reads to be able to used for analysis", true, 50); l.set_arg_pointer_back((void *)&MIN_score);
		l.add_option("edge-len", 		'e', "Additional reference around the break point in reference", true, 200); l.set_arg_pointer_back((void *)&edge_len);
		l.add_help_msg_back("Should be same as 'edge-len' used in 'sv_ref'");
		l.add_option("print-detail", 	'D', "Print original read information and assembly contigs in stderr"); l.set_arg_pointer_back((void *)&print_detail_information);
		l.add_option("depth-detail", 	'd', "Print depth coverage analysis information in stderr"); l.set_arg_pointer_back((void *)&show_coverage_depth_analysis);
		l.add_option("cluster-dis", 	'c', "Max distance of SVs to be treat as one SV cluster", true, 150); l.set_arg_pointer_back((void *)&max_cluster_distance);
		l.add_help_msg_back("SVs that has same SV type and has distance less than [cluster-dis] will be clustered into one block");
		l.add_option("output", 			'o', "VCF output file, if not set, the result will be output into stdout", true, "stdout"); l.set_arg_pointer_back((void *)&output_vcf);

		if(l.default_option_handler(argc, argv)) return 1;

		if (argc - optind < 3)	return l.output_usage();
		index_dir = argv[optind];
		input_bam_fn = argv[optind + 1];
		ori_header_fn = argv[optind + 2];
		ori_reference_fn = argv[optind + 3];

		return 0;

	}
};



#endif /* GENERATEVCF_GENERATEVCFOPTIONS_HPP_ */
