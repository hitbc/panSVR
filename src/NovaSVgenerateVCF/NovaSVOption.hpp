/*
 * NovaSVOption.hpp
 *
 *  Created on: 2021年8月18日
 *      Author: fenghe
 */

#ifndef NOVASVGENERATEVCF_NOVASVOPTION_HPP_
#define NOVASVGENERATEVCF_NOVASVOPTION_HPP_

#include "../cpp_lib/get_option_cpp.hpp"

struct NOVA_SV_PARA{
	char*	referenceFilename;
	char* 	outputFile;
	char* 	statsFileName;
	char*   bamFile;
	bool 	is_compression;
	int 	MIN_SV_len;

	bool not_output_vcf_header;

	int st_chr_ID; int st_pos;
	int ed_chr_ID; int ed_pos;


	int get_option(int argc, char *argv[]){
		char default_outputFile[] = "./novaSV.vcf";
		options_list l;
		l.show_command(stderr, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  sv_calling  [Options] [novaSignal.bam]\n");
	    l.add_title_string("  Basic:   \n");
	    l.add_title_string("    <novaSignal.bam>   FILES  sorted BAM files\n");

		//thread number
		//l.add_option("thread", 			't', "Number of threads", true, 4); l.set_arg_pointer_back((void *)&thread_n);
		l.add_option("MIN_SV_len", 		'm', "MIN length of output SVs", true, 30); l.set_arg_pointer_back((void *)&MIN_SV_len);
		l.add_option("reference",  		'r', "the reference (.fa)", false, ""); l.set_arg_pointer_back((void *)&referenceFilename);
		l.add_option("output",  		'o', "the output file (.vcf)", true, default_outputFile); l.set_arg_pointer_back((void *)&outputFile);
		l.add_option("bcf",				'b', "when set to be true, output BCF, not VCF");l.set_arg_pointer_back((void *)&is_compression);
		l.add_option("status",  		'T', "Status file name", false, ""); l.set_arg_pointer_back((void *)&statsFileName);
		l.add_help_msg_back("");
		l.add_help_msg_back("The following options used for multiple threads SV calling");
		l.add_help_msg_back("Running SV calling REGION by REGION then 'cat' the results");
		l.add_option("no-header", 		'H', "Not output VCF header"); l.set_arg_pointer_back((void *)&not_output_vcf_header);
		l.add_option("st_chr_ID", 		'S', "The start chr_ID for SV calling", true, 0); l.set_arg_pointer_back((void *)&st_chr_ID);
		l.add_option("st_pos", 			's', "The start position for SV calling", true, 0); l.set_arg_pointer_back((void *)&st_pos);
		l.add_option("ed_chr_ID", 		'E', "The end chr_ID for SV calling", true, 10000); l.set_arg_pointer_back((void *)&ed_chr_ID);
		l.add_option("ed_pos", 			'F', "The end position for SV calling(500M)", true, 500000000); l.set_arg_pointer_back((void *)&ed_pos);

		if(l.default_option_handler(argc, argv)) {
			l.show_c_value(stderr);
			return 1;
		}
		l.show_c_value(stderr);
		if (argc - optind < 1)
			return l.output_usage();
		bamFile = argv[optind];
		return 0;
	}
};


#endif /* NOVASVGENERATEVCF_NOVASVOPTION_HPP_ */
