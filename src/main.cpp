/*************************************************************************
	> File Name: main.c
	> Author: 
	> Mail: 
	> Created Time: 2018年05月08日 星期二 18时44分23秒
 ************************************************************************/

#include<stdio.h>
#include<string.h>
#include <time.h>

#include"jlra_aln.hpp"
#include"getSignalRead.hpp"
#include"getSV_ref.hpp"
#include"cpp_lib/MSA/abPOA_handler.hpp"
#include"cpp_lib/get_option_cpp.hpp"
#include"cpp_lib/Assembler/mantaAssembler.hpp"
extern "C"
{
#include "clib/utils.h"
#include "clib/desc.h"
}

int analysis_main(int argc, char *argv[]);

int classify_main(int argc, char *argv[])
{
	fprintf(stderr, "\n\n V1.21\n\n");
	//get option
	CLASSIFY_MAIN cm;
	cm.init_run(argc, argv);
	return 0;
}

int get_signal_main(int argc, char *argv[]){
	READ_SIGNAL_HANDLER sig;
	sig.init_run(argc, argv);
	return 1;
}

int SV_ref_signal_main(int argc, char *argv[]){
	 VCF_HANDLER v;
	 v.run(argc, argv);
	return 1;
}

int remapped_read_assembly(int argc, char *argv[]);

int abPOA_main(int argc, char *argv[]){
	abPOA_test();

	return 1;
}

int manta_test_main(int argc, char *argv[]){
	assembly_test();
	return 1;
}

int main(int argc, char *argv[])
{
	COMMAND_HANDLER ch;
	ch.add_function("sv_ref", "get reference from original reference and VCF file", SV_ref_signal_main);
	ch.add_function("signal", "get all signal reads from a bam file", get_signal_main);

	ch.add_function("index", "building index for all SV region", build_deBGA_index);
	ch.add_help_msg_back("Get SV region reference using [sv_ref] command before building reference index");

	ch.add_function("aln", "align all signal reads into reference", classify_main);
	ch.add_help_msg_back("Get signal reads using [signal] command before alignment ");

	ch.add_function("assembly", "analysis the coverage region of vcf and assemble all remapped read in that region", remapped_read_assembly);
	ch.add_help_msg_back("Align signal reads using [aln] command before assembly ");

	ch.add_function("analysis", "small analysis functions", analysis_main);
	ch.add_function("abPOA", "a test for the abPOA", abPOA_main);
	ch.add_function("assembly_test", "a test for manta assembly program ", manta_test_main);

	return ch.run(argc, argv);
}
