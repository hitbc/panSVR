#include<stdio.h>
#include<string.h>
#include <time.h>

#include"PanSVgenerateVCF/getSignalRead.hpp"
#include"cpp_lib/get_option_cpp.hpp"
#include"cpp_lib/Assembler/mantaAssembler.hpp"
#include "PanSVgenerateVCF/get_anchor_ref.hpp"
#include "PanSVgenerateVCF/read_realignment.hpp"
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
	deCOY_CLASSIFY_MAIN cm;
	cm.init_run(argc, argv);
	return 0;
}

int get_signal_main(int argc, char *argv[]){
	READ_SIGNAL_HANDLER *sig = (READ_SIGNAL_HANDLER *) xcalloc(1, sizeof(READ_SIGNAL_HANDLER));
	sig->init_run(argc, argv);
	return 1;
}

int SV_ref_signal_main(int argc, char *argv[]){
	 VCF_HANDLER v;
	 v.run(argc, argv);
	return 1;
}

int remapped_read_assembly(int argc, char *argv[]);
int SVLoci_main(int argc, char *argv[]);

int manta_test_main(int argc, char *argv[]){
	assembly_test();
	return 1;
}

int main(int argc, char *argv[])
{
	fprintf(stderr, "\n\n main version V2.00\n\n");
	COMMAND_HANDLER ch;
	ch.add_function("sv_calling", "used for MANTA like SV calling", SVLoci_main);

	ch.add_help_msg_back("");
	ch.add_help_msg_back("");
	ch.add_help_msg_back("The following functions used for pan-genome based SV force calling");
	ch.add_help_msg_back("");

	ch.add_function("fc_anchor_ref", "The 1st step of force calling (S1)", SV_ref_signal_main);
	ch.add_help_msg_back("Getting anchor reference from original reference and VCF file");

	ch.add_function("fc_index", "The 2ed step of force calling (S2)", build_deBGA_index);
	ch.add_help_msg_back("Building index for all anchor reference. Using [fc_anchor_ref] command before building index");

	ch.add_function("fc_signal", "The 3rd step of force calling (S3)", get_signal_main);
	ch.add_help_msg_back("Getting signal reads from bam files");

	ch.add_function("fc_aln", "The 4th step of force calling (S4)", classify_main);
	ch.add_help_msg_back("Aligning all signal reads into anchor reference");

	ch.add_function("fc_sv", "The 5th step of force calling (S5)", remapped_read_assembly);
	ch.add_help_msg_back("Assembling all re-mapped read in each anchor reference region and generated SV results");

	ch.add_help_msg_back("");
	ch.add_help_msg_back("");

	ch.add_function("tools", "Small analysis tools", analysis_main);

	ch.add_function("assembly_test", "a test for assembly program ", manta_test_main);

	return ch.run(argc, argv);
}
