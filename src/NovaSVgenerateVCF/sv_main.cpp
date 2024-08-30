/*
 * sv_main.cpp
 *
 *  Created on: 2020年4月27日
 *      Author: fenghe
 */

#include <getopt.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include "../NovaSVgenerateVCF/SveHandler.hpp"
#include "../NovaSVgenerateVCF/NovaSVOption.hpp"

int SVLoci_main(int argc, char *argv[])
{
	NOVA_SV_PARA * nova_para = (NOVA_SV_PARA *)xcalloc(1, sizeof(NOVA_SV_PARA));
	if(nova_para->get_option(argc, argv) != 0)
		return 0;

	RefHandler *refHandler = (RefHandler *)xcalloc(1,sizeof(RefHandler));
	refHandler->init(nova_para);

	BAM_STATUS_NOVA_SV *bs = (BAM_STATUS_NOVA_SV *)xcalloc(1, sizeof(BAM_STATUS_NOVA_SV));
	bs->generate_state(nova_para->referenceFilename, nova_para->bamFile); //todo::

	//bs->load(nova_para->statsFileName);
	SveHandler *sve_h = (SveHandler *)xcalloc(1, sizeof(SveHandler));
	sve_h->init(refHandler, nova_para->bamFile, bs, nova_para->outputFile, nova_para->is_compression, nova_para->MIN_SV_len, !nova_para->not_output_vcf_header);
//	//PART5: for each choromosome
	while(refHandler->load_seg_index())//get a new ref block in the reference, 2M bases per block
		sve_h->SVE_handle_region();//S1: in step1: we clear all signals in the SVE list
	sve_h->distory();
	return 0;
}
