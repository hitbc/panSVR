/*
 * SV_ref_sequence.hpp
 *
 *  Created on: 2021年5月31日
 *      Author: fenghe
 */

#ifndef GENERATEVCF_SV_REF_SEQUENCE_HPP_
#define GENERATEVCF_SV_REF_SEQUENCE_HPP_
#include <stdlib.h>
#include <vector>
#include "../deBGA_index.hpp"
extern "C"
{
#include "../clib/utils.h"
#include "../clib/bam_file.h"
#include "../clib/desc.h"
}

struct SV_cluster_item{

	SV_cluster_item(int sv_ID_, int cluster_number_, int begin_SV_ID_, int next_SV_ID_){
		set(sv_ID_, cluster_number_, begin_SV_ID_);
		next_SV_ID = next_SV_ID_;
	}

	void set(int sv_ID_, int cluster_number_, int begin_SV_ID_){
		sv_ID = sv_ID_;
		cluster_number = cluster_number_;
		begin_SV_ID = begin_SV_ID_;
	}

	void print(FILE * output){ fprintf(output, "ID: %d[ B %d; N %d; C %d] -->> ", sv_ID, begin_SV_ID, next_SV_ID, cluster_number); }

	int sv_ID;//the SV ID, used for check
	int cluster_number;//the number of SV in this cluster

	int begin_SV_ID;//the first SV in the cluster, used for traverse all the SVs in a cluster
	int next_SV_ID;//the SV next to this one in the cluster, if it is the final one, next will be -1. This used for traverse all the SVs in a cluster
};

struct SV_ref_sequence{

public:

	int init(char *index_dir, char * ori_header_fn, char * ori_reference_fn, int max_cluster_distance)
	{
		//load bin reference //buffer_ref_seq:store reference seq;load from dm file ref.seq
		result_ref_seq = (load_data_from_file(index_dir, "ref.seq", ((void ** )&buffer_ref_seq), true, 536) >> 3);

		//********************************************************************************************
		//load original sam header from file(grch37d5 or chr38)
		htsFile *ori_header_file = hts_open(ori_header_fn, "r");//open output file
		ori_header = sam_hdr_read(ori_header_file);
		hts_close(ori_header_file);

		//********************************************************************************************
		//load original reference from file(grch37d5 or chr38)
		original_ref_fai = fai_load(ori_reference_fn);
		original_ref_pointer = NULL;

		//read chr names and length from file [unipath.chr]
		char unichr[ROUTE_LENGTH_MAX] = {0};
		char chr_line_content[MAX_CHR_NAME_LENGTH];
		strcpy(unichr, index_dir);
		if (unichr[strlen(unichr) - 1] != '/') strcat(unichr, "/");
		strcat(unichr, "unipath.chr");
		FILE *fp_chr = xopen (unichr, "r" );
		uint32_t chr_line_n = 0;
		int fscanf_rst = fscanf(fp_chr,"%s",chr_line_content);
		if(fscanf_rst == EOF) fprintf(stderr, "Warning, error fscanf_rst @ %s @ line %d\n ", __func__,__LINE__);
		xassert(1, "");
		while(!feof(fp_chr)){
			if ((chr_line_n & 0X1) == 0)	{
				sv_info.emplace_back(chr_line_content, ori_header);
			}
			else				            {
				uint32_t c_chr_end_n;
				sscanf(chr_line_content, "%u", &c_chr_end_n);
				chr_end_n.emplace_back(c_chr_end_n);
			}
			fflush(stdout);
			chr_line_n++;
			int fscanf_rst = fscanf(fp_chr,"%s",chr_line_content);
			if(fscanf_rst == EOF) fprintf(stderr, "Warning, error fscanf_rst @ %s @ line %d\n ", __func__,__LINE__);
		}
		fclose(fp_chr);

		build_SV_Cluster(max_cluster_distance);
		sv_already_used = (bool *)xcalloc(sv_info.size(), sizeof(bool));
		if(false){
			printCluster();
			xassert(0, "Test END\n");
		}
		return 0;
	}

	void destory(){
		free(sv_c);
		free(buffer_ref_seq);
		fai_destroy(original_ref_fai);
	}

	char *get_chr_name_by_ID(int chr_ID){
		return ori_header->target_name[chr_ID];
	}

	void output_VCF_header(FILE* vcf_header){
		//##fileDate=2021-05-23|12:54:26PM|CST|+0800
		//##source=SVIM-asm-v1.0.2
		//##contig=<ID=1,length=249250621>
		//##contig=<ID=hs37d5,length=35477943>
		fprintf(vcf_header, "##fileformat=VCFv4.2\n");
		fprintf(vcf_header, "##source=%s-%s\n", PACKAGE_NAME, PACKAGE_VERSION);
		for(int chr_ID = 0; chr_ID < ori_header->n_targets; chr_ID++)
			fprintf(vcf_header, "##contig=<ID=%s,length=%d>\n", ori_header->target_name[chr_ID], ori_header->target_len[chr_ID]);
		//##ALT=<ID=DEL,Description="Deletion">
		//##ALT=<ID=INV,Description="Inversion">
		//##ALT=<ID=DUP,Description="Duplication">
		//##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
		//##ALT=<ID=DUP:INT,Description="Interspersed Duplication">
		//##ALT=<ID=INS,Description="Insertion">
		//##ALT=<ID=BND,Description="Breakend">
		fprintf(vcf_header, "##ALT=<ID=DEL,Description=\"Deletion\">\n");
		fprintf(vcf_header, "##ALT=<ID=INV,Description=\"Inversion\">\n");
		fprintf(vcf_header, "##ALT=<ID=DUP,Description=\"Duplication\">\n");
		fprintf(vcf_header, "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n");
		fprintf(vcf_header, "##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">\n");
		fprintf(vcf_header, "##ALT=<ID=INS,Description=\"Insertion\">\n");
		fprintf(vcf_header, "##ALT=<ID=BND,Description=\"Breakend\">\n");

		//##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
		//##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description="Genomic origin of interspersed duplication seems to be deleted">
		//##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
		//##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
		//##FILTER=<ID=not_fully_covered,Description="Tandem duplication is not fully covered by a contig">
		//##FILTER=<ID=incomplete_inversion,Description="Only one inversion breakpoint is supported">
		//##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		//##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number of tandem duplication (e.g. 2 for one additional copy)">
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample

		fprintf(vcf_header, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
		fprintf(vcf_header, "##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">\n");
		fprintf(vcf_header, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
		fprintf(vcf_header, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
		fprintf(vcf_header, "##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a contig\">\n");
		fprintf(vcf_header, "##FILTER=<ID=incomplete_inversion,Description=\"Only one inversion breakpoint is supported\">\n");
		fprintf(vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		fprintf(vcf_header, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of tandem duplication (e.g. 2 for one additional copy)\">\n");
		fprintf(vcf_header, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample\n");
	}

	//this function get the reference sequence of SV.fa by SV ID
	//meanwhile, function "get_original_ref_seq" get the sequence in original reference like hs37d5, by chr_ID, begin and end position
	uint8_t* get_sv_str_by_id(int sv_id, int &sv_length){
		xassert((uint32_t)sv_id < chr_end_n.size(), "");
		sv_length =  sv_info[sv_id].region_len;
		if((int)seq_buff.size() < sv_length)
			seq_buff.resize(sv_length);
		uint32_t start = (sv_id == 0)?0:chr_end_n[sv_id - 1];
		get_refseq(&(seq_buff[0]), sv_length, start);
		return 	&(seq_buff[0]);
	}

	//get original reference from file hs37D5 or GRCH38
	char * get_original_ref_seq(int chr_ID, int st_pos, int ed_pos){
		if(original_ref_pointer != NULL) free(original_ref_pointer);
		original_ref_pointer = NULL;
		char reg[1024];
		sprintf(reg, "%s:%d-%d", ori_header->target_name[chr_ID], st_pos, ed_pos);
		int true_region_load_len = 0;
		original_ref_pointer = fai_fetch(original_ref_fai, reg, &true_region_load_len);
		if(true_region_load_len < ed_pos - st_pos + 1)
			return NULL;
		return original_ref_pointer;
	}

	SV_chr_info* get_sv_info(int SV_ID){ return &(sv_info[SV_ID]); }
	bool get_SV_already_used(int SV_ID){ return sv_already_used[SV_ID]; }
	void set_SV_already_used(int SV_ID){ sv_already_used[SV_ID] = true; }
	SV_cluster_item*  get_cluster(int SV_ID){ return sv_c + SV_ID; }

	std::vector<SV_chr_info> * get_sv_info_list(){return &(sv_info);}

private:

	faidx_t * original_ref_fai;//index file for fasta
	char * original_ref_pointer;
	uint64_t result_ref_seq;
	uint64_t* buffer_ref_seq;//original reference
	std::vector<uint32_t> chr_end_n;
	bam_hdr_t *ori_header;
	std::vector<uint8_t> seq_buff;
	std::vector<SV_chr_info> sv_info;
	bool * sv_already_used;

	void get_refseq(uint8_t *ref, uint32_t len, uint32_t start){
		uint32_t m;
		for (m = 0; m < len; ++m)
			ref[m] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
	}

	SV_cluster_item* sv_c;

	void build_SV_Cluster(int max_cluster_distance){
		int sv_list_size = sv_info.size();
		sv_c = (SV_cluster_item * )xcalloc(sv_list_size, sizeof(SV_cluster_item));

		int cb_i = 0; // current_begin_item
		for(; cb_i < sv_list_size; cb_i++){
			if(sv_c[cb_i].cluster_number != 0) continue; //skip already clustered SVs
			//initiation a new cluster
			int cluster_number = 1;
			int chr_ID = sv_info[cb_i].chr_ID;
			int begin_pos = sv_info[cb_i].st_pos;
			std::string &sv_type = sv_info[cb_i].sv_type;
			//try to find cluster items:
			int ct_i = cb_i + 1;  // current try search item
			int old_index = cb_i;
			for(;ct_i < sv_list_size; ct_i++){
				int try_st_pos = sv_info[ct_i].st_pos;
				int dis_st_pos = ABS(try_st_pos - begin_pos);
				if(dis_st_pos > max_cluster_distance)
					break;
				if(chr_ID == (int)sv_info[ct_i].chr_ID && sv_type == sv_info[ct_i].sv_type){
					cluster_number ++;
					begin_pos = MAX(begin_pos, try_st_pos);
					sv_c[old_index].next_SV_ID = ct_i;
					old_index = ct_i;
				}
			}
			sv_c[old_index].next_SV_ID = -1;
			//set other information
			for(int next_ID = cb_i; next_ID != -1; next_ID = sv_c[next_ID].next_SV_ID)
				sv_c[next_ID].set(next_ID, cluster_number, cb_i);
		}
	}

	//debug function
	void printCluster(){
		int total_number = sv_info.size();
		for(int i = 0; i< total_number; i++){
			if(sv_c[i].begin_SV_ID != sv_c[i].sv_ID)
				continue;
			fprintf(stderr, "\n\n");
			int st_pos = sv_info[i].st_pos;
			for(int next_ID = sv_c[i].begin_SV_ID ; next_ID != -1; next_ID = sv_c[next_ID].next_SV_ID){
				fprintf(stderr, "\n");
				fprintf(stderr, "[%d %s]", (int)sv_info[next_ID].st_pos - st_pos, sv_info[next_ID].vcf_print_string.c_str());
				sv_c[next_ID].print(stderr);
			}
		}
	}

};


// ************************************************************************************************************//

#endif /* GENERATEVCF_SV_REF_SEQUENCE_HPP_ */
