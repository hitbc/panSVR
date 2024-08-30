#ifndef GET_ANCHOR_REF_HPP_
#define GET_ANCHOR_REF_HPP_

extern "C"
{
#include "../clib/utils.h"
#include "../clib/bam_file.h"
#include "../clib/vcf_file.h"
}

#include "../cpp_lib/get_option_cpp.hpp"

//------------------------------MAIN--------------------------------------//
#define DUP_MAX_LEN 2000

//'ALL' 'DEL' 'INS' 'DUP' 'TRA' 'INV' 'BND'");
struct SV_TYPE{
	static const int ALL = 0;
	static const int DEL = 1;
	static const int INS = 2;
	static const int DUP = 3;
	static const int TRA = 4;
	static const int INV = 5;
	static const int BND = 6;
	static const int UNKNOWN = 7;

	static int get_sv_type_int(char * sv_type_str){
			 if(strcmp(sv_type_str, "DEL") == 0) return SV_TYPE::DEL;
		else if(strcmp(sv_type_str, "INS") == 0) return SV_TYPE::INS;
		else if(strcmp(sv_type_str, "DUP") == 0) return SV_TYPE::DUP;
		else if(strcmp(sv_type_str, "TRA") == 0) return SV_TYPE::TRA;
		else if(strcmp(sv_type_str, "INV") == 0) return SV_TYPE::INV;
		else if(strcmp(sv_type_str, "BND") == 0) return SV_TYPE::BND;
		else if(strcmp(sv_type_str, "ALL") == 0) return SV_TYPE::ALL;
		else if(strcmp(sv_type_str, "all") == 0) return SV_TYPE::ALL;
		return SV_TYPE::UNKNOWN;
	}
};

//used to store buff into FA string, and adding a \n for each 70 DNA base
struct FA_string_BUFF{

	void buff_init(char * str_buff_){
		str_buff_ori = str_buff = str_buff_;
		total_string_len = 0;
		finish = false;
	}

	//stop storing when reach the str_ed
	void storeData(char * str_bg, char * str_ed){
		xassert(finish == false,"");
		for(; str_bg < str_ed ; str_bg++, total_string_len++){
			sprintf(str_buff, "%c", *str_bg); str_buff += strlen(str_buff);
			if(total_string_len %70 == 69)	{sprintf(str_buff, "\n");  str_buff++; }
		}
	}

	//overload stop storing when reach 0
	void storeData(char * str_end_with_0){
		xassert(finish == false,"");
		for(; *str_end_with_0 != 0 ; str_end_with_0++, total_string_len++){
			sprintf(str_buff, "%c", *str_end_with_0); str_buff += strlen(str_buff);
			if(total_string_len %70 == 69)	{sprintf(str_buff, "\n");  str_buff++; }
		}
	}

	char * get_string(){
		if(finish == false){
			if((total_string_len - 1) % 70 != 69){
				sprintf(str_buff, "\n"); str_buff++;
			}
		}
		finish = true;
		return str_buff_ori;
	}

	int get_string_len(){
		// for the final string
		if(finish == false){
			if((total_string_len - 1) % 70 != 69){
				sprintf(str_buff, "\n"); str_buff++;
			}
		}
		finish = true;
		return total_string_len;
	}

private:
	bool finish;
	int total_string_len;
	char * str_buff; //1M
	char * str_buff_ori; //1M
};

#define condition_chech(CONDITION, REASON_STRING)  if(CONDITION) { fprintf(stderr, "Skip %d because %s \n", load_ref_ID - 1, REASON_STRING); continue; }

struct VCF_HANDLER{
	char *fasta_fn = NULL;
	char *vcf_fn = NULL;
	int edge_len = 0;
	int minSV_len = 0;
	bool begin_at_0;
	//filers
	int all_sample = false;
	int all_type = false;
	int all_chrom = false;
	char *special_sample_name;
	char *special_sv_type;
	char *special_CHROM_ID;
	int special_sv_type_int;
	int special_CHROM_ID_int;

	bool discard_decoy_seq;
	bool skip_N_ref;
	bool skip_blank_ref_or_allele;

	//files
	BCF_FILE vcf_r;//vcf for read
	bcf_hdr_t *header;//vcf header
	faidx_t * fai;//index file for fasta

	int run(int argc, char *argv[]){
		fprintf(stderr, "\n\n V1.01\n\n");
		options_list l;
		add_options(l);
		l.show_command(stderr, argc + 1, argv - 1);
		//get optional parameters
		if(l.default_option_handler(argc, argv)){
			l.show_c_value(stderr); return 1;
		}
		l.show_c_value(stderr);
		if (argc - optind < 2) return l.output_usage();
		//get necessary parameters
		fasta_fn = strdup(argv[optind]);
		vcf_fn = strdup(argv[optind + 1]);

		filter_init();//filter:
		fai = fai_load(fasta_fn);//load reference
		xassert(VCF_open(&vcf_r, vcf_fn) == true, "");	header = vcf_r.header;//load vcf files

		//main function
		getSV_ref();

		//finish program
		bcf_close(vcf_r.file);
		fai_destroy(fai);
		return 0;
	}

	void add_options(options_list &l){
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  fc_anchor_ref  [Options] [ref.fa] [input.vcf]>\n");
	    l.add_title_string("  Basic:   \n");
	    l.add_title_string("    [ref.fa(.gz)]  FILE    reference files, must be the reference used to generate VCF file\n");
	    l.add_title_string("    [input.vcf]    FILE    the input vcf file to generate the reference\n");
	    l.add_title_string("                           the result output into stdout\n");
		l.add_option("edge-len", 		'e', "Additional reference around the break point", true, 500); l.set_arg_pointer_back((void *)&edge_len);
		l.add_option("minSV-len", 		'm', "min SV length, SV shorter than it will be ignored", true, 50); l.set_arg_pointer_back((void *)&minSV_len);
		l.add_option("begin_at_0", 		'b', "the position is begin at 0 in a vcf"); l.set_arg_pointer_back((void *)&begin_at_0);
			l.add_help_msg_back("set this option for pbsv, and ignore it for cute SV");
		//filters
		l.add_option("sample-name", 	'S', "Assigning a special sample for output, others will be filtered out", true, "ALL"); l.set_arg_pointer_back((void *)&special_sample_name);
			l.add_help_msg_back("To use this function, the VCF must have 'SAMPLE=XXX' tag");
		l.add_option("sv-type", 		'T', "Assigning a special SV type for output, others will be filtered out", true, "ALL"); l.set_arg_pointer_back((void *)&special_sv_type);
			l.add_help_msg_back("Normally should be one of 'ALL' 'DEL' 'INS' 'DUP' 'TRA' 'INV' 'BND'");
		l.add_option("CHROM_ID", 		'I', "Assigning a special chr_ID ['ALL', 0~9999] for output, others will be filtered out", true, "ALL"); l.set_arg_pointer_back((void *)&special_CHROM_ID);
		l.add_option("discard_decoy", 	'J', "Discard SV in decoy region, it implemented by setting a MIN_ref_length (40M), any reference shorter than it will be discard."); l.set_arg_pointer_back((void *)&discard_decoy_seq);
		l.add_help_msg_back("Besides, chr_name not start by 'c', '1' ~'9', 'X', 'Y' will be discard, too");
		l.add_help_msg_back("'discard_decoy' used to cutoff SV in 'decoy' sequence in hs37d5 or grch38, the shortest chromosome chr21 is 48M, while decoy is 35M ");
		l.add_option("skip-N-ref", 		'N', "Skip the SV when the reference is start by N"); l.set_arg_pointer_back((void *)&skip_N_ref);
		l.add_option("skip-<>-ref", 	'B', "Skip the SV the reference or the allele which is begin with '<', like <DEL> or <INS>", true, true); l.set_arg_pointer_back((void *)&skip_blank_ref_or_allele);
		l.add_help_msg_back("When not skip, the string like <INS> will be just leave Blank");
	}

	void filter_init(){
		//filter:
		if(strcmp("all", special_sample_name) == 0 || strcmp("ALL", special_sample_name) == 0)	all_sample = true;
		special_sv_type_int = SV_TYPE::get_sv_type_int(special_sv_type); if(special_sv_type_int == SV_TYPE::ALL)	all_type = true;
		if(strcmp("all", special_CHROM_ID) == 0 || strcmp("ALL", special_CHROM_ID) == 0)		all_chrom = true;
		else special_CHROM_ID_int = strtol(special_CHROM_ID,0,10);
	}

	bool try_filter(bcf1_t *c_r, int c_sv_type_int){
		char c_sample[100];
		if(all_sample == false){ vcf_get_sample(vcf_r.header, c_r, c_sample);if(*c_sample != 0 && strcmp(c_sample, special_sample_name) != 0) return false;	}
		if(all_type == false && special_sv_type_int != c_sv_type_int)	return false;
		if(all_chrom == false && special_CHROM_ID_int != c_r->rid)		return false;

		//SV type check
		if(c_sv_type_int == SV_TYPE::BND) return false;
		if(c_sv_type_int == SV_TYPE::TRA) return false;
		if(c_sv_type_int == SV_TYPE::INV) return false;
		if(c_sv_type_int == SV_TYPE::DUP)
			if(c_r->rlen > DUP_MAX_LEN)   return false;
		return true;
	}

	//this function check whether the reference string in "REF" field of VCF line(REF_field_STR) is same as the the reference string
	//in the reference FASTA file, when the string are not same, a WARNING will be output into STDERR
	void CHECK_REF_filed_VS_true_ref(char * REF_field_STR , int sv_type_int, int sv_ref_length_by_POS_END, int sv_ID, char * true_ref_string){
		//reference length check
		int32_t REF_field_STR_Len = strlen(REF_field_STR); //only used for reference check
		if(sv_type_int == SV_TYPE::DUP)	{
			if(REF_field_STR_Len != 1)
				fprintf(stderr, "Warning： Ref length not equal to VCF data, vcf show ref [%s] has length %d @ [VCF ID: %d], ",
						REF_field_STR, sv_ref_length_by_POS_END, sv_ID);
		}
		else{
			if(sv_ref_length_by_POS_END != REF_field_STR_Len){
				fprintf(stderr, "Warning： reference string length in 'REF' field of VCF line not equal to SV REF length,"
						"SV of REF length %d has 'REF' filed of string [%s] @ [VCF ID: %d], ", sv_ref_length_by_POS_END, REF_field_STR, sv_ID);
				if(REF_field_STR_Len > sv_ref_length_by_POS_END){ REF_field_STR[sv_ref_length_by_POS_END] = 0; }
				fprintf(stderr, "so the reference will be set to be %s with length %d \n", REF_field_STR, sv_ref_length_by_POS_END);
				REF_field_STR_Len = sv_ref_length_by_POS_END;
			}
		}

		if(REF_field_STR[0] == 'N' || REF_field_STR[0] == 'n'){
			fprintf(stderr, "Warning , REF string is 'N' or 'n' @ [VCF ID: %d] \n", sv_ID);
		}
		else{
			//reference string check, for string not 'N'
			//printf("ori_ref: %s_%s\n", reg, ref_seq); //debug code:
			if(sv_type_int == SV_TYPE::DUP) true_ref_string ++;
			if(strncmp(REF_field_STR, true_ref_string, REF_field_STR_Len) != 0){
				char old_char = true_ref_string[REF_field_STR_Len];
				true_ref_string[REF_field_STR_Len] = 0;
				if(REF_field_STR_Len < 1000){
					fprintf(stderr, "Warning: string in 'REF' filed is not same as the string in true reference file,\n"
						"string in 'REF' filed is %s; string in true reference is :%s  @ [VCF ID: %d] \n ", REF_field_STR, true_ref_string, sv_ID);
				}else{
					fprintf(stderr, "Warning: string in 'REF' filed is not same as the string in true reference file,\n"
											"string in 'REF' filed is %s; string in true reference is :<too long to show>  @ [VCF ID: %d] \n ", REF_field_STR, sv_ID);
				}
				true_ref_string[REF_field_STR_Len] = old_char;
			}
		}
	}

	void getSV_ref(){
		int new_ref_ID = 0;
		int load_ref_ID = 0;
		char * str_buff_ori = (char *)xmalloc(1000000); //1M
		FA_string_BUFF fa_str;

		while(VCF_next(&vcf_r)){//read one //load a new vcf data
			load_ref_ID ++;
			bcf1_t *c_r = &( vcf_r.r);
			//unpack the vcf data to get the alt string
			bcf_unpack(c_r, BCF_UN_STR);

			//SV format check
			char * SV_ref_in_vcf = c_r->d.allele[0];
			condition_chech((skip_blank_ref_or_allele && SV_ref_in_vcf[0] == '<'), "reference begin with : ['<']");//like <INS>
			condition_chech((skip_N_ref && SV_ref_in_vcf[0] == 'N') , "option: [skip_N_ref]");//skip BND

			//SV type
			char sv_type_buff[100];
			vcf_get_sv_type(header, c_r, sv_type_buff);
			int sv_type_int = SV_TYPE::get_sv_type_int(sv_type_buff);
			condition_chech((try_filter(c_r, sv_type_int) == false) , "[type filter]");//type filter

			const char * ref_chr_name = bcf_hdr_id2name(header, c_r->rid);
			if(discard_decoy_seq){
				condition_chech((faidx_seq_len(fai, ref_chr_name) < 40000000) , " option: [discard_decoy_seq]"); //40M
				char bc = ref_chr_name[0];
				condition_chech((bc != 'c' && bc != 'C' && bc != 'X' && bc != 'Y' && (bc > '9' || bc < '0')) , " option: [discard_decoy_seq]"); //40M
			}

			//get load reference string region:
			if(begin_at_0) c_r->pos++;
			int ref_pos = c_r->pos;
			int load_ref_st_pos = ref_pos - edge_len;
			int ori_ref_pos_in_load_string = edge_len;
			if(load_ref_st_pos <= 0){
				ori_ref_pos_in_load_string = ref_pos - 1;
				load_ref_st_pos = 0;
			}
			int load_ref_ed_pos = ref_pos + c_r->rlen + edge_len;

			//int faidx_seq_len(const faidx_t *fai, const char *seq);
			//const char * ref_chr_name = faidx_iseq(fai, c_r->rid);
			//get ref region string
			char reg[1024];
			sprintf(reg, "%s:%d-%d", ref_chr_name, load_ref_st_pos, load_ref_ed_pos);

			//load reference string;
			int true_region_load_len = 0;
			char *ref_seq = fai_fetch(fai, reg, &true_region_load_len);

			CHECK_REF_filed_VS_true_ref( c_r->d.allele[0], sv_type_int, c_r->rlen, load_ref_ID - 1, ref_seq + ori_ref_pos_in_load_string);

			//get all new reference string for
			for(uint32_t i = 1; i < c_r->n_allele; i++){
				int32_t allele_len = strlen(c_r->d.allele[i]);
				condition_chech((skip_blank_ref_or_allele && c_r->d.allele[i][0] == '<'), "allele begin with : ['<']");//like <INS>
				//SV length check
				condition_chech((c_r->rlen < minSV_len && allele_len < minSV_len), "option: [minSV_len]");
				fa_str.buff_init(str_buff_ori);

				if(sv_type_int == SV_TYPE::DUP){
					//get duplication length:
					int dupLength = 0;
					vcf_get_sv_LENGTH(header, c_r, &dupLength);

					fa_str.storeData(ref_seq, ref_seq + ori_ref_pos_in_load_string + c_r->rlen);
					fa_str.storeData(ref_seq + ori_ref_pos_in_load_string);
				}
				else{
					ref_seq[ori_ref_pos_in_load_string] = 0; //set END OF STRING for the first part of the reference
					//print first part of the reference string
					fa_str.storeData(ref_seq);
					fa_str.storeData(c_r->d.allele[i]);
					fa_str.storeData(ref_seq + ori_ref_pos_in_load_string + c_r->rlen);
				}
				//get new string from REF and ALT

				int st_in_ref = ref_pos - edge_len;
				int bp1_in_ref = ref_pos;
				int bp2_in_ref = ref_pos + c_r->rlen;
				int end_in_ref = ref_pos + c_r->rlen + edge_len;
				printf(">%d_%s_%d_%d_%s_%d_%d_%d_%s\n",new_ref_ID++, ref_chr_name, st_in_ref, fa_str.get_string_len(), sv_type_buff,
						bp1_in_ref, bp2_in_ref, end_in_ref, c_r->d.id ? c_r->d.id : ".");
				printf("%s",fa_str.get_string());
			}
			if(ref_seq != NULL){ free(ref_seq); ref_seq = NULL;	} 		//free reference sequence
		}
	}

};

#endif /* GET_ANCHOR_REF_HPP_ */
