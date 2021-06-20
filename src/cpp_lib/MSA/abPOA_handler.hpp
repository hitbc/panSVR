/*
 * abPOA_handler.hpp
 *
 *  Created on: 2021年3月18日
 *      Author: fenghe
 */

#ifndef ABPOA_HANDLER_HPP_
#define ABPOA_HANDLER_HPP_

#include <vector>

extern "C"
{
//#include "clib/utils.h"
#include "abPOA/abpoa.h"
}

struct MSA_READ{
	int bg_in_cons;
	int ed_in_cons;
	std::vector<uint8_t> seq;

	MSA_READ(int bg_in_cons_, const char * seq_){	set_new(bg_in_cons_, seq_);	}
	MSA_READ(int bg_in_cons_, uint8_t * seq_, int str_len){ set_new(bg_in_cons_, seq_, str_len); }

	void set_new(int bg_in_cons_, const char * seq_);
	void set_new(int bg_in_cons_, uint8_t * seq_, int str_len);
};

struct abPOA_handler{
    //options
    bool out_cons;
    bool out_msa;
    bool out_png;
    abpoa_t *ab;
    abpoa_para_t *abpt;

    //inputs
    std::vector<MSA_READ> read_list;
    int cons_len;

    //buffs
    std::vector<int> out_degree_buff;
    //output
    std::vector<char> consensus_str;


    void init(){
	    // initialize variables
	    ab = abpoa_init();
	    abpt = abpoa_init_para();

	    // output options
	    out_cons = true;
	    out_msa = true;
	    out_png = false;
	    if(out_cons)	abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
	    if(out_msa)     abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable

	    abpoa_post_set_para(abpt);

	    consensus_str.resize(10000);
	    out_degree_buff.resize(10000);
    }

	void destory() {
	    abpoa_free(ab); abpoa_free_para(abpt);
	}

	void clear(){
		cons_len = 0;
		read_list.clear();
	}

	//run MSA
    // perform abpoa-msa
    void run_msa(){
	    ab->abs->n_seq = read_list.size();
	    abpoa_res_t res;
	    for (int i = 0; i < ab->abs->n_seq; ++i) {
	    	MSA_READ &r = read_list[i];
	        res.graph_cigar = 0, res.n_cigar = 0; res.is_rc = 0;
	        abpoa_align_sequence_to_subgraph(ab, abpt, r.bg_in_cons, r.ed_in_cons, &(r.seq[0]), r.seq.size(), &res);
	        int exc_beg, exc_end;
	        if (i != 0) abpoa_subgraph_nodes(ab, r.bg_in_cons, r.ed_in_cons, &exc_beg, &exc_end);
	        else exc_beg = 0, exc_end = 0;
	        abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, &(r.seq[0]), r.seq.size(), res, i, ab->abs->n_seq);
	        if (res.n_cigar) free(res.graph_cigar);
	    }

	    if(out_cons){
	        abpoa_generate_consensus(ab, abpt, stdout, NULL, NULL, NULL, NULL);
	    }
	    	//cons_len = abpoa_generate_consensus_simple(ab, abpt, &consensus_str[0], &out_degree_buff[0]);

	    if(out_msa)
	    	abpoa_generate_rc_msa(ab, abpt, stdout, NULL, NULL);

	    if(out_png){
		    abpt->out_pog = strdup("/home/fenghe/sub_example.png"); // dump parital order graph to file
		    if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);
	    }
    }
};

void abPOA_test();

//demo codes
//void abPOA_test(){
//
//	abPOA_handler ah;
//	ah.init();
//
//	char *seqs_demo_p[6];
//	int * bg_p[6];
//	for(int i = 0; i < 6; i++){
//		seqs_demo_p[i] = seqs_demo[i];
//		bg_p[i] = beg_end_id_demo[i];
//	}
//	ah.abpt->out_msa = 1;//output MSA
//	ah.add_reads_M2(seqs_demo_p, 6, bg_p);
//	ah.run_msa();
//	ah.getMSA();
//
//}


#endif /* ABPOA_HANDLER_HPP_ */
