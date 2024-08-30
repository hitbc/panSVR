/*
 * SignalAssembly.cpp
 *
 *  Created on: 2021年4月20日
 *      Author: fenghe
 */

#include "../PanSVgenerateVCF/SignalAssembly.hpp"

#include "../PanSVgenerateVCF/samTag.hpp"
//'ALL' 'DEL' 'INS' 'DUP' 'TRA' 'INV' 'BND'");

//************************************************************************************************************//
//return read number, when no data be registered left, return -1 ,otherwise return load read number
int AssemblyHandler::read_register(Remapped_read_loader &rl, SV_ref_sequence &sf){
	SV_already_used = false;
	if(rl.get_next_sv_region(&start_index, &load_read_number, &sv_id) == false)	return -1; //the final sv region
	if(sf.get_SV_already_used(sv_id)) {SV_already_used = true; return 0; }//sv region already used
	sf.set_SV_already_used(sv_id);

	//get basic information
	//try to get next SV in the same cluster:
	cluster_used_sv_id.clear();
	SV_cluster_item * c_item = sf.get_cluster(sv_id);
	//cluster has more than one SV
	if(c_item->next_SV_ID != -1){//it has next node
		//get average read mapQ
		//get start and end index for each block
		int best_sv_id = -1; int best_score = 0;
		for(int test_sv_id = sv_id; test_sv_id != -1; test_sv_id = c_item->next_SV_ID, c_item = sf.get_cluster(test_sv_id)){
			sf.set_SV_already_used(test_sv_id);//set the SV to be used
			if(rl.get_sv_region_by_ID(&start_index, &load_read_number, test_sv_id) == false) continue; //not out of range
			if(load_read_number == 0) continue;
			//load read list
			index_st = rl.sort_l + start_index;
			read_buff = rl.read_load_buff;
			cluster_used_sv_id.emplace_back(test_sv_id, load_read_number,rl.sort_l + start_index, rl.read_load_buff, false);
			int count_read_high_mapQ = 0; int read_total_mapQ = 0; int read_total_num = 0; //counters
			for(int i = 0; i < load_read_number; i++){
				bam1_t  * br = read_buff + index_st[i].bam_id;
				int chaining_score = -1; bam_get_num_tag(br, vcf_tag.CS, &chaining_score);
				if(chaining_score == -1) continue;
				if(br->core.qual > 5)
					count_read_high_mapQ++;
				read_total_num++;
				read_total_mapQ += br->core.qual;
			}
			int score = count_read_high_mapQ * 10 + read_total_mapQ + read_total_num * 2;//todo:::
			if(score > best_score){	best_score = score; best_sv_id = test_sv_id;}

		}
		sv_id = best_sv_id;
		//store the final used SV
		for(auto& c_sv_id:cluster_used_sv_id)
			if(c_sv_id.sv_ID == sv_id)
				c_sv_id.is_main = true;

		rl.get_sv_region_by_ID(&start_index, &load_read_number, sv_id);
		if(sv_id == -1)
		{SV_already_used = true; return 0; }
	}

	//get read list, read number stored in [load_read_number], index_st stored the read index list
	//using [bam1_t  * br = read_buff + index_st[i].bam_id;] to get the i-th bam
	index_st = rl.sort_l + start_index;
	read_buff = rl.read_load_buff;

	sv_info = sf.get_sv_info(sv_id);
	sv_region_st = sv_info->st_pos;
	sv_length = sv_info->region_len;
	ori_accept_region1_bg = sv_info->st_pos - normal_read_len;
	ori_accept_region1_ed = sv_info->break_point1_pos;
	ori_accept_region2_bg = sv_info->break_point1_pos - normal_read_len;
	ori_accept_region2_ed = sv_info->ed_pos;

	fprintf(detail_output, "\n\nDetail of %s\t", sv_info->vcf_print_string.c_str());	//print basic SV info
	//with other part of SVs in the cluster
	if(cluster_used_sv_id.size() > 1){
		fprintf(detail_output, "[With other SV in the cluster, we treat as one: ");	//print basic SV info
		for(auto&  c_sv_ID:cluster_used_sv_id ){
			fprintf(detail_output, " %s %s;\t", sf.get_sv_info(c_sv_ID.sv_ID)->vcf_print_string.c_str(), c_sv_ID.is_main?"(MIAN)":"");	//print basic SV info
		}
		fprintf(detail_output, "] ");
	}

	//clear depth_counter
	memset(depth_counter, 0, sv_length * sizeof(uint32_t));
	//clear assembler
	am->clear();
	//fprintf(stderr, "load_read_number, %d\n", load_read_number);

	return load_read_number;
}

bool AssemblyHandler::get_read_filter1(bam1_t  * br){
	//get read information
	int chaining_score = -1;
	bam_get_num_tag(br, vcf_tag.CS, &chaining_score);

	if(chaining_score == -1){
		int read_pos = br->core.pos;
		if(read_pos > ori_accept_region1_bg && read_pos < ori_accept_region1_ed)
			return true;
		if(read_pos > ori_accept_region2_bg && read_pos < ori_accept_region2_ed)
			return true;
		return false;
	}
	else{
		return true;
	}
}

AssemblyHandler::read_score_filter * AssemblyHandler::get_rsf_by_pos(int pos){
	int st_idx = pos - sv_region_st;
	//left and right has additional blocks to store additional read
	int sf_idx = (st_idx >> rsf_shift_offset) + 9;
	if(sf_idx < 0 || sf_idx >= total_filter_block_num){
		return NULL;
	}
	return &(rsf_buff[sf_idx]);
}

void AssemblyHandler::addRead2depthCounter(bam1_t  * br){
	int bam_pos = br->core.pos;
	//get depth
	int read_in_ref_offset = bam_pos - sv_region_st;
	int cigar_l = br->core.n_cigar;
	uint32_t *cigar =  bam_get_cigar(br);
	for(int cigar_ID = 0;cigar_ID < cigar_l; cigar_ID++){
		int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
		int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
		switch(type){
		case 0:	for(int i = 0; i < cigar_len; i++, read_in_ref_offset++) { if(read_in_ref_offset >= 0){depth_counter[read_in_ref_offset]++;}} break;//M
		case 1:	break;//I, do nothing
		case 2:	case 3:	case 4:	read_in_ref_offset += cigar_len;	break;
		default: fprintf(detail_output, "ERROR CIGAR  %d %d ", type, cigar_len);
		}
	}
}

struct SCORE_FILTER{
	static const int SCORE_PASS = 0;
	static const int SMALL_SCORE = 1;
	static const int SAME_SCORE = 2;
	static const int XA_BIGGER_2 = 3;
	static const int XA_2 = 4;
	static const int UNKOWN = 5;

	static const char* get_str(int reason){
		switch(reason){
		case SCORE_PASS: return "SCORE_PASS"; break;
		case SMALL_SCORE: return "SMALL_SCORE"; break;
		case SAME_SCORE: return "SAME_SCORE"; break;
		case XA_BIGGER_2: return "XA_BIGGER_2"; break;
		case XA_2: return "XA_2"; break;
		default: return "UNKOWN"; break;
		}
		return "";
	}
};

#define SCORE_DIFF_L1 35
int readScoreFilter(bam1_t  * br){
	//get read original mapq filter
	//FIR mapQ:40 flag: 80 score: [238, 238, 107]
	//[OA:0,2038276,0,60,M;]
	//[MV:485_0_2037843_1105_INS_svim] [XA:(null)]
	//[RC:0,2038276,0,238,60,60,0,0,237,RNNN0,FNYN0_]126M
	//get map score:
	int score, ori_score;
	bam_get_num_tag(br, vcf_tag.AS, &score);
	bam_get_num_tag(br, vcf_tag.OS, &ori_score);

	if(score < ori_score)
		return SCORE_FILTER::SMALL_SCORE;
	if(score == ori_score)
		return SCORE_FILTER::SAME_SCORE;
	else if(score < ori_score + SCORE_DIFF_L1){
		//RC:Z:15,90270534,0,228,25,25,0,0,590,RNNN0,FNNN0_
		char* RC_tag_char = bam_get_string_tag(br, vcf_tag.RC);
		//get original mapq:
		int str_len_RC = strlen(RC_tag_char);
		char *multy_safty_strtok_ptr = NULL; char * token; int chr_ID; int mapq; int XA_num;
		token = strtok_r(RC_tag_char, ",", &multy_safty_strtok_ptr); 	chr_ID = atoi(token); //chr_ID
		token = strtok_r(NULL, ",", &multy_safty_strtok_ptr); 	//position
		token = strtok_r(NULL, ",", &multy_safty_strtok_ptr); 	//soft left
		token = strtok_r(NULL, ",", &multy_safty_strtok_ptr); 	//score
		token = strtok_r(NULL, ",", &multy_safty_strtok_ptr); 	mapq = atoi(token);//mapq
		token = strtok_r(NULL, ",", &multy_safty_strtok_ptr); 	//mate mapq
		token = strtok_r(NULL, ",", &multy_safty_strtok_ptr); 	XA_num = atoi(token);//XA_number[i]
		for(int i = 0; i < str_len_RC - 1; i++)	if(RC_tag_char[i] == 0)	RC_tag_char[i] = ',';
		if(mapq == 0 && XA_num > 2)
			return SCORE_FILTER::XA_BIGGER_2;
		if(mapq == 0 && XA_num == 2 && chr_ID < 24)
			return SCORE_FILTER::XA_2;
	}
	return SCORE_FILTER::SCORE_PASS;
}

void AssemblyHandler::output_reads(bam1_t  * br, int SCORE_FILTER){
	int read_in_ref_offset = br->core.pos - sv_region_st;
	int cigar_l = br->core.n_cigar;
	uint32_t *cigar =  bam_get_cigar(br);
	for(int i = 0; i < read_in_ref_offset; i++) fprintf(detail_output, "-");
	int seq_i = 0;
	for(int cigar_ID = 0;cigar_ID < cigar_l; cigar_ID++){
		int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
		int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
		switch(type){
		case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {if(read_in_ref_offset >= 0){fprintf(detail_output, "%c", read_str[seq_i]);}read_in_ref_offset++;}	break;//M
		case 1:	seq_i += cigar_len;	break;//I, print nothing
		case 2:	for(int i = 0; i < cigar_len; i++) {if(read_in_ref_offset >= 0) fprintf(detail_output, "-");  read_in_ref_offset++;}	break;//D, print '-'
		case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {if(read_in_ref_offset >= 0) fprintf(detail_output, "N"); read_in_ref_offset++;}	break;//N, print N
		case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {if(read_in_ref_offset >= 0) fprintf(detail_output, "-"); read_in_ref_offset++;}	break;//S, print -
		default: fprintf(detail_output, "ERROR CIGAR  %d %d ", type, cigar_len);
		}
	}
	int bam_position = br->core.pos;
	fprintf(detail_output, "pos %d offset %d %s ", bam_position, bam_position - sv_region_st, SCORE_FILTER::get_str(SCORE_FILTER));
	print_info(br);
	for(seq_i = 0; seq_i < br->core.l_qseq; seq_i++)	fprintf(detail_output, "%c", read_str[seq_i]);
	fprintf(detail_output, "\n");
}

bool AssemblyHandler::read_depth_filter(bam1_t  * br, bool read_from_main_SV, int read_idx){
	//filter 1, position filter
	if(read_from_main_SV){ if(!read_filter[read_idx] == true) return false;}
	else{ if(!get_read_filter1(br)) return false; }//filter 1, position filter

	//read filter 2: depth and score filter
	int c_read_score;bam_get_num_tag(br, vcf_tag.AS, &c_read_score);
	read_score_filter *rsf = get_rsf_by_pos(br->core.pos);
	if(rsf == NULL || c_read_score < rsf->filter_score){return false; }

	//read filter 3: read from main SV is enough, no need for more reads from other SV in the SV cluster
	if(!read_from_main_SV && rsf->passed_read_num >= max_read_in_block){return false;}

	//read filter 4: random filter
	if(c_read_score == rsf->filter_score && rsf->passed_read_num >= (max_read_in_block * 2)){
		if(rand() % rsf->passed_read_num > max_read_in_block){ return false;}
	}
	return true;
}

void AssemblyHandler::addReads2assembler(){
	//step1: get reference block minus score
	//we take 64 base pair as a reference block
	//when the block is over depth, (> 100 X coverage, or read_num * read_len > 64 * max_coverage), we count the distribution of score
	//and discard some low score read to let the coverage be less then 60;
	ab_block_size = (sv_id == 1)?10000: ab_block_size_option;//set special ab length for SV 0

	ab_size = (sv_length / ab_block_size) + 1;
	if((int)ab_list.size() < ab_size)
		ab_list.resize(ab_size);
	//clear
	for(int i = 0; i < ab_size; i++)
		ab_list[i].init();

	total_filter_block_num = (sv_length >> rsf_shift_offset) + 21;
	//reset filter block
	rsf_buff.resize(total_filter_block_num);
	for(auto & c_sf : rsf_buff){
		c_sf.read_idx_ed = 0;
		c_sf.read_idx_bg = MAX_int32t;
		c_sf.filter_score = 0;
		c_sf.passed_read_num = 0;
	}

	//get start and end index for each block
	for(int i = 0; i < load_read_number; i++){
		bam1_t  * br = read_buff + index_st[i].bam_id;
		read_score_filter *rsf = get_rsf_by_pos(br->core.pos);
		if(rsf == NULL)	 continue;
		rsf->read_idx_bg = MIN(i, rsf->read_idx_bg);
		rsf->read_idx_ed = MAX(i, rsf->read_idx_ed);
	}

	//set filter1 for each reads and get score filter for each block
	read_filter.resize(load_read_number);//to show whether a read will be used in assembly
	for(auto & c_sf : rsf_buff){
		//get read number in the block
		int read_num = 0;
		for(int i = c_sf.read_idx_bg; i < c_sf.read_idx_ed; i++){
			bam1_t  * br = read_buff + index_st[i].bam_id;
			bool pass_filter = get_read_filter1(br);
			read_filter[i] = pass_filter;
			if(pass_filter) read_num++;
		}
		//when read is not over-depth, not set score filter
		if(read_num <= max_read_in_block){
			c_sf.filter_score = 0; continue;
		}

		//count score when over depth
		score_counter.clear();
		for(int i = c_sf.read_idx_bg; i < c_sf.read_idx_ed; i++){
			if(!read_filter[i]) continue;
			bam1_t  * br = read_buff + index_st[i].bam_id;
			//get score and store score in map
			int c_read_score;
			bam_get_num_tag(br, vcf_tag.AS, &c_read_score);
			std::map<int, int>::iterator sc_it = score_counter.find(c_read_score);
			if(sc_it != score_counter.end())
				sc_it->second++;
			else
				score_counter[c_read_score] = 1;
		}
		//score is auto sorted from small to big
		int total_passed_read_count = 0;
		for(auto rit = score_counter.rbegin(); rit != score_counter.rend(); rit++){
			total_passed_read_count += rit->second;
			if(total_passed_read_count > max_read_in_block){
				c_sf.filter_score = rit->first; c_sf.passed_read_num = total_passed_read_count; break;
			}
		}
	}

	//adding read into assembly
	for(int i = 0; i < load_read_number; i++){
		bam1_t  * br = read_buff + index_st[i].bam_id;
		if(read_depth_filter(br, true, i) == false)	continue;
		//get read data
		get_bam_seq(0, br->core.l_qseq, read_str, br);
		int score_filter = readScoreFilter(br);

		addRead2depthCounter(br);//get depth
		if(output_original_read){ output_reads(br, score_filter); }		//output read
		//get assembly block ID:
		if(score_filter == SCORE_FILTER::SCORE_PASS ){
			int ab_block_id = (br->core.pos - sv_region_st)/ ab_block_size;
			ab_block_id = MAX(0, ab_block_id); ab_block_id = MIN(ab_size - 1, ab_block_id);
			ab_list[ab_block_id].add_read(br, read_str, sv_region_st, true);
		}
	}

	 //adding reads in other sv into cluster
	for(auto & c_sv_in_cluster: cluster_used_sv_id){
		if(c_sv_in_cluster.is_main) continue;

		for(int i = 0; i < c_sv_in_cluster.load_read_number; i++){
			bam1_t  * br = c_sv_in_cluster.read_buff + c_sv_in_cluster.index_st[i].bam_id;
			if(read_depth_filter(br, true, i) == false)	continue;

			get_bam_seq(0, br->core.l_qseq, read_str, br);
			//get read original mapq filter
			//FIR mapQ:40 flag: 80 score: [238, 238, 107]
			//[OA:0,2038276,0,60,M;]
			//[MV:485_0_2037843_1105_INS_svim] [XA:(null)]
			//[RC:0,2038276,0,238,60,60,0,0,237,RNNN0,FNYN0_]126M
			int score_filter = readScoreFilter(br);
			addRead2depthCounter(br);//get depth
			if(output_original_read){ output_reads(br, score_filter); }		//output read
			if(score_filter == SCORE_FILTER::SCORE_PASS ){
				//get assembly block ID:
				int ab_block_id = (br->core.pos - sv_region_st)/ ab_block_size;
				ab_block_id = MAX(0, ab_block_id); ab_block_id = MIN(ab_size - 1, ab_block_id);
				//get read data
				ab_list[ab_block_id].add_read(br, read_str, sv_region_st, false);
			}
		}
	}
}

void printContigInfo(FILE * output, int ori_read_number, AssembledContig & contig, int contig_ID, std::vector<ASS_add_reads> &ass_read_list){
	//basic part
	if(true){
		fprintf(output,
				"contig_ID: [%d] "
				"word length: [%d] "
				"CONTIG size: [%ld] "
				"supportReads [%ld] "
				"ending_reason: [%d %d]"
				"new_support_read[%d] "
				,
				contig_ID,
				contig.wordLength, contig.seq.size(), contig.supportReads.size(), contig.ending_reason[0], contig.ending_reason[1], contig.new_support_read);
	}
	//supplementary information
	if(false){
		fprintf(output,
				"contig_begin_offset: [%d] "
				"seedCount: [%d] "
				"rejectReads[%ld] "
				,
				contig.ass_begin_offset_in_contig,
				contig.seedReadCount,  contig.rejectReads.size());
		fprintf(output, "supportReads: ");
		for (auto &r : contig.supportReads){
			if((int)r < ori_read_number){
				fprintf(output, "[ %d %d %d]\t", r,
						ass_read_list[r].read_in_ref_offset,
						ass_read_list[r].br->core.qual);
			}
		}
	}

	fprintf(output, "\n");
}

#define MAX_INDEL_LEN 80
#define HAED_MIN_MATCH_BASE 20
int get_var(Contig_aligner &ca, int suggest_st_pos, uint16_t * contig_depth, int ass_part, int contig_ID, VI_list &vi_l, std::vector<GlobalDepthItem> &global_depth){
	//get variations
	uint8_t *tseq = ca.tseq;
	uint8_t *qseq = ca.qseq;
	int output_index = 0;
	int seq_i = 0;
	output_index += suggest_st_pos;
	uint32_t* bam_cigar =ca.ez.cigar;
	int cigar_number = ca.ez.n_cigar;
	//ignore the first 8 match bases, the first 8 base of a contig not used to calculate depth , nor get variation, this function is used to discard some false positive deletions
	int match_base = 0;int ref_pos_enough_match_base = 0; bool finish_head_check = false;
	for(int cigar_ID = 0;cigar_ID < cigar_number; cigar_ID++){
		int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
		int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
		switch(type){
		case 0:
		for(int i = 0; i < cigar_len; i++, seq_i++, output_index++){
			if(!finish_head_check){
				if(global_depth[output_index].ref_base == qseq[seq_i]) match_base++; else match_base--;
				if(match_base >= HAED_MIN_MATCH_BASE) {ref_pos_enough_match_base = output_index; finish_head_check = true;}
				continue;
			}
			global_depth[output_index].set_base(qseq[seq_i], ass_part, contig_depth[seq_i]);
			if(qseq[seq_i] != tseq[output_index])
				vi_l.add_data(tseq + output_index, 1, qseq + seq_i, 1, output_index, seq_i, contig_depth[seq_i], ass_part, contig_ID, VAR_TYPE::SNP);
		}	break;//M
		case 1://Insertion
		if(cigar_ID != 0 && cigar_ID != cigar_number - 1 && cigar_len < MAX_INDEL_LEN){//not at the begin or end of a cigar list
			if(!finish_head_check)
				match_base -= 2;
			else{
				vi_l.add_data(tseq + output_index, 1, qseq + seq_i, cigar_len + 1, output_index, seq_i, contig_depth[seq_i], ass_part, contig_ID, VAR_TYPE::INS);
				global_depth[output_index].set_base(5, ass_part, contig_depth[seq_i] * 2);
			}
		}
		seq_i += cigar_len;
		break;
		case 2://Deletion
		if(cigar_ID != 0 && cigar_ID != cigar_number - 1 && cigar_len < MAX_INDEL_LEN){//not at the begin or end of a cigar list
			if(!finish_head_check){ match_base -= 2;}
			else{
				vi_l.add_data(tseq + output_index, cigar_len + 1, qseq + seq_i, 1, output_index, seq_i, contig_depth[seq_i], ass_part, contig_ID, VAR_TYPE::DEL);
				for(int i = 0; i < cigar_len; i++, output_index++)
					global_depth[output_index].set_base(4, ass_part, contig_depth[seq_i]);
				output_index -= cigar_len;
			}
		}
		output_index += cigar_len;
		break;
		case 3:	case 4:
		output_index +=	cigar_len; break;//S
		default: fprintf(stderr, "ERROR CIGAR  %d %d ", type, cigar_len);
		}
	}
	return ref_pos_enough_match_base;
}

static int analysis_region(std::vector<GlobalDepthItem> & di, int st_pos, int ed_pos, int &ins_number, int &del_number, int &SNP_number, double & ave_depth, int & min_depth){
	//number count
	int blank_num = 0;
	int total_depth = 0;
	min_depth = MAX_int32t;
	int depth_counter = 0;
	xassert(st_pos >= 0, ""); xassert(ed_pos <= (int)di.size(), "");
	for(int i = st_pos; i < ed_pos; i++){
		int ei = di[i].event_info();
		switch(ei){
		case 0: case 2: break;//do nothing
		case 1: blank_num++; break;
		case 3: case 4: case 5: case 6: SNP_number ++; break;
		case 7: del_number++; break;
		case 8: ins_number++; break;
		default: di[i].print(stderr); fprintf(stderr, "di.size %ld ;idx: %d ei: %d\n", di.size(), i, ei); xassert(0, "Wrong event"); break;
		}
		if(ei != 7 && ei != 1){//not blank, not deletion
			depth_counter ++;
			total_depth += di[i].total_depth;
			min_depth = MIN(min_depth, di[i].total_depth);
		}
	}
	ave_depth = (depth_counter == 0)?0:((float)total_depth/depth_counter);
	if(min_depth == MAX_int32t) min_depth = 0;
	return blank_num;
}

struct CHECK_SV_FAIL_REASON{
	static const int filter_pass = -1;
	static const int bp1_uncovered = 0;
	static const int bp2_uncovered = 1;
	static const int ins_uncovered = 2;
	static const int del_length_not_enough = 3;
	static const int ins_length_not_enough = 4;
	static const int del_depth_change_sharply = 5;
	static const int low_total_depth = 6;
	static const int wrong_sv_type = 7;

	static const char * show_reason(int reason_ID){
		switch(reason_ID){
		case -1: return "filter_pass"; break;
		case 0: return "bp1_uncovered"; break;
		case 1: return "bp2_uncovered"; break;
		case 2: return "ins_uncovered"; break;
		case 3: return "del_length_not_enough"; break;
		case 4: return "ins_length_not_enough"; break;
		case 5: return "del_depth_change_sharply"; break;
		case 6: return "low_total_depth"; break;
		case 7: return "wrong_sv_type"; break;
		default: return "Unknown"; break;
		}
	}
};

static void analysis_full_data(FILE * summary_output, FILE * output, std::vector<GlobalDepthItem> & di, int sv_length, int break_point1, int break_point2, VI_list &vi_l, SV_chr_info * sv_info, SV_ref_sequence *sf, double ave_read_depth, float depth_bp1, float depth_bp2){
	int ins_number_bp1 = 0; int del_number_bp1 = 0; int SNP_number_bp1 = 0; double depth_bp1_contig = 0; int min_depth_bp1 = 0;
	int ins_number_bp2 = 0; int del_number_bp2 = 0; int SNP_number_bp2 = 0; double depth_bp2_contig = 0; int min_depth_bp2 = 0;
	int ins_number_ins = 0; int del_number_ins = 0; int SNP_number_ins = 0; double depth_ins_contig = 0; int min_depth_ins = 0;
	//detail around the break point1
	bool has_ins_part = (break_point2 > break_point1 + 10); int ins_part_len = break_point2 - break_point1;
	int bp_chech_region_len;
	if(has_ins_part) bp_chech_region_len = 10; else bp_chech_region_len = 20; //10 for insertion; 20 for deletion or other type

	int bp1_blank_num = analysis_region(di, break_point1 - bp_chech_region_len , break_point1 + bp_chech_region_len, ins_number_bp1, del_number_bp1, SNP_number_bp1, depth_bp1_contig, min_depth_bp1);
	int bp2_blank_num = analysis_region(di, break_point2 - bp_chech_region_len , break_point2 + bp_chech_region_len, ins_number_bp2, del_number_bp2, SNP_number_bp2, depth_bp2_contig, min_depth_bp2);

	int ins_blank_num = 0;
	if(has_ins_part)
		ins_blank_num = analysis_region(di, break_point1, break_point2, ins_number_ins, del_number_ins, SNP_number_ins, depth_ins_contig, min_depth_ins);

	int fail_reason = CHECK_SV_FAIL_REASON::filter_pass;
	if(has_ins_part){//insertion:
		//maybe only one insertion site is found, when this occurred， and at same time the insertion part is covered over 50%, the insert string will be set in the covered break point
		if(!(sv_info->sv_type[0] == 'I' ||(sv_info->sv_type[0] == 'D' && sv_info->sv_type[1] == 'U')))
			fail_reason = CHECK_SV_FAIL_REASON::wrong_sv_type;
		else if(bp1_blank_num > 0 && bp2_blank_num > 0)
			fail_reason = CHECK_SV_FAIL_REASON::bp1_uncovered;
		else{
			if((bp1_blank_num > 0 && (ins_number_bp2 + del_number_bp2) > 0) || (bp2_blank_num > 0 && (ins_number_bp1 + del_number_bp1) > 0)){
				fail_reason = CHECK_SV_FAIL_REASON::bp1_uncovered;
			}
			else if(ins_blank_num > 0.5*ins_part_len){  fail_reason = CHECK_SV_FAIL_REASON::ins_uncovered; }
			else if(del_number_ins + ins_blank_num + 30 > ins_part_len){fail_reason = CHECK_SV_FAIL_REASON::ins_length_not_enough;}
		}
	}else{//deletion
		if(!(sv_info->sv_type[0] == 'D' && sv_info->sv_type[1] == 'E'))
			fail_reason = CHECK_SV_FAIL_REASON::wrong_sv_type;
		else if(bp1_blank_num > 0 || bp2_blank_num > 0)
			fail_reason = CHECK_SV_FAIL_REASON::bp1_uncovered;
		else if(depth_bp1_contig != 0 && min_depth_bp1 * 2 < depth_bp1_contig){
			fail_reason = CHECK_SV_FAIL_REASON::del_depth_change_sharply;
		}
		else{
			int ins_len = 0;
			for(auto & vi:vi_l.vi)	if(vi.ref_position > break_point1 - 10 && vi.ref_position < break_point2 + 10 && di[vi.ref_position].event_info() == 8)	ins_len += vi.alt.size();
			int vcf_deletion_length = sv_info->getDeletionLen();
			xassert(vcf_deletion_length > 0, "");
			if(ins_len + 30 > vcf_deletion_length){fail_reason = CHECK_SV_FAIL_REASON::del_length_not_enough; }
		}
	}

	int minReadDepth = ave_read_depth * 0.1;
	minReadDepth = MAX(minReadDepth, 3);//at least 3 reads support
	int read_depth = (depth_bp1 + depth_bp2) / 2;//at least 3 perfect reads support the SV
	if(fail_reason == CHECK_SV_FAIL_REASON::filter_pass && read_depth < minReadDepth){fail_reason = CHECK_SV_FAIL_REASON::low_total_depth; }
	int ass_depth = (depth_bp1_contig + depth_bp2_contig)/2;
	if(fail_reason == CHECK_SV_FAIL_REASON::filter_pass && ass_depth < minReadDepth){fail_reason = CHECK_SV_FAIL_REASON::low_total_depth; }

	//show results
	fprintf(output, "Final result: %s reason %s\t demoVCF:{\t", (fail_reason < 0)?"PASS":"FAIL", CHECK_SV_FAIL_REASON::show_reason(fail_reason));
	if(fail_reason  == CHECK_SV_FAIL_REASON::filter_pass){
		fprintf(output, "with vcf\n");

		std::vector<char> ref; ref.clear();//fatal get deletion string from original reference files
		std::vector<char> alt; ref.clear();//generate the alt string for SV analysis
		if(has_ins_part){ //insertion
			//get ALT string
			//collect insertion actions
			std::vector<VariationInfo> insertionActions;
			for(auto & c_vi : vi_l.vi)
				if(c_vi.var_type == VAR_TYPE::INS)
					insertionActions.emplace_back(c_vi);
			int current_ins_index = insertionActions.size() - 1;
			//reverse insert deletion and insertions into alt string
			char * insert_string ; int insert_string_size;
			for(int c_pos = break_point2 - 1; c_pos >= break_point1; c_pos--){
				int ei = di[c_pos].event_info();
				switch(ei){
				case 0: case 2:  case 3: case 4: case 5: case 6:
					alt.emplace_back("ACGT"[di[c_pos].max_base]); break;
				case 1: case 7:  break; //blank and deletion: do nothing
				case 8: //search ins string
					xassert(current_ins_index != -1, "");
					for(int ins_idx = current_ins_index; ins_idx > -1; ins_idx--){
						if(c_pos != insertionActions[ins_idx].ref_position) continue;
						else {current_ins_index = ins_idx; break;}
					}
					insert_string_size = insertionActions[current_ins_index].alt.size();
					insert_string  = &(insertionActions[current_ins_index].alt[0]);
					for(int base_idx = insert_string_size - 1; base_idx > 0; base_idx--) //not insert the first base
						alt.emplace_back(insert_string[base_idx]);
					current_ins_index -= 1;//skip this base
					break;
				default: xassert(0, "Wrong event"); break;
				}
			}
			//reverse the alt string:
			int alt_size = alt.size();
			for(int base_idx = alt_size / 2; base_idx >= 0 ; base_idx --){
				uint8_t tmp_base = alt[base_idx];
				alt[base_idx] = alt[alt_size - base_idx - 1];
				alt[alt_size - base_idx - 1] = tmp_base;
			}
			//add two additional events at break point BASE:
			//event at bp 1
			//event at bp 2
		}else{
			xassert(sv_info->sv_type[0] == 'D', "");//DEL
			//get deletion start base
			uint64_t sv_st_base = sv_info->break_point1_pos;
			uint64_t sv_ed_base = sv_info->break_point2_pos;

			//get additional event at bp 1 and bp 2 to modify bp 1 and bp2
			//event at bp 1
			//event at bp 2

			char * ori_ref_base = sf->get_original_ref_seq(sv_info->chr_ID, sv_st_base, sv_ed_base);
			//get ref string
			xassert(ori_ref_base != NULL, "");
			for(int i = 0; ori_ref_base[i] != 0; i++ )
				ref.emplace_back(ori_ref_base[i]);
		}

		//todo:: adding additional base at begin of ref and alt
		bool begin_position_is_modified = false; //todo::
		char additional_base_ref_alt = "ACGT"[di[break_point1 - 1].ref_base];
		if(begin_position_is_modified){ //todo::
			//todo::load base from original reference
		}
		//1	11028	pbsv.DEL.0	GTGTGTTGCAGGAGCAAAGTCGCACGGCGCCGGGCT	G	.	PASS	SVTYPE=DEL;END=11063;SVLEN=-35;SVANN=TANDEM	GT:AD:DP	0/1:33,12:45
		bool low_depth = (depth_bp1_contig + depth_bp2_contig < 5);
		//print simple VCF:
		int st_pos = sv_info->break_point1_pos; //todo:: polish st_pos
		int length = alt.size() - ref.size() + 1; //todo: get sv length
		int endPos = st_pos + ref.size();
		std::string genotype = "./.";//todo:: get genotype
		double c_ave_depth = (depth_bp1 + depth_bp2) / 2;
		bool isHeterozygous = ((c_ave_depth) < ave_read_depth *0.45) ? true: false;
		int allele_number = 2;
		int allele_read_depth[2] = {0};//todo::
		allele_read_depth[0] = (depth_bp1 + depth_bp2) / 2;

		int depth = 0;
		for(int i = 0; i < allele_number; i++)
			depth += allele_read_depth[i];
		fprintf(summary_output, "%s\t%d\t%s\t", sf->get_chr_name_by_ID(sv_info->chr_ID), st_pos, sv_info->vcf_print_string.c_str());
		//ref and alt seq:
		ref.emplace_back(0); fprintf(summary_output, "%c%s\t", additional_base_ref_alt, &(ref[0]));
		alt.emplace_back(0); fprintf(summary_output, "%c%s\t", additional_base_ref_alt, &(alt[0]));

		fprintf(summary_output, ".\t%s\tSVTYPE=%s;END=%d;SVLEN=%d\t" "GT:DP\t" "%s:", (low_depth)?"LOW_DEPTH":"PASS" ,sv_info->sv_type.c_str(), endPos, length , (isHeterozygous)?"0/1":"1/1");
		fprintf(summary_output, "%d,%d,%d,%d", (int)depth_bp1, (int)depth_bp2, (int)depth_bp1_contig, (int)depth_bp2_contig);
		//for(int i = 0; i < allele_number; i++)
			//fprintf(summary_output, "%d,", allele_read_depth[i]);
		//fprintf(summary_output,	":%d", depth);
		fprintf(summary_output,	"\n");
	}else
		fprintf(output, "No vcf \n");
}

void AssemblyHandler::show_assembly_block_info(int ab_idx ,size_t contig_size){
	//print part basci info, include: read information and contig information
	fprintf(detail_output,	"part [%d], contigs.size[%ld]\t read list (format: [offset mapq])", ab_idx, contig_size);
	for (ASS_add_reads &r : ab_list[ab_idx].ass_read_list){
		fprintf(detail_output, "[%d %d]\t",
				r.read_in_ref_offset,
				r.br->core.qual);
	}
	fprintf(detail_output,	"\n");
}

void printAction(std::vector<AddReadAction> &actions, FILE * detail_output, int ori_read_number, std::set<int> &remove_read_set, std::vector<ASS_add_reads> &ass_read_list){
	for (auto &ca : actions)
		if(ca.read_ID < ori_read_number && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
			ca.print(detail_output, ass_read_list[ca.read_ID].is_from_main_SV);
	fprintf(detail_output, "\n");
}

bool AssemblyHandler::get_suggention_alignment_position_list(AssembledContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_add_reads> &ass_read_list ){

	const char *contig_seq = contig.seq.c_str();
	int contig_seq_len = contig.seq.size();
	//step1: get all the read positions at assembled CONTIG
	//set read position for each actions
	remove_read_set.clear();
	for (auto &ca : contig.actions)
		if(ca.read_ID < ori_read_number && ass_read_list[ca.read_ID].is_from_main_SV)
			ca.set_read_pos(read_list[ca.read_ID], contig.seq, contig.ass_begin_offset_in_contig, contig.wordLength, remove_read_set, ass_read_list[ca.read_ID].read_in_ref_offset);

	//step2: find the suggestion alignment position of CONTIG at reference from each read position
	suggest_st_pos_map.clear();
	//get contig coverage
	int align_suggention_total_used_read = 0;
	if((int)contig_depth.size() < contig_seq_len)	contig_depth.resize(contig_seq_len);
	memset(&(contig_depth[0]), 0, contig_seq_len*sizeof(uint16_t));

	int new_score_is_bigger_read_cn = 0;
	int new_score_is_smaller_read_cn = 0;

	int new_score_is_far_bigger_read_cn = 0;
	int new_score_is_not_far_bug_read_cn = 0;

	for (auto &ca : contig.actions){
		if(ca.read_ID < ori_read_number && ass_read_list[ca.read_ID].is_from_main_SV){
			if(remove_read_set.find(ca.read_ID) != remove_read_set.end()) continue;
			if(ass_read_list[ca.read_ID].new_score_is_bigger){ new_score_is_bigger_read_cn ++;} else {new_score_is_smaller_read_cn++;}
			if(ass_read_list[ca.read_ID].new_score_is_far_bigger){ new_score_is_far_bigger_read_cn ++;} else {new_score_is_not_far_bug_read_cn++;}
			//set coverage:
			const char* read_seq = read_list[ca.read_ID].c_str();
			int st_pos_ref = ca.position_in_contig - contig.ass_begin_offset_in_contig - ca.position_read;
			int ed_pos_ref = st_pos_ref + read_list[ca.read_ID].size();
			int st_pos_read = 0; if(st_pos_ref < 0){	st_pos_read -= st_pos_ref;st_pos_ref = 0;}
			ed_pos_ref = MIN(contig_seq_len, ed_pos_ref);
			ca.wrong_base = 0;
			for(int i = st_pos_ref; i < ed_pos_ref && ca.wrong_base <= 8; i++, st_pos_read++){
				if(contig_seq[i] != read_seq[st_pos_read]) ca.wrong_base++;
			}
			if(ca.wrong_base <= 8){
				st_pos_read = 0;
				for(int i = st_pos_ref; i < ed_pos_ref; i++, st_pos_read++){
					if(contig_seq[i] == read_seq[st_pos_read]) contig_depth[i] ++;
				}
				align_suggention_total_used_read++;
				std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
				if(it != suggest_st_pos_map.end())	it->second++;
				else suggest_st_pos_map[ca.suggest_contig_offset_in_ref] = 1;
			}
		}
	}
	//if no result is there, simple selected the first read position
	if(suggest_st_pos_map.empty())	suggest_st_pos_map[contig.actions[0].suggest_contig_offset_in_ref] = 1;

	//step3: find the max suggest position
	//get max suggest alignment position
	//try simple merge:
	int max_suggention = 0; int max_count = 0;
	for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++)
		if(max_count < it->second) {max_suggention = it->first; max_count = it->second;}

	//step4: find the suggestion list
	//if a suggestion covers most of read, [unique_alignment_pos] is true and use it as alignment suggestion, otherwise store all suggestions in a list [suggent_pos_list]
	// when no result or result > 1: [unique_alignment_pos] = false
	//suggest_coverage only use when suggent_pos_list.size() > 1;
	float suggest_coverage = 0;
	suggent_pos_list.clear();
	if(max_count * 2 <= align_suggention_total_used_read || align_suggention_total_used_read == 0){
		for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++)
			if(it->second*2 >= max_count || it->second > 5)
				suggent_pos_list.emplace_back(it->first);
		if(suggent_pos_list.size() > 1){ //== 1
			suggest_coverage = (float)contig_seq_len/(suggent_pos_list.back() + contig_seq_len - suggent_pos_list[0]);
			xassert(suggest_coverage <= 1, "");
		}
	}
	else
		suggent_pos_list.emplace_back(max_suggention);

	//step5: modify the depth for multy-suggestion positions
	//coverage modify
	if(suggest_coverage != 0)
		for(int i = 0; i < contig_seq_len; i++)
			contig_depth[i] = contig_depth[i] * suggest_coverage;
	//delete low coverage tail
	int low_cov_tail_len = 0; int low_cov_tail_len_max_search = MIN(20, contig_seq_len);
	for(; low_cov_tail_len < low_cov_tail_len_max_search && contig_depth[contig_seq_len - low_cov_tail_len - 1] < 2; low_cov_tail_len++);
	if(low_cov_tail_len > 0){
		contig.seq.erase(contig_seq_len - low_cov_tail_len, low_cov_tail_len);
		contig_seq_len = contig.seq.size();
	}
	//set MIN depth:
	for(int i = 0; i < contig_seq_len; i++)	contig_depth[i] = MAX(contig_depth[i], 2);

	//print suggestion alignment range information
	if(output_assembly_contig){
		//all possible ranges
		for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++) fprintf(detail_output, "[SUG: %d CN: %d]\t", it->first, it->second);
		fprintf(detail_output, "\n");
		{fprintf(detail_output, "new_score_is_bigger_read_cn [%d %d]", new_score_is_bigger_read_cn, new_score_is_smaller_read_cn); }
		//possible true ranges
		if(suggent_pos_list.empty()) 			{fprintf(detail_output, "NO suggestion\n"); xassert(0, "NO suggestion");}
		else if(suggent_pos_list.size() == 1) 	{fprintf(detail_output, "UNIQUE: [MAX_SUG: %d]\n", max_suggention);}
		else									{fprintf(detail_output, "MULTY:  [COV: %f]\n", suggest_coverage); for(auto sp: suggent_pos_list){ fprintf(detail_output, "[SG: %d]\t", sp);}}
		if(false) printAction(contig.actions, detail_output, ori_read_number, remove_read_set, ass_read_list); //print actions
	}

	if(new_score_is_bigger_read_cn > 1) //at least 66% read has far bigger score
		return true;
	else return false;
	//else if(new_score_is_far_bigger_read_cn > new_score_is_not_far_bug_read_cn) //or at least 33% read has far bigger score
	//	return true;
	//else
	//	return false;
}

void print_contig_aln_ori_cigar(ksw_extz_t &ez, FILE * detail_output, int suggest_contig_st_pos){
	uint32_t* bam_cigar = ez.cigar;
	int cigar_len = ez.n_cigar;
	fprintf(detail_output, "\tOriginal cigar : ");
	for(int i = 0; i < cigar_len; i++){
		fprintf(detail_output, "%d%c", (bam_cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(bam_cigar[i] & BAM_CIGAR_MASK)]);
	}
	fprintf(detail_output, " [ori_ref_st%d]", suggest_contig_st_pos);
}

//using "bin_contig" to store contig sequence
void AssemblyHandler::contig_alignment_and_get_var(int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len){
	//alignment contig into reference
	int suggest_ref_ed_pos = suggest_ref_st_pos + contig_seq_len;
	int suggest_contig_st_pos = 0;
	suggest_ref_st_pos -= 15;
	if(suggest_ref_st_pos < 0){
		if(suggest_ref_st_pos < -15) suggest_contig_st_pos = (- suggest_ref_st_pos - 30);
		suggest_ref_st_pos = 0;
	}
	suggest_ref_st_pos = MAX(0, suggest_ref_st_pos);//get 15 addition alignment region at begin
	suggest_ref_ed_pos += 60; suggest_ref_ed_pos = MIN(suggest_ref_ed_pos, sv_length);//get 60 additional alignment region at end
	bool contig_out_range = (suggest_ref_ed_pos < suggest_ref_st_pos + 20 || suggest_contig_st_pos > contig_seq_len); //when true region length less then 20 bp, it is out of range
	if(output_assembly_contig) fprintf(detail_output, "[aln_ref_region %d ~ %d ] [aln contig st %d ] %s", suggest_ref_st_pos , suggest_ref_ed_pos, suggest_contig_st_pos, contig_out_range?("[ contig_out_range]"):"");
	if(!contig_out_range){
		//print alignment range
		ca.align_non_splice(&(bin_contig[0]) + suggest_contig_st_pos, contig_seq_len - suggest_contig_st_pos,
				suggest_ref_st_pos, suggest_ref_ed_pos);
		if(output_assembly_contig) print_contig_aln_ori_cigar(ca.ez, detail_output, suggest_contig_st_pos);

		uint32_t n_cigar = ca.ez.n_cigar;
		suggest_ref_st_pos += cigar_adjust(&n_cigar, ca.ez.cigar, false, 15);
		ca.ez.n_cigar = n_cigar;
		//fprintf(detail_output, "[aln_ref_region %d ~ %d ] [aln contig st %d ] %s", suggest_ref_st_pos , suggest_ref_ed_pos, suggest_contig_st_pos, contig_out_range?("[ contig_out_range]"):"");
		int ref_pos_enough_match_base = get_var(ca, suggest_ref_st_pos, &contig_depth[0] + suggest_contig_st_pos, ab_idx, contig_ID, vi_l, global_depth);	//get variations
		if(output_assembly_contig){ ca.printf_alignment_detail(detail_output, suggest_ref_st_pos, &contig_depth[0] + suggest_contig_st_pos, ref_pos_enough_match_base);
		fprintf(detail_output, "ref_pos_enough_match_base: %d\t", ref_pos_enough_match_base); }
	}
}

extern uint8_t charToDna5n[];
void AssemblyHandler::store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig){
	int contig_seq_len = contig_string.size();
	const char * contig_seq = contig_string.c_str();
	xassert(nullptr != contig_seq, "");
	//store bin contig
	bin_contig.resize(contig_seq_len);
	for (int i = 0; i < contig_seq_len; ++i)
		bin_contig[i] = charToDna5n[(uint8_t)contig_seq[i]];
}


#define MIN_NEW_SUPPORT_READ 2
void AssemblyHandler::getContig(){
	//running assembly
	//new code:
	vi_l.clear();
	//init, set reference for Contig_aligner
	ca.setRef(sv_id);
	//print reference info
	if(output_assembly_contig) ca.printf_ref(detail_output);
	//global depth
	//clear global depth list
	if((int)global_depth.size() < sv_length) global_depth.resize(sv_length);
	memset(&(global_depth[0]), 0, sizeof(GlobalDepthItem) * sv_length);

	for(int i = 0; i  < sv_length; i++)
		global_depth[i].ref_base = ca.tseq[i];

	//assembly one region part each time
	for(int ab_idx = 0; ab_idx < ab_size; ab_idx++){
		ab_list[ab_idx].run_assembly(am);
		if(output_assembly_contig) show_assembly_block_info(ab_idx ,ab_list[ab_idx].contigs.size());//show results

		//get results
		std::vector<std::string> &read_list = ab_list[ab_idx].reads;
		std::vector<ASS_add_reads> &ass_read_list = ab_list[ab_idx].ass_read_list;
		int ori_read_number = ass_read_list.size();
		int contig_ID = -1;
		auto &contigs = ab_list[ab_idx].contigs;

		//handle each contig
		for (auto &contig : contigs) {
			contig_ID++;
			if(output_assembly_contig) printContigInfo(detail_output, ori_read_number, contig, contig_ID, ass_read_list);

			if(contig_ID != 0 && (contig.new_support_read <= MIN_NEW_SUPPORT_READ && contig.wordLength < 100))
			{ fprintf(detail_output, "This CONTIG is discarded, reason: 'MIN_NEW_SUPPORT_READ'\n"); continue; }

			//get suggestion alignment position
			bool enough_read_has_bigger_score = get_suggention_alignment_position_list(contig, ori_read_number, read_list, ass_read_list );
			if(!enough_read_has_bigger_score) { fprintf(detail_output, "Not enough_read_has_bigger_score than original aligner\n"); }
			//for each suggestion position:
			store_bin_contig(contig.seq, bin_contig);
			for(int suggest_ref_st_pos : suggent_pos_list)
				contig_alignment_and_get_var(ab_idx, contig_ID, suggest_ref_st_pos, contig.seq.size());
		}
	}

	//
	for(int i = 0; i  < sv_length; i++)
		global_depth[i].add_tmp_data();
	vi_l.sort_merge(detail_output, global_depth);

	fprintf(detail_output, "\n");
	GlobalDepthItem::print_full_data(detail_output, global_depth, sv_length, o->edge_len, sv_length - o->edge_len);

	//get total read coverage in full SV reference, normally, the value is higher than the contig coverage
	int break_point1 = o->edge_len; int break_point2= sv_length - o->edge_len;
	bool has_ins_part = (break_point2 > break_point1 + 10);
	int bp_chech_region_len;
	if(has_ins_part) bp_chech_region_len = 10; else bp_chech_region_len = 20; //10 for insertion; 20 for deletion or other type
	float depth_bp1 = averageDepth(break_point1 - bp_chech_region_len , break_point1 + bp_chech_region_len);
	float depth_bp2 = averageDepth(break_point2 - bp_chech_region_len , break_point2 + bp_chech_region_len);
	analysis_full_data(summary_output, detail_output, global_depth, sv_length, break_point1, break_point2, vi_l, sv_info, ca.sf, ave_read_depth, depth_bp1, depth_bp2);
}

void AssemblyHandler::print_seq(bam1_t  * br){
	int bam_pos = br->core.pos;
	int output_index = 0;
	int sv_region_st_tmp = sv_region_st;
	if((int)bam_pos >= sv_region_st_tmp){
		while(sv_region_st_tmp++ < bam_pos){
			fprintf(detail_output, "-"); output_index++;
		}

		int cigar_l = br->core.n_cigar;
		uint32_t *cigar =  bam_get_cigar(br);

		int seq_i = 0;
		for(int cigar_ID = 0;cigar_ID < cigar_l; cigar_ID++){
			int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(detail_output, "%c", read_str[seq_i]); depth_counter[output_index]++; output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(detail_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(detail_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(detail_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(detail_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
	}

	fprintf(detail_output, "pos %d offset %d  ", bam_pos, bam_pos - sv_region_st);

	print_info(br);
	int seq_i = 0;
	for(; seq_i < br->core.l_qseq; seq_i++)
		fprintf(detail_output, "%c", read_str[seq_i]);
	fprintf(detail_output, "\n");

}

void AssemblyHandler::print_info(bam1_t  * br ){
	int score, ori_score, chaining_score = -1;
	bam_get_num_tag(br, vcf_tag.AS, &score);
	bam_get_num_tag(br, vcf_tag.OS, &ori_score);
	bam_get_num_tag(br, vcf_tag.CS, &chaining_score);

	char* OA_tag_char = bam_get_string_tag(br, vcf_tag.OA);
	char* MV_tag_char = bam_get_string_tag(br, vcf_tag.MV);
	char* XA_tag_char = bam_get_string_tag(br, vcf_tag.XA);
	char* RC_tag_char = bam_get_string_tag(br, vcf_tag.RC);

	fprintf(detail_output, "%s " "%d %d ""%d "	"%s " "mapQ:%d " "flag: %d ",
		bam_qname(br),	br->core.tid,	br->core.pos,
		bam_is_fwd_strand(br),	bam_is_first(br)?"FIR":"SEC",	br->core.qual,	br->core.flag);
	fprintf(detail_output, "score: [%d, %d, %d]", score, ori_score, chaining_score);
	fprintf(detail_output, "[OA:%s] [MV:%s] [XA:%s] [RC:%s]", OA_tag_char, MV_tag_char, XA_tag_char, RC_tag_char);

	//char MV[4] = "MV"; //mate SV information
	//char XA[4] = "XA"; //ALT alignment
	//char RC[4] = "RC"; //Read commend at signal STEP

	//ERR894726.126014465	16	1	406062	14	126M	16	90270070	0	ACTCT	70<70	AS:i:224
	//OS:i:228	OA:Z:15,90270534,0,0,M;	CS:i:105	SV:Z:37_0_406021_402_DEL_svim	MV:Z:37_0_406021_402_DEL_svim
	//XA:Z:0,584223,0,210,F,svim;	RC:Z:15,90270534,0,228,25,25,0,0,590,RNNN0,FNNN0_

	//cigar:
	int cigar_l = br->core.n_cigar;
	uint32_t *cigar =  bam_get_cigar(br);
	for(int i = 0; i < cigar_l; i++){
		fprintf(detail_output, "%d%c", (cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[cigar[i] & BAM_CIGAR_MASK]);
	}
	fprintf(detail_output, "\t");
}

#define MAX_AVE_DEPTH 100
#define MIN_AVE_DEPTH 20
#define heterozygote_MAX_DEPTH 35

float AssemblyHandler::averageDepth(int begin_pos, int end_pos){
	int total_depth = 0;
	for(int i = begin_pos; i < end_pos; i++){
		total_depth += depth_counter[i];
	}
	return (begin_pos >= end_pos)?0:((float)total_depth / (end_pos - begin_pos));
}

void AssemblyHandler::coverage_analysis(){
	FILE* out_put_file = detail_output;
	fprintf(out_put_file, "\n ID [%d] Analysis %s ",sv_id, SV_tag_char);
	int break_point1 = o->edge_len;
	int break_point2 = sv_length - o->edge_len;
	bool pass_filter = true;
	bool overdepth = false;
	bool not_enough_depth = false;
	bool is_heterozygote = false;

	float break_point1_average_depth =  averageDepth(break_point1 - 80, break_point1);
	fprintf(out_put_file, "BP1AveDepth: %f ", break_point1_average_depth);
	if(break_point1_average_depth <= MAX_AVE_DEPTH && break_point1_average_depth >= MIN_AVE_DEPTH ){
		if(break_point1_average_depth < heterozygote_MAX_DEPTH){
			is_heterozygote = true;
		}
	}else{
		if(break_point1_average_depth > MAX_AVE_DEPTH ) overdepth = true;
		if(break_point1_average_depth < MIN_AVE_DEPTH ) not_enough_depth = true;
		pass_filter = false;
	}

	float break_point2_average_depth =  averageDepth(break_point2, break_point2 + 80);
	fprintf(out_put_file, "BP2AveDepth: %f ", break_point2_average_depth);
	if(break_point2_average_depth <= MAX_AVE_DEPTH && break_point2_average_depth >= MIN_AVE_DEPTH ){
		if(break_point2_average_depth < heterozygote_MAX_DEPTH){
			is_heterozygote = true;
		}
	}else{
		if(break_point2_average_depth > MAX_AVE_DEPTH ) overdepth = true;
		if(break_point2_average_depth < MIN_AVE_DEPTH ) not_enough_depth = true;
		pass_filter = false;
	}

	float MIN_region_depth = MIN(break_point1_average_depth, break_point2_average_depth);
	float MAX_region_depth = MAX(break_point1_average_depth, break_point2_average_depth);

	float INS_part_average_depth =  averageDepth(break_point1 + 1, break_point2 - 1);
	fprintf(out_put_file, "IPAveDepth: %f ", INS_part_average_depth);
	if(break_point1 + 10 < break_point2 - 10){ //when with INS region:
		if(INS_part_average_depth <= MAX_AVE_DEPTH && INS_part_average_depth >= MIN_AVE_DEPTH ){
			if(INS_part_average_depth < heterozygote_MAX_DEPTH){
				is_heterozygote = true;
			}
		}else{
			if(INS_part_average_depth > MAX_AVE_DEPTH ) overdepth = true;
			if(INS_part_average_depth < MIN_AVE_DEPTH ) not_enough_depth = true;
			pass_filter = false;
		}
		MIN_region_depth = MIN(MIN_region_depth, INS_part_average_depth);
		MAX_region_depth = MAX(MAX_region_depth, INS_part_average_depth);
	}

	float Normal_part1_average_depth = averageDepth(100, break_point1 - 100);
	fprintf(out_put_file, "NP1AveDepth: %f ", Normal_part1_average_depth);

	float Normal_part2_average_depth = averageDepth(break_point2 + 100, sv_length - 100);
	fprintf(out_put_file, "NP2AveDepth: %f ", Normal_part2_average_depth);

	fprintf(out_put_file, "MIN_depth: %f ", MIN_region_depth);
	fprintf(out_put_file, "MAX_depth: %f ", MAX_region_depth);

	fprintf(out_put_file, "Rst: %s %s ", pass_filter?("PASS"):("FAIL"), is_heterozygote?"hete":"homo" );
	fprintf(out_put_file, "%s %s ", overdepth?("overDp"):("normal"), not_enough_depth?"minuDp":"normal");
	fprintf(out_put_file, ">15%c >10%c ", (MIN_region_depth > 15)?('Y'):('N'), (MIN_region_depth > 10)?('Y'):('N'));
	fprintf(out_put_file, "\t");

	int print_all_read_depth = false;
	if(print_all_read_depth){
		fprintf(out_put_file, "detail of read depth");
		for(int i = 0; i < sv_length; i++){
			fprintf(out_put_file, "index: %d read depth: %d\n", i, depth_counter[i]);
		}
	}
	fprintf(out_put_file, "\n");
}

//*********************************************FUNC MIAN loader***************************************************//

int remapped_read_assembly(int argc, char *argv[]){
	fprintf(stderr, "version 1.00\n");
	SignalAssemblyOption o;
	int option_rst = o.get_option(argc, argv);
	if(option_rst != 0)
		return option_rst;

	SV_ref_sequence sf;
	sf.init(o.index_dir, o.ori_header_fn, o.ori_reference_fn, o.max_cluster_distance);

	Remapped_read_loader rl;
	rl.init(&o, sf.get_sv_info_list());

	AssemblyHandler ah;

	ah.init(&o, &sf);

	while(ah.read_register(rl, sf) != -1){
		ah.run_assembly_analysis();
	}
	rl.destory();
	ah.destroy();
	return 0;
}
