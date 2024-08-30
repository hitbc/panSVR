/*
 * SveHandler.cpp
 *
 *  Created on: 2020年4月29日
 *      Author: fenghe
 */

#include "../NovaSVgenerateVCF/SveHandler.hpp"

//******************************************filters************************************************/
#define MIN_MAPQ 20
static bool pass_mapQ_filter(bam1_t *br) {
	if (br->core.qual >= MIN_MAPQ)
		return true;
	if (bam_mate_map_qual(br) >= MIN_MAPQ && br->core.tid == br->core.mtid
			&& ABS(br->core.isize) < 20000)
		return true;
	return false;
}

//tools
inline bool is_mate_fwd(uint16_t flag) {
	return (!((flag & BAM_MATE_STRAND) != 0));
}

bool inline pos_within(int pos, int r_st, int r_ed) {
	if (pos >= r_st && pos <= r_ed)
		return true;
	return false;
}

void debug_print_flag(uint8_t flag) {
	std::cerr << "BAM_PAIRED " << (flag & 0x001) << "\t" << "BAM_PROPER_PAIR "
			<< (flag & 0x002) << "\t" << "BAM_UNMAPPED " << (flag & 0x004)
			<< "\t" << "BAM_MATE_UNMAPPED " << (flag & 0x008) << "\t"
			<< "BAM_STRAND " << (flag & 0x010) << "\t" << "BAM_MATE_STRAND "
			<< (flag & 0x020) << "\t" << "BAM_FIRST_READ " << (flag & 0x040)
			<< "\t" << "BAM_SECOND_READ " << (flag & 0x080) << "\t"
			<< "BAM_SECONDARY " << (flag & 0x100) << "\t" << "BAM_FILTER "
			<< (flag & 0x200) << "\t" << "BAM_DUPLICATE " << (flag & 0x400)
			<< "\t" << "BAM_SUPPLEMENTARY " << (flag & 0x800) << "\t" << "\t";
}

//*************************************Handle SA signals**********************************************/
void SveHandler::storeClipSignals(bool isClipAtRight, uint32_t ori_pos, uint8_t read_mapq) {
	//SVE c_sve;
	RefRegion ori_r(ref->get_chr_ID(), ori_pos, ori_pos);
	if(isClipAtRight)
		ori_r.ed_pos += 10;
	else
		ori_r.st_pos -= 10;
	//when unmapped
	//if (rst.empty()) {	//clip not mapped, insertion
	RST::T isBegin = (isClipAtRight) ? RST::BEGIN :RST::END;
	sve[SIG::SH][SV::INS].emplace_back(isBegin, true, MIN(15, read_mapq), ori_r, ori_r);
}

#define MIN_MISMATCH_QUAL 10
#define MINI_KMER_LEN 10
extern uint64_t kmerMask[33];
void SveHandler::storeMismatchSignals(bam1_t *br, READ_record &c_r) {
	//s1: get high quality NM number:(filter)
	uint8_t * tseq = ref->getRefStr(br->core.pos);
	uint8_t * qseq = read.getReadStr(c_r);
	uint8_t * qqual = read.getQualStr(c_r);

	mismatch_position.clear();

	int seq_i = 0;
	int ref_i = 0;

	uint32_t* bam_cigar = bam_get_cigar(br);
	uint32_t n_cigar = br->core.n_cigar;
	for (uint i = 0; i < n_cigar; ++i)
	{
		int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		switch (c_type){
		case CIGAR_MATCH: case CIGAR_SEQ_MATCH:
			for(int i = 0; i < c_size; i++, seq_i++, ref_i++)
				if((qseq[seq_i] != tseq[ref_i]) && (qqual[seq_i] >= MIN_MISMATCH_QUAL))
					mismatch_position.emplace_back(seq_i);
			break;
		case CIGAR_INSERT:
			for(int i = 0; i < c_size; i++, seq_i++)
				mismatch_position.emplace_back(seq_i);
			break; //do nothing
		case CIGAR_DELETE:
			for(int i = 0; i < c_size; i++, ref_i++)
				mismatch_position.emplace_back(seq_i);
			break;
		case CIGAR_SOFT_CLIP:
		case CIGAR_HARD_CLIP:
			return; //not handle reads with hard/soft clip
			break;
		default:	break;
		}
	}
	//check independence event number
	int independence_event_number = 1;
	int mismatch_position_size = mismatch_position.size();
	for(int i = 0; i < mismatch_position_size - 1 ; i++)
		if(mismatch_position[i] + 1 != mismatch_position[i + 1] && mismatch_position[i] != mismatch_position[i + 1])
			independence_event_number ++;

	if(independence_event_number < 3)
		return;

	//s2: get position of minimizer
	uint32_t kmer = bit2_nextKmer_init(qseq, MINI_KMER_LEN);
	uint32_t min_kmer  = MAX_uint32_t; int min_kmer_idx = -1;
	int kmer_number = br->core.l_qseq - MINI_KMER_LEN + 1;
	uint64_t MASK = kmerMask[MINI_KMER_LEN];
	for(int i = 0; i < kmer_number; i++){
		kmer     = bit2_nextKmerMASK( qseq + i, kmer, MINI_KMER_LEN);
		if(min_kmer > kmer){
			min_kmer = kmer;
			min_kmer_idx = i;
		}
	}
	int global_clip_point = min_kmer_idx + br->core.pos;

	//S3: decide clip point and direction
	int mis_before_minimizer = 0;
	int total_mis_size = mismatch_position.size();
	for(;mis_before_minimizer < total_mis_size && mismatch_position[mis_before_minimizer] < min_kmer_idx; mis_before_minimizer++);
	bool clip_after_minimizer = (mis_before_minimizer + mis_before_minimizer < total_mis_size);

	RefRegion ori_r(ref->get_chr_ID(), global_clip_point, global_clip_point);
	if(clip_after_minimizer)
		ori_r.ed_pos += 10;
	else
		ori_r.st_pos -= 10;
	RST::T isBegin = (clip_after_minimizer) ? RST::BEGIN :RST::END;
	sve[SIG::SH][SV::INS].emplace_back(isBegin, true, MIN(3, br->core.qual), ori_r, ori_r);
}

#define MIN_CLIP_LEN 1
void SveHandler::handleSASignal(READ_record &c_r, bam1_t *br) {
	//uint8_t clip_string[MAX_SH_SIGNAL_STR_LEN + 1];//to store clip string
	//PART1： get left clip string
	if (c_r.soft_left >= MIN_CLIP_LEN) {
		uint32_t ori_pos = c_r.position;
		storeClipSignals(false, ori_pos, br->core.qual);
	}
	//PART1： get right clip string
	if (c_r.soft_right >= MIN_CLIP_LEN) {
		uint32_t ori_pos = c_r.position + c_r.read_l - c_r.soft_right - c_r.soft_left;
		storeClipSignals(true, ori_pos, br->core.qual);
	}
	if(c_r.NM_NUM >= 3 && c_r.soft_left == 0 && c_r.soft_right == 0){
		storeMismatchSignals(br, c_r);
	}
}

//combine signals in the try list, and get the max possible break point
int SveHandler::getTopPossibilityIdx(int r_min, int r_max, SVE_L &l, bool isR1, bool is_forward,
		float &max_poss, std::vector<float> &dis_bp_percent) {
	//clear possibility
	int dis_bp_percent_size = dis_bp_percent.size();
	int p_r_size = r_max - r_min + 2;
	p_r_size = MIN(5000, p_r_size);
	for (int i = 0; i < p_r_size; i++)
		possibility_r[i] = 0;

	if(is_forward){
		for (auto idx : cmb_try_list) {	//calculation
			int region_st = ((isR1)?(l[idx].r1.st_pos):(l[idx].r2.st_pos)) - r_min;
			//fprintf(stderr, " %d %d-%d \n", idx,  ((isR1)?(l[idx].r1.st_pos):(l[idx].r2.st_pos)),  ((isR1)?(l[idx].r1.ed_pos):(l[idx].r2.ed_pos)));
			for (int i = 0; i < dis_bp_percent_size; i++)
				possibility_r[region_st + i] += (dis_bp_percent[i]);
		}
	}else{
		for (auto idx : cmb_try_list) {	//calculation
			int region_ed = ((isR1)?(l[idx].r1.ed_pos):(l[idx].r2.ed_pos)) - r_min;
			//fprintf(stderr, " %d %d-%d \n", idx,  ((isR1)?(l[idx].r1.st_pos):(l[idx].r2.st_pos)),  ((isR1)?(l[idx].r1.ed_pos):(l[idx].r2.ed_pos)));
			for (int i = 0; i < dis_bp_percent_size; i++){
				//xassert(region_ed - i >= 0, "");
				if(region_ed - i < 0) continue;
				possibility_r[region_ed - i] += (dis_bp_percent[i]);
			}
		}
	}

	float max_possibility = 0;
	int max_index = 0;
	//debug code
	for (int i = 0; i < p_r_size; i++){
		//debug code:
		if (max_possibility < possibility_r[i]) {
			max_possibility = possibility_r[i];
			max_index = i;
		}
	}
	max_poss = possibility_r[max_index];
	return (r_min + max_index);
}

//combined signals to begin/end SV type
void SveHandler::single_type_sve_combine(SVE_L &l, int min_score_cutoff, SIG::T sigT, SV::T svt)//method 2:
		{
	float MIN_ACCEPT_POSS = MIN_ACCEPT_POSSIBILITY_SH;
	int MAX_ACCEPT_REGION;
	int BP_REGION;
	std::vector<float> *breakpoint_distribution;
	if (sigT == SIG::DR) {
		MIN_ACCEPT_POSS = MIN_ACCEPT_POSSIBILITY_DR;
		MAX_ACCEPT_REGION = DR_bp_distribution.size();
		BP_REGION = 200;
		breakpoint_distribution = &(DR_bp_distribution);
	} else if (sigT == SIG::SH) {
		MIN_ACCEPT_POSS = MIN_ACCEPT_POSSIBILITY_SH;
		MAX_ACCEPT_REGION = 8;
		BP_REGION = 200;
		breakpoint_distribution = &(SH_bp_distribution);
	}
	std::sort(l.begin(), l.end(), SVE::cmp_by_position);


	auto sve_bg = l.begin();
	auto sve_ed = l.end();	//for each sve
	//for each sve
	cmb_store_tmp.clear();
	for (auto sve = sve_bg; sve < sve_ed; sve++) {
		if (sve->info.is_solid == RST::UNKNOWN)
			continue;	//already combined sve
		cmb_try_list.clear();
		RST::T is_solid = sve->info.is_solid;
		int max_score = 0;
		int r1_min = sve->r1.st_pos, r1_max = r1_min + 1;
		int r2_min = sve->r2.st_pos, r2_max = r2_min + 1;
		for (auto sve_try = sve;
				sve_try < sve_ed && sve_try->r1.region_overlap(r1_min, r1_max);
				sve_try++) {	//use sve as main, than try to combine
			if (is_solid != sve_try->info.is_solid)
				continue;	//condition 1
			if (sve_try->r2.chr_ID != sve->r2.chr_ID)
				continue;	//condition 3
			if (!sve_try->r2.region_overlap(r2_min, r2_max))
				continue;	//condition 3
			int try_idx = sve_try - sve_bg;
			cmb_try_list.emplace_back(try_idx);
			r1_min = MIN(r1_min, sve_try->r1.st_pos);
			r1_max = MAX(r1_max, sve_try->r1.ed_pos);
			r2_min = MIN(r2_min, sve_try->r2.st_pos);
			r2_max = MAX(r2_max, sve_try->r2.ed_pos);
			max_score += sve_try->info.getScore();
			sve_try->info.is_solid = RST::UNKNOWN;
		}
		//try to combine
		float max_possibility_r1, max_possibility_r2;
		//at least 3 normal clip signal or 8 multiple SNPs signals
		if (cmb_try_list.size() > 2 && max_score > sig_para.SVE_combine_min_score_step1 && (r1_max - r1_min < 5000) && (r2_max - r2_min < 5000)){
//
			bool search_forward = true;//handling R1 and is forward, for deletion, INV_1 and TRA
			if(svt == SV::INV_2) search_forward = false; //handling R1 and is reverse, for INV_2
			int r1_bp_position = getTopPossibilityIdx(r1_min, r1_max, l, true, search_forward,
					max_possibility_r1, *breakpoint_distribution);
			int max_accecpt_r1 = r1_bp_position;
			int min_accecpt_r1 = max_accecpt_r1 - MAX_ACCEPT_REGION;

			//handling R2 and reverse, for deletion
			search_forward = false;//handling R2 and reverse, for deletion, TRA and INV2
			if(svt == SV::INV_1) search_forward = true; //handling R2 and forward, for INV2
			int r2_bp_position = getTopPossibilityIdx(r2_min, r2_max, l, false, search_forward,
					max_possibility_r2, *breakpoint_distribution);
			int min_accecpt_r2 = r2_bp_position;
			int max_accecpt_r2 = min_accecpt_r2 + MAX_ACCEPT_REGION;

			if (max_possibility_r1 < MIN_ACCEPT_POSS
					&& max_possibility_r2 < MIN_ACCEPT_POSS)
				continue;
			int sve_n = 0;
			for (auto idx : cmb_try_list) {
				bool region_check_pass = false;
				if(svt != SV::INS){
					region_check_pass = (pos_within(l[idx].r1.st_pos, min_accecpt_r1, max_accecpt_r1)
							&& pos_within(l[idx].r2.ed_pos, min_accecpt_r2,	max_accecpt_r2));
				}
				else{
					region_check_pass = (is_solid == RST::BEGIN)? pos_within(l[idx].r1.st_pos, min_accecpt_r1, max_accecpt_r1) :
							pos_within(l[idx].r2.ed_pos, min_accecpt_r2, max_accecpt_r2);
				}
				if (region_check_pass)
					sve_n++;
				else
					l[idx].info.is_solid = is_solid;
			}
			//get new score:
			int score = (MAX(max_possibility_r1, max_possibility_r2))*2 / MIN_ACCEPT_POSS;
			if(score < min_score_cutoff) continue;
			//store new sve into
			cmb_store_tmp.emplace_back(sve->r1.chr_ID, r1_bp_position, sve->r2.chr_ID,
					r2_bp_position, sve_n, score, is_solid, BP_REGION);
		}
	}
	l.swap(cmb_store_tmp);
}

#define MIN_REALIGNMENT_SCORE 50
void SveHandler::handleUMSignal(bam1_t *br) {

	uint8_t *query = UMQueryBuff;
	const int read_len = br->core.l_qseq;
	if (get_bam_seq_bin(0, read_len, query, br) == false)//get string failed, return
		return;
	if (!BAM_handler::pass_compact_filter(query, read_len))// || !BAM_handler::passComplexFilter(query, read_len))//filter for string
		return;
	//store to unmapped read list
	read.storeReadUM(br, query);
}

//just equal : !bam_is_DR_signal(); but is more understandable
bool bamIsNormalDrPair(bam1_t *br, int i_size_min, int i_size_MAX) {
	bool direction = bam_is_fwd_strand(br);
	bool mate_dir = bam_is_mate_fwd_strand(br);
	int ABSisize = ABS(br->core.isize);

	if (br->core.tid == br->core.mtid && direction != mate_dir
			&& ABS(ABSisize) < i_size_MAX && ABS(ABSisize) > i_size_MAX) {
		if (direction == FORWARD && br->core.pos <= br->core.mpos)
			return true;	//normal condition 1
		if (direction == REVERSE && br->core.pos <= br->core.mpos)
			return true;	//normal condition 2
	}
	return false;
}

bool bam_aligned_analysis(bam1_t *b, int *clip_left, int *clip_right, int *gap_mismatch_inside){
	if(b->core.n_cigar == 0)//CIGAR un-available
		return false;
	int gap_number = 0;
	uint32_t* bam_cigar = bam_get_cigar(b);
	for (uint i = 0; i < b->core.n_cigar; ++i){
		int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		if(type == CIGAR_INSERT || type == CIGAR_DELETE){
			gap_number += (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		}
	}
	int NM = 0;
	bam_get_num_tag(b, "NM", &NM);
	*gap_mismatch_inside = MAX(gap_number, NM);
	//NM = MIS + INS + DEL
	int begin_type = (int)(1 + (bam_cigar[0] & BAM_CIGAR_MASK));
	int end_type   = (int)(1 + (bam_cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK));
	*clip_left = 0; *clip_right = 0;
	if(	begin_type == CIGAR_SOFT_CLIP || begin_type == CIGAR_HARD_CLIP)
		*clip_left = (bam_cigar[0] >> BAM_CIGAR_SHIFT);
	if(	end_type == CIGAR_SOFT_CLIP ||	end_type == CIGAR_HARD_CLIP)
		*clip_right = (bam_cigar[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);

	return (*clip_left > 0 || *clip_right > 0 || *gap_mismatch_inside >= 2);
}

void SveHandler::get_original_signals_from_reads(){
		Bam_file *c_b = &read.file;
		R_region region;
		ref->get_cur_region()->toR_region(region);
		resetRegion_ID(c_b, &region);	//reset region
		read_counter = 0;

		int clip_left, clip_right, gap_mismatch_inside; //the results of analysis alignment
		rdc.set_dc_st(region.st_pos);
		rdc.clear_read_depth_list();

		read.tr_loader.clear();

		while (bam_next(c_b)) {
			bam1_t *br = &(c_b->_brec);
			rdc.add_read_depth_item(read_counter, br->core.pos);			//depth counter:
			//debug code:
			read_counter++;
			if(read_counter % 10000 == 0) fprintf(stderr, "%d\n", read_counter);
			//basic filter
			if(bam_is_secondary(br))
				continue;
			if(bam_is_supplementary(br))
				continue;
			//WARNING: the unmapped reads will always not pass the MAPQ filter, therefore it should be before the mapQ filter
			if(bam_is_unmapped(br)){
				handleUMSignal(br);
				continue;
			}//unmapped read
			if (!pass_mapQ_filter(br))
				continue;
			//add new signal type:
			//SR signal
			if (bam_aligned_analysis(br, &clip_left, &clip_right, &gap_mismatch_inside)) {
				if(read.storeReadSR(br, clip_left, clip_right, gap_mismatch_inside))//store signal reads
					handleSASignal(read.read_list.back(), br);
			}
			//DR signal
			int middle_size = br->core.l_qseq - clip_left - clip_right;	//length of match bases
			if (bam_is_DR_signal(br, sig_para.insert_size_min, sig_para.insert_size_max)){
				handleDRSignal(&(br->core), middle_size);
			}

			//other reads pair not normal paired
			if(br->core.mtid != br->core.tid || ABS(br->core.pos - br->core.mpos) > 200000){
				read.tr_loader.trans_list.emplace_back(br);
			}
		}
		return;
	}

void SveHandler::handleDRSignal(bam1_core_t *core, int middle_size) {
	bool dir = (((core->flag & BAM_STRAND) == 0)), mate_dir = ((core->flag & BAM_MATE_STRAND) == 0);
	int absIsize = (int) ABS(core->isize);
	SV::T t = SV::UNKNOWN;
	RST::T is_begin = RST::UNKNOWN;
	if (core->tid == core->mtid && dir != mate_dir && absIsize < sig_para.max_del_dup_length) {//insertion/ tandem duplication/ deletion
		is_begin = (dir == FORWARD) ? RST::BEGIN : RST::END;
		bool normalOri = ((dir == FORWARD && core->pos <= core->mpos) || (dir == REVERSE && core->pos >= core->mpos));
		if (normalOri){
			if(absIsize > sig_para.insert_size_max)			t = SV::DEL;
			else if(absIsize < sig_para.insert_size_min)	t = SV::DUP;
			else 											t = SV::UNKNOWN;
		} else 												t = SV::DUP;
	} else if (core->tid == core->mtid && dir == mate_dir) {
		is_begin = (core->isize > 0) ? RST::BEGIN : RST::END;
		t = (dir == FORWARD) ? SV::INV_1 : SV::INV_2;
	} else if (core->tid != core->mtid) {
		is_begin = (core->tid < core->mtid) ? RST::BEGIN : RST::END;
		t = (dir != mate_dir) ? SV::TRA : SV::TRA_INV;
	}
	if (t != SV::UNKNOWN)
		sve[SIG::DR][t].emplace_back(core, dir, mate_dir,
				sig_para.insert_region_len, middle_size, is_begin);
}

//**********************************MAIN*********************************************************/

//combine begin and end of each signal types, to be solid or unstable; 	//when it is solid SV, is begin set to be "SV_SOLID"
void sve_begin_end_combine(SVE_L & l, int MIN_SOLID_S, int MIN_READ_NUM){
	std::sort(l.begin(), l.end(), SVE::cmp_by_position);
	int store_index = 0;
	auto sve_ed = l.end();
	for(auto sve = l.begin();sve != sve_ed; sve++){
		if(sve->info.is_solid == RST::UNKNOWN) continue;
		//search pair:
		for(auto sve_try = sve + 1; sve_try < sve_ed && (sve->r1.region_overlap(sve_try->r1)); sve_try++){//use sve as main, than try to combine
			if(sve->info.is_solid == sve_try->info.is_solid) continue;//condition 1: only combine begin+end
			if(!sve->r2.region_overlap(sve_try->r2))	continue;//condition 3
			sve->r1.Combine(sve_try->r1, true);
			sve->r2.Combine(sve_try->r2, true);
			sve->info.is_solid = RST::SOLID;
			sve->info.combine_score(sve_try->info);
			sve_try->info.is_solid = RST::UNKNOWN;
			break;
		}
		bool pass_filter = true;
		if(pass_filter && sve->info.is_solid  < RST::SOLID && (sve->info.getScore() < MIN_SOLID_S*2 || sve->info.getReadNum() < MIN_READ_NUM*2)) pass_filter = false;
		if(pass_filter && sve->info.is_solid == RST::SOLID && (sve->info.getScore() < MIN_SOLID_S || sve->info.getReadNum() < MIN_READ_NUM)) pass_filter = false;//condion 1: solid but skip TODO:: 10% of depth

		if(!pass_filter){
			std::cerr << "Not pass filter:\t" << *sve;
			continue;
		}

		l[store_index++] = sve[0];
	}
	l.resize(store_index);
}

void SveHandler::combine_duplication(SVE_L & l) {
	auto sve_bg = l.begin();
	auto sve_ed = l.end();	//for each sve
	//for each sve
	cmb_store_tmp.clear();
	for (auto sve = sve_bg; sve < sve_ed; sve++) {
		if (sve->info.is_solid == RST::UNKNOWN)
			continue;	//already combined sve
		cmb_try_list.clear();
		for (auto sve_try = sve + 1; sve_try < sve_ed && sve->r1.region_overlap(sve_try->r1); sve_try++) {	//use sve as main, than try to combine
			if(sve->SVE_combine_UNION(*sve_try)){
				sve_try->info.is_solid = RST::UNKNOWN;
			}
		}
		cmb_store_tmp.emplace_back(*sve);
	}
	l.swap(cmb_store_tmp);
}

void SveHandler::cluster_and_combine_original_signals() {
	//combine the DR signals, signals combined in two steps, firstly, signals will be clustered by overlap, then, the break points will be calculated using those signals
	std::cerr << SIG::STR[SIG::DR] << std::endl;
	for (int svt = 0; svt < SV::LEN; svt++) {
		SVE_L &svl = sve[SIG::DR][svt];
		if (svl.size() == 0)
			continue;
		//std::cerr << SV::STR[svt] << std::endl;
		single_type_sve_combine(svl, 2, SIG::DR, (SV::T)svt);//10% * read depth

		if (svt == SV::DEL) //DEL
			sve_begin_end_combine(svl, sig_para.SVE_MIN_SOLID_SCORE, sig_para.SVE_MIN_READ_NUM);
		else //DUP and BND
			sve_begin_end_combine(svl, sig_para.SVE_MIN_SOLID_SCORE*1.5, sig_para.SVE_MIN_READ_NUM*1.5);
		for(unsigned int i = 0; i < svl.size(); i++)
			svl[i].setSignalType((SV::T)svt, SIG::DR);
	}

	//combined the SH signals, signals combined in two steps, firstly, signals will be clustered by overlap, then, the break points will be calculated using those signals
	std::cerr << SIG::STR[SIG::SH] << std::endl;
	for (int svt = 0; svt < SV::LEN; svt++) {
		SVE_L &svl = sve[SIG::SH][svt];
		if (svl.size() == 0)
			continue;
		//std::cerr << SV::STR[svt] << std::endl;
		single_type_sve_combine(svl, 2, SIG::SH, (SV::T)svt);
		sve_begin_end_combine(svl, sig_para.SVE_MIN_SOLID_SCORE, sig_para.SVE_MIN_READ_NUM);
		for(unsigned int i = 0; i < svl.size(); i++)
			svl[i].setSignalType((SV::T)svt, SIG::SH);
		//delete duplication:
		combine_duplication(svl);
	}
}

void printContig(FILE * output, int ori_read_number, AssemblyContig & contig, int contig_ID, std::vector<ASS_reads_info> &ass_read_list){
	//basic part
	if(true){
		fprintf(output,
				"\ncontig_ID: [%d] "
				"word length: [%d] "
				"CONTIG size: [%ld]\n "	"CONTIG seq:  [%s]\n "
				"supportReads [%ld] "
				"ending_reason: [%d %d]"
				"new_support_read[%d] "
				,
				contig_ID,
				contig.wordLength, contig.seq.size(), contig.seq.c_str(), contig.supportReads.size(), contig.ending_reason[0], contig.ending_reason[1], contig.new_support_read);
	}
	//supplementary information
	if(true){
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
						ass_read_list[r].read_list_index);
			}
		}
	}

	fprintf(output, "\n");
}
//
#define MAX_WRONG_BASE 8
void SveHandler::get_suggention_alignment_position_list(AssemblyContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_reads_info> &ass_read_list, FILE* log_f ){

	const char *contig_seq = contig.seq.c_str();
	int contig_seq_len = contig.seq.size();
	//step1: get all the read positions at assembled CONTIG
	//set read position for each actions
	remove_read_set.clear();
	for (auto &ca : contig.actions)
		if(ca.read_ID < ori_read_number)
			ca.set_read_pos(read_list[ca.read_ID], contig.seq, contig.ass_begin_offset_in_contig, contig.wordLength, remove_read_set, ass_read_list[ca.read_ID].read_in_ref_offset);

	//step2: find the suggestion alignment position of CONTIG at reference from each read position
	suggest_st_pos_map.clear();

	//get contig coverage
	if((int)contig_depth.size() < contig_seq_len)	contig_depth.resize(contig_seq_len);
	memset(&(contig_depth[0]), 0, contig_seq_len*sizeof(uint16_t));

	for (auto &ca : contig.actions){
		ca.wrong_base = 999;
		if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
		if(ass_read_list[ca.read_ID].signal_type != Read_type::SR) continue;
		if(remove_read_set.find(ca.read_ID) != remove_read_set.end()) continue;
		//set coverage:
		const char* read_seq = read_list[ca.read_ID].c_str();
		int st_pos_ref = ca.position_in_contig - contig.ass_begin_offset_in_contig - ca.position_read;
		int ed_pos_ref = st_pos_ref + read_list[ca.read_ID].size();
		int st_pos_read = 0; if(st_pos_ref < 0){	st_pos_read -= st_pos_ref;st_pos_ref = 0;}
		ed_pos_ref = MIN(contig_seq_len, ed_pos_ref);
		ca.wrong_base = 0;
		for(int i = st_pos_ref; i < ed_pos_ref && ca.wrong_base <= MAX_WRONG_BASE; i++, st_pos_read++){
			if(read_seq[st_pos_read] != 'N' && contig_seq[i] != read_seq[st_pos_read])
				ca.wrong_base++;
		}
		if(ca.wrong_base <= MAX_WRONG_BASE){
			st_pos_read = 0;
			for(int i = st_pos_ref; i < ed_pos_ref; i++, st_pos_read++){
				if(contig_seq[i] == read_seq[st_pos_read]) contig_depth[i] ++;
			}
			std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
			if(it != suggest_st_pos_map.end())	it->second++;
			else suggest_st_pos_map[ca.suggest_contig_offset_in_ref] = 1;
		}else{
			std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
			if(it != suggest_st_pos_map.end())	it->second--;
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
	SUGGEST_pos_list.clear();

	for(std::map<int, int >::value_type &sug : suggest_st_pos_map){
		if(sug.second >= 2){
			SUGGEST_pos_list.emplace_back(sug.first);
			SUGGEST_POS_LIST_ITEM & sug = SUGGEST_pos_list.back();
			int c_suggest_pos = sug.suggest_pos;
			for(auto &ca : contig.actions){
				if(ca.read_ID < ori_read_number && ca.suggest_contig_offset_in_ref == c_suggest_pos)
					sug.add_read_start_pos(ca.position_in_contig - ca.position_read, ca.wrong_base);
			}
			sug.count_ave_read_start();
			if(sug.low_wrong_base_read_number < 1 || sug.high_wrong_base_read_number > sug.low_wrong_base_read_number * 4)//only 0 support reads or less than 20% right supported reads
				SUGGEST_pos_list.erase(SUGGEST_pos_list.end() - 1);
		}
	}

	std::sort(SUGGEST_pos_list.begin(), SUGGEST_pos_list.end(), SUGGEST_POS_LIST_ITEM::cmp_by_read_position);
	//std::sort(SUGGEST_pos_list.begin(), SUGGEST_pos_list.end(), SUGGEST_POS_LIST_ITEM::cmp_by_read_count);

	//print suggestion alignment range information
	if(false){
		//all possible ranges
		fprintf(log_f, "All suggestions:\t");
		for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++) fprintf(log_f, "[SUG: %d CN: %d]\t", it->first, it->second);
		fprintf(log_f, "\n Used suggestions: \t");
		for(SUGGEST_POS_LIST_ITEM & sug : SUGGEST_pos_list) sug.printf(log_f);
		//possible true ranges
		if(SUGGEST_pos_list.empty()) 			{fprintf(log_f, "\nNO suggestion\n"); }
		else if(SUGGEST_pos_list.size() == 1) 	{fprintf(log_f, "\nUNIQUE: [MAX_SUG: %d]\n", max_suggention);}
		else									{fprintf(log_f, "\nMULTY\n"); }

		if(true){
			for (auto &ca : contig.actions)
				if(ca.read_ID < ori_read_number)// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
				{
					ass_read_list[ca.read_ID].print(log_f, 0);
					ca.print(log_f, true);
				}
			fprintf(log_f, "\n");
		}
	}
}

extern uint8_t charToDna5n[];
void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig){
	int contig_seq_len = contig_string.size();
	const char * contig_seq = contig_string.c_str();
	xassert(nullptr != contig_seq, "");
	//store bin contig
	bin_contig.resize(contig_seq_len);
	for (int i = 0; i < contig_seq_len; ++i)
		bin_contig[i] = charToDna5n[(uint8_t)contig_seq[i]];
}

//using "bin_contig" to store contig sequence
void SveHandler::alignment_and_get_var(int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len, int ref_region_length){
	//alignment contig into reference
	int suggest_ref_ed_pos = suggest_ref_st_pos + contig_seq_len + region_addition_load;
	int suggest_contig_st_pos = 0;
	if(suggest_ref_st_pos < 0){	suggest_contig_st_pos = (- suggest_ref_st_pos);	suggest_ref_st_pos = 0;}
	suggest_ref_st_pos = MAX(0, suggest_ref_st_pos);//get 15 addition alignment region at begin
	suggest_ref_ed_pos += 60; suggest_ref_ed_pos = MIN(suggest_ref_ed_pos, ref_region_length);//get 60 additional alignment region at end
	bool contig_out_range = (suggest_ref_ed_pos < suggest_ref_st_pos + 20 || suggest_contig_st_pos > contig_seq_len); //when true region length less then 20 bp, it is out of range
	//if(output_assembly_contig) fprintf(detail_output, "[aln_ref_region %d ~ %d ] [aln contig st %d ] %s", suggest_ref_st_pos , suggest_ref_ed_pos, suggest_contig_st_pos, contig_out_range?("[ contig_out_range]"):"");
	if(!contig_out_range){
		//print alignment range
		ca.align_non_splice(&(bin_contig[0]) + suggest_contig_st_pos, contig_seq_len - suggest_contig_st_pos,
				suggest_ref_st_pos, suggest_ref_ed_pos);
		//if(output_assembly_contig) print_contig_aln_ori_cigar(ca.ez, detail_output, suggest_contig_st_pos);
		int ref_adj_size = ca.adjustCIGAR();
		suggest_ref_st_pos += ref_adj_size;
		fprintf(stderr, "[aln_ref_region %d ~ %d ] [aln contig st %d ] %s", suggest_ref_st_pos , suggest_ref_ed_pos, suggest_contig_st_pos, contig_out_range?("[ contig_out_range]"):"");
		if(true){
			ca.printf_alignment_detail(stderr, suggest_ref_st_pos, &contig_depth[0] + suggest_contig_st_pos, 0);
			fprintf(stderr, "ref_pos_enough_match_base: %d\t", 0);
		}
		ca.get_canditate_SVs(SVE_SVs, MIN_sv_len, suggest_ref_st_pos, suggest_SV_length, region_ref_global_position);
	}
}

bool SveHandler::SV_region_depth_filter(SVE & c_sve){
	float r1_max_depth, r2_max_depth;
	float r1_depth = rdc.get_ave_depth(c_sve.r1.st_pos, c_sve.r1.ed_pos, &r1_max_depth);
	float r2_depth = rdc.get_ave_depth(c_sve.r1.st_pos, c_sve.r1.ed_pos, &r2_max_depth);
	fprintf(stderr, "r1:AVE(max) depth: [%f, (%f)], r2:AVE(max) depth: [%f, (%f)] ", r1_depth, r1_max_depth, r2_depth, r2_max_depth);
	std::cerr << "r1: " <<  c_sve.r1 << ", r2:" << c_sve.r2 << std::endl;
	if(r1_max_depth > 15 || r2_max_depth > 15){
		fprintf(stderr, "Assembly will be skipped because over depth in this region\n");
		return false;
	}
	if(r1_depth > 5 || r2_depth > 5){
		fprintf(stderr, "Assembly will be skipped because over depth in this region\n");
		return false;
	}
	return true;
}

bool SveHandler::assembly_load_read(SVE &sv, RefRegion &main, RefRegion &supp){//buff to store contig
	//am->clear();//clear
	RefRegion r = *(ref->get_cur_region());
	if(sv.isInRegion(r) == false)return false;//check
	int region_st_pos = r.st_pos;

	//step1: get SV regions
	int MAX_load_read_number = 2000;
	//am->reads.reserve(MAX_load_read_number*4);//malloc for reads list
	//for one SV region, search all reads from position: [st_pos - edge_len] to [ed_pos + edge_len]
	bool with_supp = true;
	if(main.region_overlap(supp)){ main.Combine(supp, true); with_supp = false; }
	//step2: load reference for contig aligner
	//step3: load reads
	//load reads:
	ass_block_r1.clear();
	//load read signals S1: load from SR read list
	read.search_reads(Read_type::SR, ass_block_r1, main.st_pos, main.ed_pos, MAX_load_read_number, region_st_pos);
	if(with_supp) read.search_reads(Read_type::SR, ass_block_r1, supp.st_pos, supp.ed_pos, MAX_load_read_number, region_st_pos);

	int UM_st_pos = main.st_pos - 1000; UM_st_pos = MAX(UM_st_pos, 0);
	int UM_ed_pos = main.ed_pos + 1000;
	//load read signals S2: load from UM read list
	read.search_reads(Read_type::UM, ass_block_r1, UM_st_pos, UM_ed_pos, MAX_load_read_number, region_st_pos);
	//load read signals S3: load from TR read list
	read.tr_loader.load_read(main.chr_ID, UM_st_pos, UM_ed_pos,ass_block_r1, ref);

	return true;
}

//return a list:
void SveHandler::getSuggestSVlength(){
	int sum_deletions = 0; int deletion_num = 0;
	int sum_insertions = 0; int insertion_num = 0;
	suggest_SV_length.clear();
	int suggest_list_size = SUGGEST_pos_list.size();
	for(int i = 1; i < suggest_list_size; i++){
		int suggest_SV_len = SUGGEST_pos_list[i - 1].suggest_pos - SUGGEST_pos_list[i].suggest_pos;
		if(suggest_SV_len < 0){
			deletion_num ++; sum_deletions += suggest_SV_len;
		}else{
			insertion_num ++; sum_insertions += suggest_SV_len;
		}
		suggest_SV_length.emplace_back(suggest_SV_len);
	}
	if(deletion_num > 1)suggest_SV_length.emplace_back(sum_deletions);
	if(insertion_num > 1)suggest_SV_length.emplace_back(sum_insertions);
	int sum_total = sum_deletions + sum_insertions;
	int sum_number_total = deletion_num + insertion_num;
	if(sum_number_total > 1)suggest_SV_length.emplace_back(sum_total);
}

bool SveHandler::assembly_variations(SVE &sv){//buff to store contig

	int edge_len = sig_para.MaxReadLen;

	RefRegion main(sv.r1.chr_ID, sv.r1.st_pos - edge_len, sv.r1.ed_pos);
	RefRegion supp(sv.r2.chr_ID, sv.r2.st_pos - edge_len, sv.r2.ed_pos);

	main.st_pos = MAX(main.st_pos, 0);
	supp.st_pos = MAX(supp.st_pos, 0);

	//set ori_reference
	int ref_region_length;
	if(supp.ed_pos - main.st_pos < 50000)
		ref_region_length = supp.ed_pos - main.st_pos;
	else
		ref_region_length = main.getLen();
	region_ref_global_position = main.st_pos;
	ca.setRef(ref->getRefStr(region_ref_global_position), ref_region_length, ref->get_cur_region()->chr_ID, main.st_pos);
	//show reference string
	if(true){ fprintf(stderr, "C_REF_STRING: [%d-%d, len: %d]\n", main.st_pos, main.st_pos + ref_region_length, ref_region_length); for(int i = 0; i < ref_region_length; i++) fprintf(stderr, "%c", "ACGT"[ ca.getTseq(i)]); }

	if( !assembly_load_read(sv, main, supp)) return false;

	if(false){//show reads list
		fprintf(stderr, "Current read list:\n");
		int read_num = ass_block_r1.ass_read_list.size();
		int first_read_ID = (read_num > 0)?ass_block_r1.ass_read_list[0].read_list_index:0;
		for(int i = 0; i < read_num; i++){
			ass_block_r1.ass_read_list[i].print(stderr, first_read_ID);
			std::cerr << ass_block_r1.reads[i] << std::endl;
		}
	}

	ass_block_r1.run_assembly(am);
	if(ass_block_r1.contigs.empty()){  fprintf(stderr, "No assembly results\n"); return false; }

	std::vector<std::string> &read_list = ass_block_r1.reads;
	std::vector<ASS_reads_info> &ass_read_list = ass_block_r1.ass_read_list;
	int ori_read_number = ass_read_list.size();
	int contig_ID = -1;
	auto &contigs = ass_block_r1.contigs;
	FILE * detail_output = stderr;
	//handle each contig
	SVE_SVs.clear();
	for (auto &contig : contigs) {
		contig_ID++;
		if(true) printContig(detail_output, ori_read_number, contig, contig_ID, ass_read_list);
		if(contig.seq.size() < 100) continue;
		bool support_read_filter_fail = (contig_ID != 0 && (contig.new_support_read <= 2 && contig.wordLength < 100));
		if(true && support_read_filter_fail)  fprintf(detail_output, "This CONTIG is discarded, reason: 'MIN_NEW_SUPPORT_READ'\n");
		if(support_read_filter_fail) continue;

		//get suggestion alignment position
		get_suggention_alignment_position_list(contig, ori_read_number, read_list, ass_read_list, stderr);

		if(SUGGEST_pos_list.size() <= 1){
			fprintf(stderr, "No SV found in this contig, skip realignment.");
			continue;
		}

		//for each suggestion position:
		store_bin_contig(contig.seq, bin_contig);
		getSuggestSVlength();
		fprintf(stderr, "Suggest SV length list for this contig:\t");
		for(int sug_len: suggest_SV_length)
			fprintf(stderr, "%d\t", sug_len);
		fprintf(stderr, "\n");

		if(!SUGGEST_pos_list.empty())
			alignment_and_get_var(0, contig_ID, SUGGEST_pos_list[0].suggest_pos - main.st_pos, contig.seq.size(), ref_region_length);
	}

	//output candidate SVs
	uint32_t candidateSVSize = SVE_SVs.size();
	if(candidateSVSize > 0){//select best result and remove duplication results
		NOVA_SV_FINAL_RST_item::resultFilter(SVE_SVs, MIN_sv_len);
		for(uint32_t i = 0; i < candidateSVSize; i++){
			SVE_SVs[i].writeVCF(&vcfBuffStr, read.file._hdr);
			fprintf(stderr, "%s: %s", (SVE_SVs[i].will_be_output_to_vcf)?"With_VCF":"Without VCF", vcfBuffStr.s);

			if(SVE_SVs[i].will_be_output_to_vcf){
				region_SVs.emplace_back();
				std::swap(SVE_SVs[i], region_SVs.back());
			}
		}
	}
	return true;
}

bool SveHandler::assembly_variations_BND(SVE &sv, bool main_read_before_BP, bool main_read_forward, bool supp_read_before_BP, bool supp_read_forward, std::vector<NOVA_SV_FINAL_RST_item> &region_SVs, bam_hdr_t * header){//buff to store contig
	//BNDs will be treated as long deletions
	int main_tid = sv.r1.chr_ID;
	int main_pos = (sv.r1.st_pos + sv.r1.ed_pos)/2;

	int supp_tid = sv.r2.chr_ID;
	int supp_pos = (sv.r2.st_pos + sv.r2.ed_pos)/2;

	//region check
	RefRegion block_region = *(ref->get_cur_region());
	if(block_region.pos_within_same_chr_region(main_tid, main_pos)){
		//do nothing
	}else if( block_region.pos_within_same_chr_region(supp_tid, supp_pos)){
		std::swap(main_tid, supp_tid);
		std::swap(main_pos, supp_pos);
		std::swap(main_read_before_BP, supp_read_before_BP);
		std::swap(main_read_forward, supp_read_forward);
	}
	else
		return false;

	Bam_file *bam_f = &read.file;
	int normal_read_length = sig_para.MaxReadLen;
	int max_isize = sig_para.insert_size_max;

	bool afte_read_before_BP = (main_read_before_BP)?false:true;

	RefRegion withInMain(main_tid, main_pos, main_pos);			bool main_dir = main_read_forward;
	RefRegion withInSUPP(supp_tid, supp_pos, supp_pos);			bool supp_dir = supp_read_forward;
	RefRegion withInAfterMain(main_tid, main_pos, main_pos);	bool afte_dir = afte_read_before_BP;

	if(main_read_before_BP){
		withInMain.st_pos -= (max_isize - normal_read_length);
		withInAfterMain.ed_pos += (max_isize - normal_read_length - normal_read_length);
	}else{
		withInMain.ed_pos += (max_isize - normal_read_length - normal_read_length);
		withInAfterMain.st_pos -= (max_isize - normal_read_length);
	}

	if(supp_read_before_BP){
		withInSUPP.st_pos -= (max_isize - normal_read_length);
	}else{
		withInSUPP.ed_pos += (max_isize - normal_read_length - normal_read_length);
	}

	fprintf(stderr, "BND searching: withInMain\t");  withInMain.print(stderr);
	fprintf(stderr, "BND searching: withInSUPP\t");  withInSUPP.print(stderr);
	fprintf(stderr, "BND searching: withInAfterMain\t");  withInAfterMain.print(stderr);

	//only consider reads in the first region:
	R_region bam_load_region;
	bam_load_region.chr_ID = main_tid;
	bam_load_region.st_pos = main_pos - max_isize;
	bam_load_region.ed_pos = main_pos + normal_read_length + normal_read_length;

	//Analyzer
	int read_num_support_BND = 0;
	int64_t total_isize_BND = 0;
	int read_num_support_REF = 0;
	int64_t total_isize_REF = 0;

	//
	resetRegion_ID(bam_f, &bam_load_region);
	while (bam_next(bam_f)) {
		bam1_t *br = &(bam_f->_brec);
		if(bam_is_secondary(br))		continue;
		if(bam_is_supplementary(br))	continue;

		//get iSIZE
		int isize = br->core.isize;

		int tid = br->core.tid;
		int pos_read_ed = br->core.pos + br->core.l_qseq;
		int pos_read_bg = br->core.pos;
		int pos_read = (main_read_before_BP)?pos_read_ed:pos_read_bg;
		bool read_forward = bam_is_fwd_strand(br);

		int mtid = br->core.mtid;
		int mpos = br->core.mpos;
		int mpos_read = (main_read_before_BP)?mpos + br->core.l_qseq:mpos;
		bool mate_forward = bam_is_mate_fwd_strand(br);

		int suppor_int = 0;
		if(	withInMain.pos_within_same_chr_region(tid, pos_read) && withInSUPP.pos_within_same_chr_region(mtid, mpos_read) && main_dir == read_forward && supp_dir == mate_forward){
			read_num_support_BND++;
			total_isize_BND += isize;
			suppor_int += 1;
		}

		mpos_read = (afte_read_before_BP)?mpos + br->core.l_qseq:mpos;
		if(withInMain.pos_within_same_chr_region(tid, pos_read) && withInAfterMain.pos_within_same_chr_region(mtid, mpos_read) && main_dir == read_forward && afte_dir == mate_forward){//for read pairs that support the reference
			read_num_support_REF++;
			total_isize_REF += isize;
			suppor_int += 2;
		}

		fprintf(stderr, "[tid %d, pos_read_bg %d pos_read_end %d dir: %d] , [mtid %d, mpos_bg %d mpos_ed %d dir %d ] suppor_int %d\n", tid, pos_read_bg, pos_read_ed , read_forward, mtid, mpos, mpos + br->core.l_qseq, mate_forward, suppor_int);
	}

	fprintf(stderr, "Read_num_support_deletion %d %f\n",
			read_num_support_BND, (float)total_isize_BND/read_num_support_BND);
	fprintf(stderr, "Read_num_support_reference %d %f\n",
			read_num_support_REF, (float)total_isize_REF/read_num_support_REF);

	int genotype = 0;
	if((read_num_support_BND >= read_num_support_REF * 4))			genotype = 2;
	else if((read_num_support_BND * 1.5 >= read_num_support_REF))	genotype = 1;
	else														 	genotype = 0;

	if(read_num_support_BND < 15) //at least 15 read support
		genotype = 0;

	fprintf(stderr, "BND genotype: %d\n", genotype);

	if(genotype > 0){
		char BND_str[1024];
		for(int mode = 0; mode < 2; mode ++){
			if(mode == 1){
				if( block_region.pos_within_same_chr_region(supp_tid, supp_pos) && sv.info.readNumB >= 15 && sv.info.readNumE >= 15){
					std::swap(main_tid, supp_tid);
					std::swap(main_pos, supp_pos);
					std::swap(main_read_before_BP, supp_read_before_BP);
					std::swap(main_read_forward, supp_read_forward);
				}
				else
					continue;
			}
			if(sv.info.type_ID == SV::TRA){
				if(main_read_forward)
					sprintf(BND_str, "N[%s:%d[", header->target_name[supp_tid], supp_pos);
				else
					sprintf(BND_str, "]%s:%d]N", header->target_name[supp_tid], supp_pos);
			}
			else if(sv.info.type_ID == SV::INV_1)
				sprintf(BND_str, "N]%s:%d]", header->target_name[supp_tid], supp_pos);
			else if(sv.info.type_ID == SV::INV_2){
				sprintf(BND_str, "[%s:%d[N", header->target_name[supp_tid], supp_pos);
			}else
				continue;
			NOVA_SV_FINAL_RST_item::store_BND(region_SVs, main_tid, main_pos, "N", BND_str, SV::STR[sv.info.type_ID].c_str(), read_num_support_BND, read_num_support_REF, 0, genotype);
		}
	}

	return true;
}
