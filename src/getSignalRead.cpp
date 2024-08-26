/*
 * getSignalRead.cpp
 *
 *  Created on: 2021年3月17日
 *      Author: fenghe
 */

#include "getSignalRead.hpp"

int sam_cmp_by_name(const void*a_,const void*b_)
{
	bam1_t *a = (bam1_t *)a_;
	bam1_t *b = (bam1_t *)b_;

	int name_cmp_rst = strcmp(bam_get_qname(a), bam_get_qname(b));
	if(name_cmp_rst != 0)
		return (name_cmp_rst < 0)?-1:1;
	return (bam_is_first(a))?-1:1;
}

#define MAX_ACCECPT_READ_LEN 2048
void bam2fastqWrite_additional_str_gz(FILE * output_file, bam1_t *br, char* additionalStr1){
	//char format_String[3000];
	char 	 seq_buff[MAX_ACCECPT_READ_LEN];
	uint8_t qual_buff[MAX_ACCECPT_READ_LEN];
	const int read_len = br->core.l_qseq;
	get_bam_seq(0, read_len, seq_buff, br);//store in binary format
	get_bam_quality_str(0, read_len, qual_buff, br);
	if(!bam_is_unmapped(br) && !bam_is_fwd_strand(br)){//reverse the string and qual when bam is reversed
		getReverseStr_char(seq_buff, read_len);
		getReverseStr_qual(qual_buff, read_len);
	}
	fprintf(output_file, ""
			"@%s %s\n"
			"%s\n"
			"+\n"
			"%s\n",
			bam_qname(br), additionalStr1,
			seq_buff,
			qual_buff);
}

uint32_t READ_SIGNAL_HANDLER::getScoreByCigar(bam1_t* b){
		int score = 0;
		uint32_t* bam_cigar = bam_get_cigar(b);
		int gap_len = 0;
		for (uint i = 0; i < b->core.n_cigar; ++i)
		{
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			int penalty,penalty1,penalty2;
			switch (type)
			{
			case CIGAR_MATCH:
			case CIGAR_SEQ_MATCH:
				score += (length * match); break;
			case CIGAR_INSERT:
			case CIGAR_DELETE:
			case CIGAR_SOFT_CLIP:
			case CIGAR_HARD_CLIP:
				if(type == CIGAR_INSERT || type == CIGAR_DELETE)
					gap_len +=length;
				penalty1 = gap_open + length*gap_ex;
				penalty2 = gap_open2 + length*gap_ex2;
				penalty  = MIN(penalty1, penalty2);
				score -= penalty; break;
			case CIGAR_SKIP:
				break;
			case CIGAR_PAD:
				break;
			case CIGAR_SEQ_MISMATCH:
				break;
			default:
				break;
			}
		}
		int NM = 0;
		bam_get_num_tag(b, "NM", &NM);
		//NM = MIS + INS + DEL
		int misMatchNum = (NM-gap_len);
		score -= (mismatch + match)*misMatchNum;
		score = MAX(0, score);
		return score;
	}

//return 0 when no XA, 1~5 when has <= 5 XA , 6 when has more than 6 XA
//when mapq == 0, BWA-MEM return at most 5 alt alignment when score are same, other wise XA is blank
static int32_t get_XA_number(bam1_t* b){
	if(b->core.qual > 0)
		return 0;
	const char *XA_Str = bam_get_string_tag(b, "XA");
	if(XA_Str == NULL)
		return 6;
	int SA_count = 0;
	for(int i = 0; i < (int)strlen(XA_Str); i++)
		if(XA_Str[i] == ';')
			SA_count ++;
	return SA_count;
}

void READ_SIGNAL_HANDLER::all_signal_records_read_pair(bam1_t &read1, bam1_t &read2, bool will_Be_used){
	if(!will_Be_used) return;
	bam1_t * b[2]; b[0] = &read1; b[1] = &read2;
	bool unmapped[2];
	uint8_t mapq[2];
	int32_t isize_read[2];
	bool direction[2];
	int low_quality_segment_len[2];
	int soft_left[2]; int soft_right[2]; int clip[2];//clip:
	int INDEL_NM[2];//SNP and INDEL
	int32_t score[2];//score
	int32_t XA_number[2];
	int tid[2];//TID

	//collect informations
	for(int i = 0; i < 2; i++){
		unmapped[i] = bam_is_unmapped(b[i]);
		mapq[i] = b[i]->core.qual; //mapQ
		isize_read[i] = b[i]->core.isize;
		direction[i] = bam_is_fwd_strand(b[i]);
		low_quality_segment_len[i] = get_bam_low_quality_num(0, 100000, '/', b[i]);
		bam_has_SH_cigar(b[i], soft_left + i, soft_right + i);	clip[i] = soft_left[i] + soft_right[i];
		INDEL_NM[i] = bam_has_INDEL_NM(b[i]);
		score[i] = getScoreByCigar(b[i]);
		tid[i] = b[i]->core.tid;
		XA_number[i] = get_XA_number(b[i]);
	}

	if(discard_both_full_match && (score[0] == (b[0]->core.l_qseq *match)) && (score[1] == (b[1]->core.l_qseq *match))) return;

	//store string
	int32_t isize = ABS(isize_read[0]);
	if(isize_read[0] + isize_read[1] != 0 && isize_read[0] != isize_read[1]){ fprintf(stderr, " wrong ISIZE: %d  %d \n",isize_read[0], isize_read[1] );}
	//read pair orientation
	if(b[0]->core.pos > b[1]->core.pos){ std::swap(direction[0], direction[1]); }
	//change direction for too short insert size read pair
	if(isize == b[0]->core.l_qseq && isize == b[1]->core.l_qseq && direction[0] == false && direction[1] == true){
		std::swap(direction[0], direction[1]);
	}

	char filter_fail_reason[2][1024];
	int reason_len[2] = {0};

	//
	for(int i = 0; i < 2; i++){
		sprintf(filter_fail_reason[i] + reason_len[i], "%d_%d_%d_%d_%d_%d_%d_%d_%d_", tid[i], b[i]->core.pos, soft_left[i], score[i], mapq[i], mapq[1- i], XA_number[i], XA_number[1- i], isize);
		reason_len[i] += strlen(filter_fail_reason[i] + reason_len[i]);
	}

	char flags[2][10];
	for(int i = 0; i < 2; i++){
		int flagLen = 0;
		flags[i][flagLen++] = (bam_is_fwd_strand(b[i]) == FORWARD)?'F':'R';
		flags[i][flagLen++] = unmapped[i]?'Y':'N';
		flags[i][flagLen++] = (INDEL_NM[i] > 8)?'Y':'N';
		flags[i][flagLen++] = (clip[i] > 10)?'Y':'N';
		flags[i][flagLen++] = '0';
	}
	for(int i = 0; i < 2; i++){
		sprintf(filter_fail_reason[i] + reason_len[i], "%s_%s_", flags[i], flags[1 - i]); reason_len[i] += strlen(filter_fail_reason[i] + reason_len[i]);
	}
	//set filters
	uint32_t reasonFlag = 0;
	bool pass_filter = true;
	for(int i = 0; i < 2; i++){
		clip[i] -= low_quality_segment_len[i]; if(clip[i] < 0){low_quality_segment_len[i] = - clip[i]; clip[i] = 0;}
		low_quality_segment_len[i] >>= 1; //1 NM or INDEL in 2 low_quality_segment
		INDEL_NM[i] -= low_quality_segment_len[i]; if(INDEL_NM[i] < 0){low_quality_segment_len[i] = - INDEL_NM[i]; INDEL_NM[i] = 0;}
	}

	if(mapq[0] < 10 && mapq[1] < 10){	pass_filter = false;	reasonFlag += 1; }
	if(unmapped[0] || unmapped[1])  {	pass_filter = false;	reasonFlag += 2; }
	if(isize > 1000)				{	pass_filter = false;	reasonFlag += 4; }
	if(direction[0] != true || direction[1] != false)
									{	pass_filter = false;	reasonFlag += 8; }
	if(INDEL_NM[0] + INDEL_NM[1] > 15){	pass_filter = false;	reasonFlag += 16;}
	if(clip[0] + clip[1] > 10)		{	pass_filter = false; 	reasonFlag += 32;}
	if(tid[0] != tid[1] || tid[0] > maxTid || tid[1] > maxTid)
									{	pass_filter = false;	reasonFlag += 64;}

	if((! pass_filter) || NOT_USING_FILTER){
		//read 1 or 2 flag
		xassert(bam_is_first(b[0]), "");
		xassert(bam_is_second(b[1]), "");
		reasonFlagCounter[reasonFlag] ++;

		xassert(b[0]->core.l_qseq < MAX_ACCECPT_READ_LEN, "");

		bam2fastqWrite_additional_str_gz(output_file1, b[0], filter_fail_reason[0]);
		bam2fastqWrite_additional_str_gz(output_file2, b[1], filter_fail_reason[1]);
	}
}

#define MAX_FIAL_READ_NUM 1000
void mate_fail_handler(int &total_unpaired_read_number, bam1_t *cr, int i, int true_load_size, const char *func, const int LINE){
	total_unpaired_read_number++;
	xassert(total_unpaired_read_number < MAX_FIAL_READ_NUM, "");

	fprintf(stderr, "NUM: [%d]: Mate failed, index:[%d] in [%d] size block, tid: [ %d ], pos [%d] name [%s]\n",
			total_unpaired_read_number, i, true_load_size, cr->core.tid, cr->core.pos, cr->data );
	fprintf(stderr, "fun: [%s], line: [%d]\n", func, LINE);
}

bool READ_SIGNAL_HANDLER::signal_be_used(){
	if(sample_rate < 0.9999){
		if(rand() > sample_max_number_int)
			return false;
	}
	return true;
}


//input is a bam/cram file sort by read name
#define SAM_LOAD_BUFF_SIZE 1000000 //100k
#define SEARCH_STEP 64 //
#define SEARCH_POS_INDEX_SIZE 1000000 //1M, used for max 64*1M = 64M bp region
#define SEARCH_REGION_MAX 60000000 // 60M

//input file is SAM/BAM/CRAM sort by read start position
void READ_SIGNAL_HANDLER::SURVIVOR_SV_region_get_all_signal_records_SORT_BY_pos(char *input_bam_fn)
{
	uint64_t readNum = 0;
	//unpaired tmp files
	char tmp_fn [1024];	sprintf(tmp_fn, "%s.tmp.bam", input_bam_fn);
	htsFile *unpaired_tmp_file = hts_open(tmp_fn, "wb");//open output file
	xassert(sam_hdr_write(unpaired_tmp_file, hdr) >=0, "write header wrong!");//write header
	int total_tmp_write_num = 0;

	//get memory for all sam records
	bam1_t * sam_load_buff  	= (bam1_t *)xcalloc(SAM_LOAD_BUFF_SIZE, sizeof(bam1_t));
	uint32_t *mateIndex 		= (uint32_t *)xmalloc(SAM_LOAD_BUFF_SIZE*sizeof(uint32_t));
	uint32_t *pairSearchIndex	= (uint32_t *)xmalloc(SEARCH_POS_INDEX_SIZE*sizeof(uint32_t));

	//reset region
	//resetRegion_ID(&c_b, &r);
	//uint32_t reasonFlagCounter[1024] = {0};

	int total_unpaired_read_number = 0;

	while(1){
		//part 1: load sam records
		int true_load_size = 0;

		for(; true_load_size < SAM_LOAD_BUFF_SIZE; ){
			int sam_rst = sam_read1(c_b._hfp, hdr, sam_load_buff + true_load_size);
			if(sam_rst < 0)			break;
			if(bam_is_secondary(sam_load_buff + true_load_size))   	 continue;
			if(bam_is_supplementary(sam_load_buff + true_load_size)) continue;
			true_load_size++;
			if(sam_load_buff[true_load_size - 1].core.tid != sam_load_buff[0].core.tid) break;//break when tid changed
			if(sam_load_buff[true_load_size - 1].core.pos - sam_load_buff[0].core.pos > SEARCH_REGION_MAX) break; //break when reach the end of max search region
		}

		//part 2: finding pairs for all records
		//clear all records and building index for all records
		int pos_index_size = 0;
		int32_t buff_start_pos = sam_load_buff[0].core.pos;
		int32_t buff_final_pos = sam_load_buff[true_load_size - 2].core.pos; //not the final node
		for(int i = 0; i < true_load_size; i++){
			mateIndex[i] = MAX_uint32_t;
			int c_record_pos = sam_load_buff[i].core.pos - buff_start_pos;
			int pos_index = c_record_pos / SEARCH_STEP;
			while(pos_index >= pos_index_size){
				pairSearchIndex[pos_index_size++] = i;
				xassert(pos_index_size < SEARCH_POS_INDEX_SIZE - 1, "");
			}
		}
		pairSearchIndex[pos_index_size] = true_load_size;// store final block

		//finding pair for all records
		for(int i = 0; i < true_load_size - 1; i++){ //skip the final one
			//for all the reads:
			if(mateIndex[i] != MAX_uint32_t) continue; // skip already found pairs
			//try to find pairs
			bam1_t * c_r = sam_load_buff + i;
			if(c_r->core.tid != c_r->core.mtid) continue;
			char * read_name = bam_get_qname(c_r);
			//for reads that both unmapped
			bam1_t * true_mate = NULL;
			if(c_r->core.tid == -1){
				xassert(c_r->core.mtid == -1, "");
				if(i == 0 || i > true_load_size - 1 - 1) continue;//reads in the edge of buff
				if		(strcmp(bam_get_qname(c_r + 1), read_name) == 0)			true_mate = c_r + 1;
				else if	(strcmp(bam_get_qname(c_r - 1), read_name) == 0)			true_mate = c_r - 1;
				//store
				if(true_mate == NULL)
					mate_fail_handler(total_unpaired_read_number, c_r, i, true_load_size, __func__, __LINE__);
				else{
					int c_mate_index = true_mate - sam_load_buff;
					xassert(c_mate_index != i, "");
					mateIndex[i] = c_mate_index;
					mateIndex[c_mate_index] = i;
				}
			}
			else{
				//get mate position
				int mpos = c_r->core.mpos;
				if(mpos <= buff_start_pos || mpos >= buff_final_pos) continue;//skip mate read out of region
				int pos = c_r->core.pos;
				//search mate in list
				int c_record_pos = c_r->core.mpos - buff_start_pos;
				int pos_index = c_record_pos / SEARCH_STEP;
				int st_search_index = pairSearchIndex[pos_index];
				int ed_search_index = pairSearchIndex[pos_index + 1];
				for(int k = st_search_index; k < ed_search_index; k++){
					bam1_t * mate_try = sam_load_buff + k;
					if(mate_try->core.pos < mpos) continue;
					if(mate_try->core.pos > mpos) break;
					if(mate_try->core.mpos != pos) continue;
					if(strcmp(bam_get_qname(mate_try), read_name) == 0 && c_r != mate_try)	{//successfully found mate
						true_mate = mate_try; break;
					}
				}
				//store
				xassert(true_mate != NULL, "");
				int c_mate_index = true_mate - sam_load_buff;
				xassert(c_mate_index != i, "");
				mateIndex[i] = c_mate_index;
				mateIndex[c_mate_index] = i;
			}
		}

		//part 3: dump all reads that have no pair
		int total_unmated_read = 0;
		for(int i = 0; i < true_load_size; i++){
			if(mateIndex[i] != MAX_uint32_t) continue; // skip already found pairs
			total_unmated_read++;
			xassert(sam_write1(unpaired_tmp_file, hdr, sam_load_buff + i) >= 0, "");//write record
		}
		total_tmp_write_num += total_unmated_read;
		if(total_unmated_read * 20 > true_load_size)
			fprintf(stderr, "WARNING, too many total_unmated_read![ %d %d %f]\n",
					total_unmated_read, true_load_size, (float)total_unmated_read / true_load_size);

		//part 4: handle all paired reads
		for(int i = 0; i < true_load_size - 1; i++){ //skip the final one
			if(mateIndex[i] == MAX_uint32_t) continue; // skip already found pairs
			bam1_t * b1 = sam_load_buff + i;
			if(bam_is_second(b1)) continue; //skip second read
			bam1_t * b2 = sam_load_buff + mateIndex[i];
			readNum += 2;
			if(readNum % 100000 == 0)
				fprintf(stderr, "%ld\r", readNum);
			xassert(bam_is_first(b1), "");
			xassert(bam_is_second(b2), "");

			all_signal_records_read_pair(b1[0], b2[0], signal_be_used());
		}
		//part5: final, ending of blocks
		if(true_load_size == 0)
			break;//end of all record
	}

	//free the sam load buff
	for(int i = 0; i < SAM_LOAD_BUFF_SIZE; i++){
		if(sam_load_buff[i].data != NULL)
			free(sam_load_buff[i].data);
	}

	if(sam_load_buff != NULL){ free(sam_load_buff); 	sam_load_buff = NULL; }
	if(mateIndex != NULL)    { free(mateIndex); 		mateIndex = NULL;     }
	if(pairSearchIndex != NULL){free(pairSearchIndex); 	pairSearchIndex = NULL;}

	//close file
	hts_close(unpaired_tmp_file);

	//phase 2:
	fprintf(stderr, "phase 2:\n");
	//handle all unpaired reads in the first pass
	//load all unpaired read
	xassert(total_tmp_write_num % 2 == 0, "");

	fprintf(stderr, "Open unpaired_tmp_file:\n");
	unpaired_tmp_file = hts_open(tmp_fn, "rb");//open output file
	bam_hdr_t* unpaired_file_hdr = sam_hdr_read(unpaired_tmp_file);
	sam_load_buff = (bam1_t *)xcalloc(total_tmp_write_num, sizeof(bam1_t));
	int true_load_size = 0;

	fprintf(stderr, "Load records from unpaired_tmp_file:\n");
	for(; true_load_size < total_tmp_write_num; ){
		//fprintf(stderr, "Load records from unpaired_tmp_file: %d\n", true_load_size);
		int sam_rst = sam_read1(unpaired_tmp_file, unpaired_file_hdr, sam_load_buff + true_load_size);
		if(sam_rst < 0)			break;
		true_load_size++;
	}

	xassert(true_load_size == total_tmp_write_num, "");
	//sort all records by name
	fprintf(stderr, "Sort all unpaired records: [%d] in total\n", true_load_size);
	qsort(sam_load_buff, true_load_size, sizeof(bam1_t), sam_cmp_by_name);

	for(int i = 0; i < true_load_size; i+=2){
		readNum += 2;
		if(readNum % 100000 == 0)
			fprintf(stderr, "%ld\r", readNum);
		bam1_t * b1 = sam_load_buff + i;
		bam1_t * b2 = sam_load_buff + i + 1;
		xassert(bam_is_first(b1), "");
		xassert(bam_is_second(b2), "");

		all_signal_records_read_pair(b1[0], b2[0], signal_be_used());
	}

	//free the sam load buff
	for(int i = 0; i < SAM_LOAD_BUFF_SIZE; i++){
		if(sam_load_buff[i].data != NULL)
			free(sam_load_buff[i].data);
	}

	if(sam_load_buff != NULL){
		free(sam_load_buff); 	sam_load_buff = NULL;
	}
	hts_close(unpaired_tmp_file);
}

//input is a bam/cram file sort by read name
void READ_SIGNAL_HANDLER::SURVIVOR_SV_region_get_all_signal_records_SORT_BY_NAME()
{
	uint64_t readNum = 0;

	bam1_t b1 = {0};//BAM record for the first read in a pair
	bam1_t b2 = {0};//BAM record for the second read in a pair

	int sam_rst1 = 0; int sam_rst2 = 0;
	while (1){
		//load SAM 1 & 2
		do{	sam_rst1 = sam_read1(c_b._hfp, hdr, &b1); } while( (sam_rst1 >= 0) && (bam_is_secondary(&b1) || bam_is_supplementary(&b1)));
		do{	sam_rst2 = sam_read1(c_b._hfp, hdr, &b2);	} while( (sam_rst2 >= 0) && (bam_is_secondary(&b2) || bam_is_supplementary(&b2)));
		if(sam_rst1 < 0 || sam_rst2 < 0)		break;
		//SAME read name
		bool namesAreSame = (0 == strcmp((char *)bam_qname(&b1), (char *)bam_qname(&b2))); xassert(namesAreSame, "");
		xassert(bam_is_first(&b1), "");
		xassert(bam_is_second(&b2), "");
		readNum += 2;
		if(readNum % 100000 == 0)
			fprintf(stderr, "%ld\r", readNum);

		all_signal_records_read_pair(b1, b2, signal_be_used());
	}

	//free
	free(b1.data);
	free(b2.data);
}
