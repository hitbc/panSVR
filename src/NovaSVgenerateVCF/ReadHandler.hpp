/*
 * BamHandler.hpp
 *
 *  Created on: 2020年4月28日
 *      Author: fenghe
 */

#ifndef BAMHANDLER_HPP_
#define BAMHANDLER_HPP_

extern "C"
{
	#include "../clib/bam_file.h"
	#include "../clib/vcf_file.h"
extern uint64_t kmerMask[33];
}

#include"../NovaSVgenerateVCF/RefHandler.hpp"
#include "../cpp_lib/Assembler/assembler.hpp"
#include"../NovaSVgenerateVCF/sve.hpp"
#include<vector>
#include<algorithm>

struct READ_record{
	READ_record(uint16_t flag_, uint8_t direction_, uint64_t position_, int read_l_, int cigar_l_,
			int read_index_, int quality_index_, int cigar_index_,	int soft_left_, int soft_right_, int NM_NUM_) :
			flag(flag_), direction(direction_), position(position_), read_l(read_l_), cigar_l(cigar_l_),
			read_index(read_index_), quality_index(quality_index_), cigar_index(cigar_index_),
			soft_left(soft_left_), soft_right(soft_right_), NM_NUM(NM_NUM_) {
	}
	//basic
	uint16_t 	flag;
	uint8_t 	direction; // the direction of the read, and will be direction of mate when reads are unmapped
	int		 	position; // the position of the read, and will be position of mate when reads are unmapped
	int 		read_l;
	int 		cigar_l;
	int 	 	read_index;
	int 		quality_index;
	int  		cigar_index;
	int soft_left;
	int soft_right;
	int NM_NUM;
};

namespace Read_type{
	enum T {SR, DR, UM, TL, unknown};
}

struct ASS_reads_info{
	ASS_reads_info(int read_in_ref_offset_, int read_st_pos_, int read_ed_pos_, int read_list_index_, int soft_left_, int soft_right_, int number_NM_, Read_type::T signal_type_){
		signal_type = signal_type_;
		read_in_ref_offset = read_in_ref_offset_;
		read_list_index = read_list_index_;
		read_st_pos = read_st_pos_;
		read_ed_pos = read_ed_pos_;
		soft_left = soft_left_;
		soft_right = soft_right_;
		number_NM = number_NM_;
	}
	int read_list_index;
	int read_in_ref_offset;
	int read_st_pos;//read_in_ref_offset - left_clip
	int read_ed_pos;//read_in_ref_offset + middle_size
	int soft_left;
	int soft_right;
	int number_NM;
	Read_type::T signal_type;
	void print(FILE* log_file, int first_read_ID){
		fprintf(log_file, "[type: ");
		switch (signal_type) {
			case Read_type::SR: fprintf(log_file, "SR"); break;
			case Read_type::DR: fprintf(log_file, "DR"); break;
			case Read_type::UM: fprintf(log_file, "UM"); break;
			case Read_type::TL: fprintf(log_file, "TL"); break;
			default:  fprintf(log_file, "unknown"); break;
		}
		fprintf(log_file, " id: %d(%d) {%d~%d}, ref pos %d, {%d, %d, %d}]\t",
				read_list_index - first_read_ID, read_list_index, read_st_pos, read_ed_pos, read_in_ref_offset, soft_left, soft_right, number_NM);
	}
};

struct Ass_Block{
	//input read info
	std::vector<ASS_reads_info> ass_read_list;
	//input read data
	std::vector<std::string> reads;
	//output: read result
	std::vector<AssemblyContig> contigs;

	void clear(){
		ass_read_list.clear();
		reads.clear();
	}

#define MIN_BASE_QUAL 6
	void add_read_bin(Read_type::T t, uint8_t * read_str, uint8_t * qual_str, int read_len, int read_in_ref_offset_, int left_clip, int right_clip, int NM_num, int read_list_index_){
		int read_st_pos = read_in_ref_offset_ + left_clip;
		int read_ed_pos =  read_in_ref_offset_ + read_len - right_clip;
		ass_read_list.emplace_back(read_in_ref_offset_, read_st_pos, read_ed_pos, read_list_index_, left_clip, right_clip, NM_num, t);
		reads.emplace_back(); std::string &s_store = reads.back(); //string to store: string pointer type
		s_store.resize(read_len);
		for(int i = 0; i < read_len; i++){
			s_store[i] = ("ACGTN"[read_str[i]]);//store binary to acgt
			if(qual_str[i] < MIN_BASE_QUAL)
				s_store[i] = 'N';
		}
	}

	//adding trans located reads(or trans-located reads) into assembly read list, compared with normal reads, those data has no position or cigar information
	//besides, those reads are stored in CHAR
	void add_read_TL(char * read_str, uint8_t * qual_str, int read_len, bool store_in_reverse, int read_list_index_){
		for(int i = 0; i < read_len; i++)
			if(qual_str[i] < MIN_BASE_QUAL)
				read_str[i] = 'N';
		if(store_in_reverse)
			getReverseStr_char(read_str, read_len);

		ass_read_list.emplace_back(-1, -1, -1, read_list_index_, -1, -1, -1, Read_type::TL);
		reads.emplace_back(read_str);
	}

	void run_assembly(MainAssemblyHandler *am){
		std::swap(am->reads, reads);
		am->assembley();
		std::swap(am->contigs, contigs);
		std::swap(am->reads, reads);//swap back the read list
	}
};

//used to store trans-location read pairs:
struct trans_Read_item{
	int32_t tid;
	int32_t mtid;
	int32_t pos;
	int32_t mpos;
	uint8_t flag;

	trans_Read_item(	int32_t tid_, int32_t mtid_, int32_t pos_, int32_t mpos_, uint8_t flag_){ tid = tid_, mtid = mtid_, pos = pos_, mpos = mpos_, flag = flag_; }
	trans_Read_item(	bam1_t *br){ tid = br->core.tid, mtid = br->core.mtid, pos = br->core.pos, mpos = br->core.mpos, flag = br->core.flag; }
	//copy:
	trans_Read_item(const trans_Read_item &b){ memcpy(this, &b, sizeof(trans_Read_item)); }
	//sort mtid:
	static inline int cmp_by_mtid(const trans_Read_item &a, const trans_Read_item &b){
		if(a.mtid == b.mtid)	return a.mpos < b.mpos;
		else							return a.mtid < b.mtid;
	}

	void prinf(FILE * log){
		fprintf(stderr, "[cur:%d:%d ; mate: %d:%d]\n", tid, pos, mtid, mpos );
	}
};

struct TransReadLoader{

	void search_read_bg_idx(int st_pos){
		int begin_idx = last_read_end_index;
		cur_bg_idx =  MAX_int32t;
		if(trans_list.empty()) return;
		int trans_list_size =  (int)trans_list.size();
		if(begin_idx >= trans_list_size){ begin_idx = trans_list_size - 1; }
		bool search_forward = (trans_list[begin_idx].pos < st_pos);
		if(search_forward){
			for(;begin_idx < trans_list_size && trans_list[begin_idx].pos < st_pos; begin_idx++);
		}else
			for(;begin_idx >= 0 && trans_list[begin_idx].pos > st_pos; begin_idx--);
		cur_bg_idx = begin_idx;
	}

	void search_read_ed_idx(int ed_pos){
		int end_idx = cur_bg_idx;
		int trans_list_size =  trans_list.size();
		for(;end_idx < trans_list_size && trans_list[end_idx].pos < ed_pos; end_idx++);
		cur_ed_idx = end_idx;
	}

	void load_read(int tid, int st_pos, int ed_pos, Ass_Block & ab, RefHandler *ref){
		std::vector<std::string> read_list;
		//search reads in reads list:
		search_read_bg_idx(st_pos);
		search_read_ed_idx(ed_pos);
		//copy and sort
		region_trans_list.clear();
		for(int i = cur_bg_idx; i < cur_ed_idx; i++)
			region_trans_list.emplace_back(trans_list[i]);
		std::sort(region_trans_list.begin(), region_trans_list.end(), trans_Read_item::cmp_by_mtid);
		//cluster:
		int region_size = region_trans_list.size();
		for(int bg_idx = 0; bg_idx < region_size;){
			trans_Read_item &ctr = region_trans_list[bg_idx];
			int try_idx = bg_idx;
			for(;try_idx < region_size; try_idx++){
				if(ctr.mtid != region_trans_list[try_idx].mtid || ctr.mpos +3000 <  region_trans_list[try_idx].mpos)
					break;
			}
			int cluster_read_number = try_idx - bg_idx;

			int UM_read_id = 0;
			if(cluster_read_number >= 3){
				fprintf(stderr, "[ cluster_read_number:%d]\n", cluster_read_number);
				R_region region;
				region.chr_ID = ctr.mtid;
				region.st_pos = ctr.mpos + 1;
				region.ed_pos = region_trans_list[try_idx - 1].mpos + 1;
				int load_ref_length = 0;
				char * c_reference = ref->load_ref_by_region(region.chr_ID, region.st_pos, region.ed_pos + 100, &load_ref_length);

				std::map<char, int> char_count; bool low_complex = false;
				for (int i = 0; i < load_ref_length; i++)
					char_count[c_reference[i]]++;
				for(std::map<char, int>::value_type & c: char_count){
					if(c.second > (load_ref_length * 0.8 + 1))
					low_complex = true;
				}
				free(c_reference);

				if(!low_complex){
					resetRegion_ID(bam, &region);	//reset region
					//reference check:
					while (bam_next(bam)) {
						bam1_t *br = &(bam->_brec);
						if(bam_is_secondary(br))		continue;
						if(bam_is_supplementary(br))	continue;
						if(br->core.mtid == tid && br->core.mpos > st_pos && br->core.mpos < ed_pos){
							//load reads
							int read_len = br->core.l_qseq;
							get_bam_seq(0, read_len, storeReadBuff, br);
							//if dir == mdir: reverse the string:
							bool b_forward = bam_is_fwd_strand(br);
							bool m_forward = bam_is_mate_fwd_strand(br);
							bool store_in_reverse = (b_forward == m_forward);
							ab.add_read_TL(storeReadBuff, bam_get_qual(br),read_len, store_in_reverse, UM_read_id++);
							if(false){
								//debug Code:
								char * bam_name = (char *)bam_qname(br);//string values(;, ref)
								fprintf(stderr, "[TL_LOADING: :%s %d:%d ; mate: %d:%d dir:(%c, %c) %s]\n", bam_name, br->core.tid, br->core.pos, br->core.mtid, br->core.mpos, b_forward?'F':'R', m_forward?'F':'R', storeReadBuff);
							}
						}
					}
				}
			}
			bg_idx += cluster_read_number;
		}
	}

	void init(Bam_file *bam_){
		storeReadBuff = (char *)xcalloc(2000, 1);
		bam = bam_;
	}

	void clear(){
		last_read_end_index = 0;
		trans_list.clear();
	}
	int cur_bg_idx;
	int cur_ed_idx;
	int last_read_end_index;
	std::vector<trans_Read_item> trans_list;
	std::vector<trans_Read_item> region_trans_list;
	//pointer
	Bam_file *bam;

	char * storeReadBuff;
};


struct BAM_handler{
public:
	typedef std::vector<READ_record> READ_LIST;
	void init(char * bamFileName, char * ref_file_name){
		const int MAX_READ_INDEX_SIZE = SEGMENT_LEN/100 + 1;
		sr_read_index = new int [MAX_READ_INDEX_SIZE];
		um_read_index = new int [MAX_READ_INDEX_SIZE];
		storeReadBuff = new uint8_t[5000];
		kv_init(base);
		kv_init(quality);
		kv_init(cigar);
		bam_file_open(bamFileName, ref_file_name, NULL, &file);
		tr_loader.init(&file);
	}

	void destroy(){
		delete[]sr_read_index;
		delete[]um_read_index;
		delete[]storeReadBuff;
		bam_file_close(&file);
	}
	bool clip_low_quality_Filter(bam1_t *br, int read_len, int soft_left, int soft_right);
	bool storeReadSR(bam1_t *br, int soft_left, int soft_right, int gap_mismatch_inside);// store normal SR signals
	void storeReadUM(bam1_t *br, uint8_t *query);
	void storeReadCore(	Read_type::T t, int readLen, uint8_t *seq, uint8_t * qual,	int cigarLen,
			uint32_t* bam_cigar, uint16_t flag, int32_t pos,	int soft_left, int soft_right, int gap_mismatch_inside);

	void getReadStr(uint8_t * to, READ_record& c_r, int bg, int len){
		memcpy(to, base.a + c_r.read_index + bg, len);
	}
	inline uint8_t* getReadStr(READ_record& c_r){ return base.a + c_r.read_index;}
	inline uint8_t* getQualStr(READ_record& c_r){ return quality.a + c_r.quality_index;}

	bool isSameHeader(const BAM_handler &B) const;
	static bool pass_compact_filter(uint8_t *s, const int len);
	//bool clip_AAA_TailFilter(uint8_t *query, const int read_len);//when a clip string clip to an "AAAAAAAA..." tail, return true
	bool clip_AAA_TailFilter(uint8_t *read_bin, const int read_len, int soft_left, int soft_right);

	void clear(){
		read_list.clear();
		um_read_list.clear();
		base.n = 0;
		quality.n = 0;
		cigar.n = 0;
		index_already_built = false;
	}

	void build_read_index_SINGLE(READ_LIST &l, int *index, int region_st_index){
		for(int i = 0; i < READ_INDEX_SIZE; i++)
			index[i] = -1;
		int read_ID = 0;
		for(auto r: l){
			int c_index = (r.position - region_st_index)/ 100;
			if(c_index >= 0 && c_index < READ_INDEX_SIZE && index[c_index] == -1)
				index[c_index] = read_ID;
			read_ID++;
		}
	}

	inline void build_read_index(int region_st_index){
		if(index_already_built) return;
		build_read_index_SINGLE(read_list, sr_read_index, region_st_index);
		build_read_index_SINGLE(um_read_list, um_read_index, region_st_index);
		index_already_built = true;
	}

	int search_reads(Read_type::T t, Ass_Block & ab, int st_pos, int ed_pos, int max_load, int region_st_index){
		//init:
		if(t != Read_type::SR && t != Read_type::UM) return 0;
		if(index_already_built == false) build_read_index(region_st_index);
		if(st_pos > ed_pos) return 0;

		//select read and index list
		READ_LIST *c_read_l = NULL; int *index = NULL;
		if(t == Read_type::SR){ c_read_l = &read_list; index = sr_read_index; }
		else				  {	c_read_l = &um_read_list; index = um_read_index; }
		READ_LIST &l = *c_read_l;

		//loading read signals
		int c_index = (st_pos - region_st_index) / 100;
		if(c_index < 0) c_index = 0;
		if(c_index >= READ_INDEX_SIZE) c_index = READ_INDEX_SIZE - 1;

		while(index[c_index] == -1)
			c_index++;
		int search_st = index[c_index];
		int list_size = l.size();
		if(list_size == 0) return 0;
		if(search_st >= list_size) return 0;
		//search the beginning
		for(;search_st < list_size; search_st ++)
			if(l[search_st].position >= st_pos)
				break;

		//search for the end
		int total_load_num = 0;
		for(auto r = l.begin() + search_st; r->position <= ed_pos && search_st < list_size; r++, search_st ++){
			ab.add_read_bin(t, getReadStr(*r), getQualStr(*r), r->read_l, r->position - r->soft_left,  r->soft_left,  r->soft_right ,r->NM_NUM,r - l.begin());
			if(++total_load_num >= max_load) break;
		}
		return total_load_num;
	}

	READ_LIST read_list;
	READ_LIST um_read_list;

	TransReadLoader tr_loader;
	//std::vector<trans_Read_item> trans_read_list;//used to store trans-location read pairs:

	Bam_file file;
private:
	BAM_handler(const BAM_handler &b);//can`t be copy

	//read lists and index
	bool index_already_built = false;
	const static int READ_INDEX_SIZE = SEGMENT_LEN/100 + 1;
	int * sr_read_index;
	int * um_read_index;
	uint8_t * storeReadBuff;

	//buffs
	kvec_T(uint8_t, BASE_STR);//store read string
	kvec_T(uint8_t, QUALITY_STR);

	BASE_STR base;
	QUALITY_STR quality;
	path_t cigar;
};

#endif /* BAMHANDLER_HPP_ */
