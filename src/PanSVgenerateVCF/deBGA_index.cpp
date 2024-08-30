#include "deBGA_index.hpp"

#define POS_N_MAX 500
#define POS_N_MAX_LEVEL2 8000
#define RANDOM_NUM 500

#define KMER_LEN_FIRST_LEVEL 14
#define UPPER_PART_SHIFT ((LEN_KMER - KMER_LEN_FIRST_LEVEL) << 1) //（20 - 14）X 2 = 12
#define LOWER_PART_MASK 0Xfff //k = 12

//load from file, return the true load data size(in byte)
uint64_t load_data_from_file(const char * path_name, const char * fn, void ** data, bool useCalloc, int additional_byte){
	char full_fn[ROUTE_LENGTH_MAX] = {0};
	strcpy(full_fn, path_name);
	if(full_fn[strlen(full_fn) - 1] != '/')	 strcat(full_fn, "/");
	strcat(full_fn, fn);
	FILE *fp_us_b = xopen (full_fn, "rb" );

	fseek(fp_us_b, 0, SEEK_END);// non-portable
	uint64_t file_size = ftell(fp_us_b);
	rewind(fp_us_b);

	if(useCalloc)
		(*data) = (uint64_t* ) xcalloc (file_size + additional_byte, 1);
	else
		(*data) = (uint64_t* ) xmalloc (file_size + additional_byte);

	xread ((*data), 1, file_size, fp_us_b);
	fclose(fp_us_b);
	return file_size;
}

int deBGA_INDEX::load_index_file(char *index_dir)
{
	fprintf(stderr, "Begin loading index @%s\n", index_dir);
    //        buffer_ref_seq:store reference seq;load from dm file ref.seq
    result_ref_seq = (load_data_from_file(index_dir, "ref.seq", ((void ** )&buffer_ref_seq), true, 536) >> 3);
    //        buffer_seq:store unipath seq(concatenate all unipaths's seq);load from dm file unipath.seqb
    result_seq = (load_data_from_file(index_dir, "unipath.seqb", ((void ** )&buffer_seq), false, 0) >> 3);
//        buffer_seqf:store unipath's offset on unipath seq;load from dm file unipath.seqfb
    result_seqf = (load_data_from_file(index_dir, "unipath.seqfb", ((void ** )&buffer_seqf), false, 0) >> 3);
    //buffer_p: store every unipath's set of positions on reference;load from dm file unipath.pos
    result_p = (load_data_from_file(index_dir,  "unipath.pos", ((void ** )&buffer_p), false, 0) >> 3);
    //buffer_pp: unipath's pointer to array buffer_p;load from dm file unipath.posp
    result_pp = (load_data_from_file(index_dir,   "unipath.posp", ((void ** )&buffer_pp), false, 0) >> 3);
    //buffer_hash_g: store kmer's hash part; load from dm file unipath_g.hash
    result_hash_g = (load_data_from_file(index_dir,  "unipath_g.hash", ((void ** )&buffer_hash_g), false, 0) >> 3);
// buffer_kmer_g:store kmer's kmer part;load from dm file unipath_g.kmer
    result_kmer_g = (load_data_from_file(index_dir,  "unipath_g.kmer", ((void ** )&buffer_kmer_g), false, 0) >> 2);
//buffer_off_g: store kmer's offset on unipath seq; load from dm file unipath_g.offset
    result_off_g = (load_data_from_file(index_dir,  "unipath_g.offset", ((void ** )&buffer_off_g), false, 0) >> 2);

    //********************************************************************************************
    //read chr names and length from unipath.chr
    char unichr[ROUTE_LENGTH_MAX] = {0};
    strcpy(unichr, index_dir);
    strcat(unichr, "unipath.chr");

    FILE *fp_chr = xopen (unichr, "r" );
    uint32_t chr_line_n = 0;
    fscanf(fp_chr,"%s",chr_line_content);
    while(!feof(fp_chr))
    {
        if ((chr_line_n & 0X1) == 0)	strcpy(chr_names[chr_file_n],chr_line_content);
        else				            sscanf(chr_line_content, "%u", &chr_end_n[chr_file_n++]);
        fflush(stdout);
        chr_line_n++;
        fscanf(fp_chr,"%s",chr_line_content);
    }
    chr_end_n[0] = START_POS_REF + 1; //START_POS_REF = 0 record the start position of reference
    strcpy(chr_names[chr_file_n], "*");
    reference_len = chr_end_n[chr_file_n - 1];
    fclose(fp_chr);

    //building search index for CHR
    building_chr_index();
    building_bam_header();
    fprintf(stderr, "End loading index\n");
    return 0;
}

//search kmer in the index, the search result "index of kmer" stored in range, range[0] is the start index, and range[1] is the end index
//return false when search failed, otherwise return 0
bool deBGA_INDEX::search_kmer(int len_k, uint64_t kmer, int64_t range[2], int8_t seed_offset){
	if (len_k == KMER_LEN_FIRST_LEVEL){ //k = 14
		range[0] = buffer_hash_g[kmer];
		range[1] = buffer_hash_g[kmer + 1] - 1;
		if (range[1] < range[0])
			return false;
	}
	else{
		uint64_t seed_kmer = (kmer & LOWER_PART_MASK);
		uint64_t seed_hash = (kmer >> UPPER_PART_SHIFT);

		int result = binsearch_range(seed_kmer, buffer_kmer_g + buffer_hash_g[seed_hash], buffer_hash_g[seed_hash + 1] - buffer_hash_g[seed_hash], range, seed_offset<<1);
		if (result == -1)	return false;
		range[0] += buffer_hash_g[seed_hash];
		range[1] += buffer_hash_g[seed_hash];
	}
	return true;
}

//function: given an index of kmer, search the MEM within a UNITIG
//return 0 when successfully stored data, -1 when failed
int deBGA_INDEX::UNITIG_MEM_search(uint64_t kmer_index, std::vector<vertex_MEM>& vertexm_v, uint64_t *read_bit, uint32_t read_off, uint32_t read_length, int len_k, uint32_t & max_right_i){

	uint64_t kmer_pos_uni = buffer_off_g[kmer_index];//this kmer's offset on unipath seq
	//find the UID of this kmer
	int64_t seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
	uint64_t ref_pos_n = buffer_pp[seed_id_r + 1] - buffer_pp[seed_id_r];

	uint32_t uni_offset_s_l = kmer_pos_uni - buffer_seqf[seed_id_r];
	uint32_t uni_offset_s_r = buffer_seqf[seed_id_r + 1] - (kmer_pos_uni + len_k);
	//extend the kmer to a exact match
	uint32_t left_i, right_i;
	for(left_i = 1; (left_i <= uni_offset_s_l) && (left_i <= read_off); left_i++)
	{
		if(((buffer_seq[(kmer_pos_uni - left_i) >> 5] >> ((31 - ((kmer_pos_uni - left_i) & 0X1f)) << 1)) & 0X3)
				!= ((read_bit[(read_off - left_i) >> 5] >> ((31 - ((read_off - left_i) & 0X1f)) << 1)) & 0X3)
				) break;
	}

	for(right_i = 1; (right_i <= uni_offset_s_r) && (right_i <= read_length - read_off - len_k); right_i++)
	{
		if(((buffer_seq[(kmer_pos_uni + len_k - 1 + right_i) >> 5] >> ((31 - ((kmer_pos_uni + len_k - 1 + right_i) & 0X1f)) << 1)) & 0X3)
				!= ((read_bit[(read_off + len_k - 1 + right_i) >> 5] >> ((31 - ((read_off + len_k - 1 + right_i) & 0X1f)) << 1)) & 0X3)
				) break;
	}

	uint32_t read_pos = read_off + 1 - left_i;
	uint32_t mem_length = len_k + left_i + right_i - 2;

	vertexm_v.emplace_back();
	vertex_MEM & vertexm = vertexm_v.back();
	vertexm.uid = seed_id_r;
	vertexm.seed_id = vertexm_v.size() - 1;
	vertexm.read_pos = read_pos;
	vertexm.uni_pos_off = uni_offset_s_l + 1 - left_i;
	vertexm.length = mem_length;
	vertexm.pos_n = ref_pos_n;

	if (right_i > max_right_i)
		max_right_i = right_i;

	return 0;
}

#define waitingLen 3 //at most 2 bases in gap, disable multy-SNP
#define	Eindel 1// disable INDLE in SEED merging
//function: merge all MEMs within a unipath, the result stored in vertexu_v
void deBGA_INDEX::merge_seed_in_unipath(std::vector<vertex_MEM>& vertexm_v, std::vector<vertex_U>& vertexu_v){
	uint32_t mem_i = vertexm_v.size();
	vertexu_v.clear();
	//merge seeds in the same unipath
	uint32_t s1, e1;
	if(mem_i == 0){} // do nothing
	else if (mem_i == 1)
	{
		vertexu_v.emplace_back();
		vertex_U & vertexu = vertexu_v.back();
		vertex_MEM & vertexm = vertexm_v.back();
		vertexu.uid = vertexm.uid;
		vertexu.read_pos = vertexm.read_pos;
		vertexu.uni_pos_off = vertexm.uni_pos_off; //whether to set uni_pos_off to the leftest position
		vertexu.pos_n = vertexm.pos_n;
		vertexu.length1 = vertexm.length;
		vertexu.length2 = vertexm.length;
		vertexu.cov = vertexm.length;
	}else{
		qsort(&(vertexm_v[0]), mem_i, sizeof(vertex_MEM), vertex_MEM::cmp);
		uint64_t uni_id_temp = vertexm_v[0].uid;
		uint32_t j = 0;
		uint32_t cov = 0;
		while (j < mem_i)
		{
			s1 = j;
			cov = vertexm_v[s1].length;
			j++;
			while ((uni_id_temp == vertexm_v[j].uid) && (vertexm_v[j].uni_pos_off > vertexm_v[j-1].uni_pos_off) && (j < mem_i))
			{
				int diff = (int)(vertexm_v[j].read_pos - vertexm_v[j-1].read_pos - vertexm_v[j-1].length);
				if (diff > waitingLen)
					break;
				int c_eindel = (vertexm_v[j].uni_pos_off - vertexm_v[j-1].uni_pos_off) - (vertexm_v[j].read_pos - vertexm_v[j-1].read_pos);
				if (std::abs(c_eindel) < Eindel)
				{
					cov += (diff > 0)? vertexm_v[j].length : (diff + vertexm_v[j].length);
					++j;
				}
				else
					break;
			}
			e1 = j - 1;
			//store vertexu_v
			vertexu_v.emplace_back();
			vertex_U & vertexu = vertexu_v.back();
			vertexu.uid = vertexm_v[s1].uid;
			vertexu.read_pos = vertexm_v[s1].read_pos;
			vertexu.uni_pos_off = vertexm_v[s1].uni_pos_off; //whether to set uni_pos_off to the leftest position
			vertexu.pos_n = vertexm_v[s1].pos_n;
			vertexu.cov = cov;
			cov = 0;

			//cal length
			if (s1 == e1)
			{
				vertexu.length1 = vertexm_v[s1].length;
				vertexu.length2 = vertexm_v[s1].length;
			}else{
				vertexu.length1 = vertexm_v[e1].read_pos + vertexm_v[e1].length - vertexm_v[s1].read_pos;
				vertexu.length2 = vertexm_v[e1].uni_pos_off + vertexm_v[e1].length - vertexm_v[s1].uni_pos_off;
			}

			uni_id_temp = vertexm_v[j].uid;
		}
	}
}

void deBGA_INDEX::expand_seed(std::vector<vertex_U>& vertexu_v, std::vector<UNI_SEED>& uniseed_v, random_data* rand_buff)
{
	for(uint32_t i = 0; i < vertexu_v.size(); i++){
		vertex_U & vertexU = vertexu_v[i];
		if(vertexU.pos_n > POS_N_MAX){
			if(vertexU.pos_n > POS_N_MAX_LEVEL2)  return;
			for(uint32_t ri = 0; ri < RANDOM_NUM; ri++){
				int rand_num;
				random_r(rand_buff, &rand_num);
				uint32_t m = (rand_num % vertexU.pos_n);
				uniseed_v.emplace_back();
				UNI_SEED & uniseed = uniseed_v.back();
				uniseed.seed_id = i;
				uniseed.read_begin = vertexU.read_pos;
				uniseed.read_end = vertexU.read_pos + vertexU.length1 - 1;
				uniseed.ref_begin =  buffer_p[m + buffer_pp[vertexU.uid]] + vertexU.uni_pos_off - 1;
				uniseed.ref_end = uniseed.ref_begin + vertexU.length2 - 1;
				uniseed.cov = vertexU.cov;
			}
		}else{
			for(uint32_t m = 0; m < vertexU.pos_n; m++){
				uniseed_v.emplace_back();
				UNI_SEED & uniseed = uniseed_v.back();
				uniseed.seed_id = i;
				uniseed.read_begin = vertexU.read_pos;
				uniseed.read_end = vertexU.read_pos + vertexU.length1 - 1;
				uniseed.ref_begin =  buffer_p[m + buffer_pp[vertexU.uid]] + vertexU.uni_pos_off - 1;
				uniseed.ref_end = uniseed.ref_begin + vertexU.length2 - 1;
				uniseed.cov = vertexU.cov;
			}
		}
	}
}
//
////store seed in uintig into seed in reference
//void deBGA_INDEX::expand_seed(std::vector<vertex_U>& vertexu_v, std::vector<UNI_SEED>& uniseed_v, unsigned int rand_num)
//{
//	for(uint32_t i = 0; i < vertexu_v.size(); i++){
//		vertex_U & vertexU = vertexu_v[i];
//		if(vertexU.pos_n < POS_N_MAX){
//			for(uint32_t m = 0; m < vertexU.pos_n; m++){
//				uniseed_v.emplace_back();
//				UNI_SEED & uniseed = uniseed_v.back();
//				uniseed.seed_id = i;
//				uniseed.read_begin = vertexU.read_pos;
//				uniseed.read_end = vertexU.read_pos + vertexU.length1 - 1;
//				uniseed.ref_begin =  buffer_p[m + buffer_pp[vertexU.uid]] + vertexU.uni_pos_off - 1;
//				uniseed.ref_end = uniseed.ref_begin + vertexU.length2 - 1;
//				uniseed.cov = vertexU.cov;
//			}
//		}
//		else if(vertexU.pos_n < POS_N_MAX_LEVEL2){
//			int AVG_step_len = (vertexU.pos_n * 2) / RANDOM_NUM;
//			uint32_t m =  (rand() % vertexU.pos_n) ;
//			for(uint32_t ri = 0; ri < RANDOM_NUM; ri++){
//				rand_num = rand_r(&rand_num);
//				int c_rand = (rand_num % AVG_step_len) + 1;
//				m += c_rand;
//				if(m > vertexU.pos_n) m -= vertexU.pos_n;
//				//fprintf(stderr, "%d:%d:%d:%d:%d[%d]\t", ri, m, c_rand, vertexU.pos_n, AVG_step_len, rand_num);
//				uniseed_v.emplace_back();
//				UNI_SEED & uniseed = uniseed_v.back();
//				uniseed.seed_id = i;
//				uniseed.read_begin = vertexU.read_pos;
//				uniseed.read_end = vertexU.read_pos + vertexU.length1 - 1;
//				uniseed.ref_begin =  buffer_p[m + buffer_pp[vertexU.uid]] + vertexU.uni_pos_off - 1;
//				uniseed.ref_end = uniseed.ref_begin + vertexU.length2 - 1;
//				uniseed.cov = vertexU.cov;
//			}
//		}
//		else if(vertexU.pos_n < POS_N_MAX_LEVEL3){
//			uint32_t m =  (rand() % vertexU.pos_n) ;
//			for(uint32_t ri = 0; ri < RANDOM_NUM; ri++){
//				m += 20; if(m > vertexU.pos_n) m -= vertexU.pos_n;
//				//fprintf(stderr, "XXX%d:%d:%d:%d:%d\t", ri, m, 20, vertexU.pos_n, 20 );
//				uniseed_v.emplace_back();
//				UNI_SEED & uniseed = uniseed_v.back();
//				uniseed.seed_id = i;
//				uniseed.read_begin = vertexU.read_pos;
//				uniseed.read_end = vertexU.read_pos + vertexU.length1 - 1;
//				uniseed.ref_begin =  buffer_p[m + buffer_pp[vertexU.uid]] + vertexU.uni_pos_off - 1;
//				uniseed.ref_end = uniseed.ref_begin + vertexU.length2 - 1;
//				uniseed.cov = vertexU.cov;
//			}
//		}
//	}
//}

void deBGA_INDEX::get_refseq(uint8_t *ref, uint32_t len, uint32_t start)
{
	uint32_t m;

    for (m = 0; m < len; ++m)
    {
        ref[m] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

void deBGA_INDEX::get_ref_jump(uint8_t *ref, uint32_t start1, uint32_t len1, uint32_t start2, uint32_t len2)
{
	uint32_t m;
	uint32_t k = 0;

    for (m = 0; m < len1; ++m)
    {
        ref[k++] = (buffer_ref_seq[(m + start1) >> 5] >> ((31 - ((m + start1) & 0X1f)) << 1)) & 0X3;
    }

    for (m = 0; m < len2; ++m)
    {
        ref[k++] = (buffer_ref_seq[(m + start2) >> 5] >> ((31 - ((m + start2) & 0X1f)) << 1)) & 0X3;
    }
}

void deBGA_INDEX::get_ref_onebyone(uint8_t *ref, uint32_t start, uint32_t len, uint32_t pre_pos)
{
	uint32_t m;
	uint32_t k = pre_pos;
    for (m = 0; m < len; ++m)
    {
        ref[k++] = (buffer_ref_seq[(m + start) >> 5] >> ((31 - ((m + start) & 0X1f)) << 1)) & 0X3;
    }
}

void deBGA_INDEX::free_memory(){
	if(buffer_ref_seq 	!= NULL){ free(buffer_ref_seq); buffer_ref_seq = NULL; 	} //original reference
	if(buffer_seq 		!= NULL){ free(buffer_seq); 	buffer_seq = NULL; 		}// UNITIG sequence
	if(buffer_seqf != NULL)		{ free(buffer_seqf); 	buffer_seqf = NULL; }//start offset of [UINTIG] in [buffer_seq]
	if(buffer_off_g != NULL)	{ free(buffer_off_g); 	buffer_off_g = NULL; }//start offset of [KMER] in [buffer_seq]
	if(buffer_p != NULL)		{ free(buffer_p); 		buffer_p = NULL; }//start offset of [UINTIG] in [buffer_ref_seq]
	if(buffer_pp != NULL)		{ free(buffer_pp); 		buffer_pp = NULL; }//start index of [UINTIG] in [buffer_p], used to search offset of [UINTIG] in [buffer_ref_seq]
	if(buffer_hash_g != NULL)	{ free(buffer_hash_g); 	buffer_hash_g = NULL; }//index for the first 14 base pair, used to search kmer in the hash table
	if(buffer_kmer_g != NULL)	{ free(buffer_kmer_g); 	buffer_kmer_g = NULL; }//used to store other part of a kmer(except the first 14 byte)
}

#define CHR_SEARCH_IDX_STEP 0x4000 //16384 = 2 ^ 14
void deBGA_INDEX::building_chr_index(){
	chr_search_index_size = (reference_len >> 14) + 2;
	chr_search_index = (uint32_t *)xcalloc(chr_search_index_size, sizeof(uint32_t));
	uint32_t pos_index_size = 0;
	for(int i = 0; i < chr_file_n; i++){
		int pos_index = chr_end_n[i] / CHR_SEARCH_IDX_STEP;
		while(pos_index >= pos_index_size){
			chr_search_index[pos_index_size++] = i;
		}
	}
	chr_search_index[pos_index_size] = chr_file_n;// store final block
}

//get the chr_ID that a 'position' belong to.
int deBGA_INDEX::get_chromosome_ID(uint32_t position)
{
	int file_n = 0;
	int pos_index = position / CHR_SEARCH_IDX_STEP;
	int low = chr_search_index[pos_index];
	int high = chr_search_index[pos_index + 1];
	int mid;
	int pos = position + 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(pos < (chr_end_n[mid] - 1))
		{
			high = mid - 1;
		}
		else if(pos > (chr_end_n[mid] - 1))
		{
			low = mid + 1;
		}
		else
		{
			return mid;
		}
		file_n = low;
	}
	return file_n;
}

void deBGA_INDEX::building_bam_header(){
	header = (bam_hdr_t *)xcalloc(1, sizeof(bam_hdr_t));
	header->n_targets = chr_file_n;
	header->target_name = (char **)xcalloc(header->n_targets, sizeof(char *));//sheader->n_targetschr_names;
	header->target_len = (uint32_t *)xcalloc(header->n_targets, sizeof(uint32_t));

	header->target_name[0] = chr_names[0];
	header->target_len[0] = chr_end_n[0];
	for(int i = 1; i < header->n_targets; i++){
		header->target_name[i] = chr_names[i];
		header->target_len[i] = chr_end_n[i] - chr_end_n[i - 1];
	}

	//0_1_136054_1148_INS
	char target_name_tmp[1024];
	for(int i = 0; i < header->n_targets; i++){
		strcpy(target_name_tmp, chr_names[i]);
		char * token = strtok(target_name_tmp, "_"); 	int id = atoi(token);
		token = strtok(NULL, "_"); 				int char_ID = bam_name2id(ori_header, token);		//char name
		token = strtok(NULL, "_"); 				uint32_t pos = atoi(token);

		token = strtok(NULL, "_"); 				int region_len = atoi(token);
		token = strtok(NULL, "_"); 				char * sv_type = token;

		token = strtok(NULL, "_"); 				int break_point1_pos =  atoi(token);//position of break point 1 in original reference
		token = strtok(NULL, "_"); 				int break_point2_pos =  atoi(token);//position of break point 2 in original reference
		token = strtok(NULL, "_"); 				int ed_pos =  atoi(token);//ending position of new reference in original reference
		token = strtok(NULL, "_"); 				char* vcf_id = token;//a name of vcf record in vcf file, like  "pbsv.INS.2"

		sv_info.emplace_back(id, char_ID, pos ,
				region_len, sv_type,
				break_point1_pos ,break_point2_pos ,ed_pos , vcf_id);
	}
}

int deBGA_index_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:	%s fc_index <ref.fa> <index_route>\n", PACKAGE_NAME);
	fprintf(stderr, "		build deBGA index file with default 22-kmer. You can get more deBGA information from https://github.com/HongzheGuo/deBGA");
	fprintf(stderr, "\n");
	return 1;
}

int get_bin_dir(char *bin, char *dir)
{
	char *end = strrchr(bin, '/');
	if (end == NULL)
		return 1;

	bin[end-bin] = '\0';
	strcpy(dir, bin);
	return 0;
}

int build_deBGA_index_core(char *dir, char *ref_fa, char *index_route)
{
	char cmd[1024];
	sprintf(cmd, "%sdeBGA index %s %s", dir, ref_fa, index_route);
	fprintf(stderr, "[deBGA_index] Executing deBGA index ...\n");
	if (system(cmd) != 0)
	{
		fprintf(stderr, "\n[deBGA_index] Indexing undoing, deBGA index exit abnormally. \n");
		exit(1);
	}
	fprintf(stderr, "[deBGA_index] Done!\n");

    return 1;
}

int build_deBGA_index(int argc, char *argv[])
{
	// clock_t t = clock();
	char *ref_fa = 0;
	char *index_route = 0;

    char dir[1024];
    char desalt_path[1024];
    int r;

	if(optind + 2 > argc)	return deBGA_index_usage();

	ref_fa = strdup(argv[optind]);
	index_route = strdup(argv[optind + 1]);

    r = readlink("/proc/self/exe", dir, 2048);
    if (r < 0 || r >= 2048)
    {
        fprintf(stderr, "Failed, could not find the program of deSALT!\n");
    }
    dir[r] = '\0';

    if (!get_bin_dir(dir, desalt_path))
    {
        strcat(desalt_path, "/");
    }

    char path_deBGA[1024];
    strcpy(path_deBGA, desalt_path);
    strcat(path_deBGA, "deBGA");

    if((access(path_deBGA, F_OK)) == -1)
    {
        fprintf(stderr, "[Wrong!] %s not exist, please check!\n", path_deBGA);
        exit(1);
    }

    build_deBGA_index_core(desalt_path, ref_fa, index_route);

	return 0;
}

