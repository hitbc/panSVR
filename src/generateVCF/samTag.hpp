/*
 * samTag.hpp
 *
 *  Created on: 2021年5月31日
 *      Author: fenghe
 */

#ifndef GENERATEVCF_SAMTAG_HPP_
#define GENERATEVCF_SAMTAG_HPP_

struct VCF_TAG{
	char AS[4] = "AS";	//alignment score
	char OS[4] = "OS"; //original alignment score
	char OA[4] = "OA"; //original alignment
	char CS[4] = "CS"; //chaining score
	char SV[4] = "SV"; // SV information
	char MV[4] = "MV"; //mate SV information
	char XA[4] = "XA"; //ALT alignment
	char RC[4] = "RC"; //Read commend at signal STEP
};

#ifndef vcf_tag
	static VCF_TAG vcf_tag;
#endif

#endif /* GENERATEVCF_SAMTAG_HPP_ */
