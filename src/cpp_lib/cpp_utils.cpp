/*
 * cpp_utils.cpp
 *
 *  Created on: 2021年10月27日
 *      Author: fenghe
 */

#include "cpp_utils.hpp"
#include <cstring>

void split_string(std::vector<std::string> &item_value, char * temp, char * split_line, const char *split_str){
	strcpy(temp, split_line);
	char * token_value = NULL; item_value.clear();
	token_value = strtok(temp, split_str); item_value.emplace_back(token_value);
	for(int item_idx = 1; item_idx < 10000000; item_idx++){
		token_value = strtok(NULL, split_str);
		if(token_value == NULL)
			break;
		item_value.emplace_back(token_value);
	}
}
