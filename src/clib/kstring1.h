/*
 * kstring1.h
 *
 *  Created on: 2019年12月28日
 *      Author: fenghe
 */

#ifndef LIB_KSTRING1_H_
#define LIB_KSTRING1_H_

	//string
	void ks_resize(kstring_t *s, size_t size);//resize a string
	int kputsn(const char *p, int l, kstring_t *s);//store a string, == memcat
	void str_enlarge(kstring_t *s, int l);//enlarge a string by l
	int kputs(const char *p, kstring_t *s);//store a string
	int kputc(int c, kstring_t *s);//store a char
	int kputw(int c, kstring_t *s);//store an int
	int kputuw(unsigned c, kstring_t *s);//store an unsigned-int
	int kputl(long c, kstring_t *s);//store a long int
	void kstrcpy(kstring_t *s, const char *st, const char *en);
	int ksprintf(kstring_t *s, const char *fmt, ...);
	void sprintf_lite(kstring_t *s, const char *fmt, ...);
	int big_or_small_endian();

#endif /* LIB_KSTRING1_H_ */
