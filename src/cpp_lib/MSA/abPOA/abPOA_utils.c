/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include "abPOA_utils.h"

#include "abPOA_kseq.h"

/********************
 * System utilities *
 ********************/

int err_func_printf(const char *func, const char *format, ...)
{
    fprintf(stderr, "[%s] ", func);
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stderr, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stderr)", strerror(saveErrno));
	return done;
}


int stdout_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

void err_fgets(char *buff, size_t s, FILE *fp)
{
    if (fgets(buff, s, fp) == NULL) {
        err_fatal_simple("fgets error.\n");
    }
}

/*********
 * alloc *
 *********/
void *err_malloc(const char *func, size_t s)
{
    void *ret = (void*)malloc(s);
    if (ret == NULL) err_fatal_core(func, "Malloc fail!\nSize: %lld\n", s);
    else return ret;
}

void *err_calloc(const char *func, size_t n, size_t s)
{
    void *ret = (void*)calloc(n, s);
    if (ret == NULL) err_fatal_core(func, "Calloc fail!\nN: %d\tSize: %lld\n", n, s);
    else return ret;
}

void *err_realloc(const char *func, void *p, size_t s)
{
    void *ret = (void*)realloc(p, s);
    if (ret == NULL) err_fatal_core(func, "Realloc fail!\nSize: %lld\n", s);
    else return ret;
}

/*********
 * Timer *
 *********/
void usr_sys_cputime(double *usr_t, double *sys_t)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
    *usr_t = r.ru_utime.tv_sec + 1e-6 * r.ru_utime.tv_usec;
	*sys_t = r.ru_stime.tv_sec + 1e-6 * + r.ru_stime.tv_usec;
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

void get_cur_time(const char *prefix)
{
    time_t now = time(0);
    struct tm ts; char buf[1024];
    ts = *localtime(&now);
    err_printf("[%s] ", prefix);
    strftime(buf, sizeof(buf), "%Y-%m-%d-%s", &ts);
}

void print_format_time(FILE *out)
{
    time_t rawtime;
    struct tm *info;
    char buffer[80];

    time(&rawtime);
    info = localtime( &rawtime );
    strftime(buffer,80,"%m-%d-%Y %X", info);
    fprintf(out, "=== %s === ", buffer);
}

int err_func_format_printf(const char *func, const char *format, ...)
{
    print_format_time(stderr);
    fprintf(stderr, "[%s] ", func);
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stderr, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stderr)", strerror(saveErrno));
	return done;
}
