
#ifndef __THREADUTIL__
#define __THREADUTIL__

void *threaddata_create(int nbthreads, long size0);
void *threaddata_release(void *_threaddata);
int threaddata_get_nbthreads(void *_threaddata);
long threadata_get_size(void *_threaddata, int no);
long threaddata_get_offset(void *_threaddata, int no);

#endif