#include <malloc.h>
#include <algorithm>
#include <ThreadUtil.h>


typedef struct _THREADDATA {
  int nbthreads;
  long *offset;
  long *size;
}THREADDATA;

void *threaddata_create(int nbthreads, long size0)
{
  THREADDATA *data = (THREADDATA*)calloc(1, sizeof(THREADDATA));
  data->nbthreads = nbthreads;
  if (data->nbthreads > size0)
  {
    data->nbthreads = size0;
  }
  long s = ( size0 - 1 ) / data->nbthreads + 1;
  data->offset = (long*)calloc(data->nbthreads, sizeof(long));
  data->size = (long*)calloc(data->nbthreads, sizeof(long));
  for (int i = 0; i < data->nbthreads; i++)
  {
    data->offset[i] = std::min((long)i * s, size0-1);
    data->size[i] = std::max(std::min(s, size0 - s * i), (long)0);
  }
  return data;
}

void *threaddata_release(void *_threaddata)
{
  THREADDATA *threaddata = (THREADDATA*)_threaddata;
  if ( threaddata == NULL ) return NULL;
  free(threaddata->offset);
  free(threaddata->size);
  free(threaddata);
  return nullptr;
}

int threaddata_get_nbthreads(void *_threaddata)
{
  THREADDATA *threaddata = (THREADDATA*)_threaddata;
  if ( threaddata == nullptr ) return 0;
  return threaddata->nbthreads;
}

long threadata_get_size(void *_threaddata, int no)
{
  THREADDATA *threaddata = (THREADDATA*)_threaddata;
  if ( threaddata == nullptr || no >= threaddata->nbthreads ) return 0;
  return threaddata->size[no];
}

long threaddata_get_offset(void *_threaddata, int no)
{
  THREADDATA *threaddata = (THREADDATA*)_threaddata;
  if ( threaddata == nullptr || no >= threaddata->nbthreads ) return 0;
  return threaddata->offset[no];
}