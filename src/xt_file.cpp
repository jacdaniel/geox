#include <malloc.h>
#include <xt_file.h>


static void endian_short_swap(short *data, long size)
{
    for (long add=0; add<size; add++)
    {
        data[add] = ((data[add] & 0xff) << 8) | ((data[add] & 0xff00) >> 8);
    }
}

#ifndef SAWP
#define SWAP(_a_, _b_, _temp_) { \
	_temp_ = _b_; \
	_b_ = _a_; \
	_a_ = _temp_; }
#endif 

static void endianSwap(void *_data, long size, int length)
{
  long i, j;
  unsigned char temp, *data = (unsigned char*)_data;
  
  for (i=0; i<size; i++)
    for (j=0; j<length/2; j++)
      SWAP(data[i*length+j],  data[i*length+(length-1-j)], temp)
}



XT_FILE::XT_FILE()
{
    this->pFile = NULL;
}

XT_FILE::~XT_FILE()
{
    if ( this->pFile ) fclose(this->pFile);
}

void XT_FILE::openForRead(char *filename)
{
    this->pFile = fopen(filename, "r");
    header_read();
}

void XT_FILE::openForWrite(char *filename)
{
    this->pFile = fopen(filename, "r+");
    header_read();
}

void XT_FILE::header_read()
{
    fseek((FILE*)this->pFile, 4, SEEK_SET);
    fread((void*)&this->dimx, sizeof(int), 1, this->pFile); endianSwap(&this->dimx, 1, sizeof(int));
    fread((void*)&this->dimy, sizeof(int), 1, this->pFile); endianSwap(&this->dimy, 1, sizeof(int));
    fread((void*)&this->dimz, sizeof(int), 1, this->pFile); endianSwap(&this->dimz, 1, sizeof(int));
    fseek((FILE*)this->pFile, 32, SEEK_SET);
    fread((void*)&this->format, sizeof(int), 1, this->pFile); endianSwap(&this->format, 1, sizeof(int));
    fread((void*)&this->_sizeof, sizeof(int), 1, this->pFile); endianSwap(&this->_sizeof, 1, sizeof(int));
    fseek((FILE*)this->pFile, 72, SEEK_SET);
    fread((void*)&this->offset, sizeof(int), 1, this->pFile); endianSwap(&this->offset, 1, sizeof(int));
    // fprintf(stderr, "offset :%d\n", this->offset);
}

void XT_FILE::header_write(int format, int dimx, int dimy, int dimz)
{
    int tmp;
    fseek(this->pFile, 4, SEEK_SET);
    tmp = dimx; endianSwap(&tmp, 1, sizeof(int)); int a = fwrite((void*)&tmp, sizeof(int), 1, this->pFile);
    tmp = dimy; endianSwap(&tmp, 1, sizeof(int)); fwrite((void*)&tmp, sizeof(int), 1, pFile);
    tmp = dimz; endianSwap(&tmp, 1, sizeof(int)); fwrite((void*)&tmp, sizeof(int), 1, pFile);
    fseek(pFile, 32, SEEK_SET);
    tmp = format; endianSwap(&tmp, 1, sizeof(int)); fwrite((void*)&tmp, sizeof(int), 1, pFile);
    tmp = SIZEOF[format]; endianSwap(&tmp, 1, sizeof(int)); fwrite((void*)&tmp, sizeof(int), 1, pFile);
}

void XT_FILE::createNew(char *filename_src, char *filename_dst)
{
    int offset0;
    FILE *pfile_src = fopen(filename_src, "r");
    fseek(pfile_src, 72, SEEK_SET);
    fread(&offset0, sizeof(int), 1, pfile_src);
    endianSwap(&offset0, 1, sizeof(int));
    char *h = (char*)calloc(offset0, sizeof(char));
    fseek(pfile_src, 0, SEEK_SET);
    fread(h, sizeof(char), offset0, pfile_src);
    fclose(pfile_src);
    this->pFile = fopen(filename_dst, "w+");
    fwrite(h, sizeof(char), offset0, this->pFile);
    header_read();
    free(h);
}


void XT_FILE::createNew(char *filename_src, char *filename_dst, int dimx, int dimy, int dimz, int format)
{
    createNew(filename_src, filename_dst);
    header_write(format, dimx, dimy, dimz); 
}

int XT_FILE::get_dimx()
{
    return this->dimx;
}

int XT_FILE::get_dimy()
{
    return this->dimy;
}

int XT_FILE::get_dimz()
{
    return this->dimz;
}

int XT_FILE::get_format()
{
    return this->format;
}

int XT_FILE::get_sizeof()
{
    return this->_sizeof;
}



float XT_FILE::inlineRead(int no, void *data)
{
    long sizeof0 = this->_sizeof/8;
    size_t offset0 = (size_t)this->offset + (size_t)this->dimx*this->dimy*sizeof0*no;
    long size0 = (long)dimx*dimy;
    fseek(this->pFile, offset0, SEEK_SET);
    fread(data, sizeof0, size0, this->pFile);
    if ( this->format == 2 || this->format == 8)
        endianSwap(data, size0, sizeof0);
    return 0.0f;
}


float XT_FILE::inlineManyRead(int no, int nbre, void *data)
{
    long sizeof0 = this->_sizeof/8;
    size_t offset0 = (size_t)this->offset + (size_t)this->dimx*this->dimy*sizeof0*no;
    long size0 = (long)dimx*dimy*nbre;
    fseek(this->pFile, offset0, SEEK_SET);
    fread(data, sizeof0, size0, this->pFile);
    if ( this->format == 2 || this->format == 8)
        endianSwap(data, size0, sizeof0);
    return 0.0f;
}

float XT_FILE::inlineManyReadShort(int no, int nbre, void *data)
{
	long sizeof0 = this->_sizeof/8;
	size_t offset0 = (size_t)this->offset + (size_t)this->dimx*this->dimy*sizeof0*no;
	fseek(this->pFile, offset0, SEEK_SET);

	if ( this->format == 2 )
	{
		long size0 = (long)dimx*dimy*nbre;
		// fseek(this->pFile, offset0, SEEK_SET);
		fread(data, sizeof0, size0, this->pFile);
		endianSwap(data, size0, sizeof0);
		return 0.0f;
	}
	else
	{
		size_t size0 = (size_t)this->dimx*this->dimy;
		void *inlinedata = calloc((size_t)this->dimx*this->dimy, sizeof0);
		if ( inlinedata == NULL ) { return 0.0f; }
		for (int z=0; z<nbre; z++)
		{
			fread(inlinedata, sizeof0, size0, this->pFile);
			if ( this->format == 1 )
			{

			}
			else if ( this->format == 0 )
			{
				endianSwap(inlinedata, size0, sizeof0);
				for (size_t add=0; add<(size_t)this->dimx*this->dimy; add++)
				{
					((short*)data)[(size_t)z*this->dimx*this->dimy+add] = (short)((float*)inlinedata)[add];
				}
			}
		}
		if ( inlinedata != NULL ) { free(inlinedata); inlinedata = NULL; }
		return 0.0f;
	}
}


float XT_FILE::inlineWrite(int no, void *data)
{
    long sizeof0 = this->_sizeof/8;
    size_t offset0 = (size_t)this->offset + (size_t)this->dimx*this->dimy*sizeof0*no;
    long size0 = (long)dimx*dimy;

    if ( this->format == 2 || this->format == 8)
         endianSwap(data, size0, sizeof0); 
    fseek(this->pFile, offset0, SEEK_SET);
    fwrite(data, sizeof0, size0, this->pFile);
    if ( this->format == 2 || this->format == 8)
         endianSwap(data, size0, sizeof0);
    return 0.0f;
}

float  XT_FILE::inlineIncrementalWrite(int no, void *data)
{
long sizeof0 = this->_sizeof/8;
    size_t offset0 = (size_t)this->offset;// + (size_t)this->dimx*this->dimy*sizeof0*no;
    long size0 = (long)dimx*dimy;

    if ( this->format == 2 || this->format == 8)
         endianSwap(data, size0, sizeof0);
    if ( no == 0 )
        fseek(this->pFile, offset0, SEEK_SET);
    fwrite(data, sizeof0, size0, this->pFile);
    if ( this->format == 2 || this->format == 8)
         endianSwap(data, size0, sizeof0);
    return 0.0f;
}



float XT_FILE::inlineManyWrite(int no, int nbre, void *data)
{
    int sizeof0 = this->_sizeof/8;
    size_t offset0 = (size_t)this->offset + (size_t)this->dimx*this->dimy*sizeof0*no;
    long size0 = (long)dimx*dimy*nbre;
    if ( this->format == 2 || this->format == 8)
        endianSwap(data, size0, sizeof0);    
    fseek(this->pFile, offset0, SEEK_SET);    
    fwrite(data, sizeof0, size0, this->pFile);
    if ( this->format == 2 || this->format == 8)
        endianSwap(data, size0, sizeof0);
    return 0.0f;
}
