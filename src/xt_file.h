

#ifndef __XT_FILE__
#define __XT_FILE__

#include <stdio.h>

class XT_FILE
{
    enum FORMAT{CHAR, UCHAR, INT16, UINT16, INT32, UINT32, INT64, UINT64, FLOAT, DOUBLE};
    const int SIZEOF[9] = {1,1,16,16,32,32,64,32,64};
    public:
    XT_FILE();
    ~XT_FILE();
    void openForRead(char *filename);
    void openForWrite(char *filename);
    void createNew(char *filename_src, char *filename_dst);
    void createNew(char *filename_src, char *filename_dst, int dimx, int dimy, int dimz, int format);

    int get_dimx();
    int get_dimy();
    int get_dimz();
    int get_format();
    int get_sizeof();
    float inlineRead(int no, void *data);
    float inlineManyRead(int no, int nbre, void *data);
    float inlineManyReadShort(int no, int nbre, void *data);
    float inlineWrite(int no, void *data);
    float inlineIncrementalWrite(int no, void *data);
    float inlineManyWrite(int no, int nbre, void *data);


    private:
        int dimx, dimy, dimz, format, _sizeof;
        int offset;
        FILE *pFile;

        void header_read();
        void header_write(int format, int dimx, int dimy, int dimz);




};



#endif
