// #pragma once

#ifndef __XDATA__
#define __XDATA__

class XData
{
private:
	void* data;
	int* size, type;
public:
	enum TYPE { CHAR, UCHAR, SHORT, USHORT, INT, UINT, LONG, ULONG, FLOAT32, DOUBLE };
	XData();
	~XData();
	int setData(void* data, int* size, int type);
	void* getData();
	int* getSize();
	int getType();
};

#endif
