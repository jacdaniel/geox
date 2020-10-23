
#include <XData.h>

XData::XData()
{
}

XData::~XData()
{
}

int XData::setData(void* data, int* size, int type)
{
	this->data = data;
	this->size = size;
	this->type = type;
	return 0;
}

void* XData::getData()
{
	return this->data;
}

int* XData::getSize()
{
	return this->size;
}

int XData::getType()
{
	return this->type;
}