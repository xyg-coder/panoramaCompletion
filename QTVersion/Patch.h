#pragma once
#include <Windows.h>
class Patch
{
public:
	Patch(void);
	~Patch(void);

public:
	int _x;
	int _y;
	double *_patchWHT;
	int _node_index;
};

