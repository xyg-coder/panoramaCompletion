#pragma once
#include <Windows.h>
#include "Patch.h"
class Node
{
public:
	Node(void);
	~Node(void);

public:
	Patch *_patchArray;
	int _patchArray_count;
	bool _flag;
};

