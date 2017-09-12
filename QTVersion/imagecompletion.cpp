#include "imagecompletion.h"

ImageCompletion::ImageCompletion(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	_WHT = NULL;
	_imgRGB = NULL;
	_imgMask = NULL;
	_mask = NULL;
	_imgwidth = 0;
	_imgheight = 0;
	_ukpixels_count = 0;
	_imgYCbCr = NULL;
	_patchsize = 8;
	_imgPatch = NULL;
	_patch_count = 0;
	_x0 = _y0 = 0;
	_x1 = _y1 = 0;
	_nodeArray = NULL;
	_nodes_count = 0;
	_index_frompixeltopatch = NULL;
	_flagmatched = NULL;
	_deltx = NULL;
	_delty = NULL;
	_minoffsetx = _minoffsety = _maxoffsetx = _maxoffsety = 0;
	_histogram = NULL;
	_peaks_count = 0;
	_peaksvalue = NULL;
	_peaksx = NULL;
	_peaksy = NULL;
	_candidates_count = 0;
	_offsetx_candidates = NULL;
	_offsety_candidates = NULL;
	_optimization_result = NULL;
	_deltx2 = NULL;
	_delty2 = NULL;

	ui.setupUi(this);
	connect(ui.actionOpen_Image, SIGNAL(triggered()), this, SLOT(openRGBImage()));
	connect(ui.actionOpen_Mask, SIGNAL(triggered()), this, SLOT(openMaskImage()));
	connect(ui.actionTo_YCbCr, SIGNAL(triggered()), this, SLOT(toYCbCrChannel()));
	connect(ui.actionCompute_WHT, SIGNAL(triggered()), this, SLOT(generateWHTVector()));
	connect(ui.actionGenerate_KD_Tree, SIGNAL(triggered()), this, SLOT(generateKDTree()));
	connect(ui.actionStatistics_of_Patch_Offset, SIGNAL(triggered()), this, SLOT(statisticsPatchOffset()));
	connect(ui.actionGCoptimization, SIGNAL(triggered()), this, SLOT(gcOptimization()));
	connect(ui.actionCompletion, SIGNAL(triggered()), this, SLOT(completionImage()));
}

ImageCompletion::~ImageCompletion()
{

}

void ImageCompletion::openRGBImage()
{
	QString filename;
	filename = QFileDialog::getOpenFileName(this,
		tr("Open Image"),
		"",
		tr("Images (*.png *.bmp *.jpg *.tif *.GIF )"));
	if(filename.isEmpty())
	{
		return;
	}
	else
	{
		if(_imgRGB != NULL)
		{
			delete _imgRGB;
		}
		_imgRGB = new QImage;

		if(!( _imgRGB->load(filename))) //¼ÓÔØÍ¼Ïñ
		{
			QMessageBox::information(this,
				tr("Open Fail"),
				tr("Open Fail!"));
			delete _imgRGB;
			return;
		}
		_imgwidth = _imgRGB->width();
		_imgheight = _imgRGB->height();
		ui.label->setPixmap(QPixmap::fromImage(*_imgRGB));
	}
}

void ImageCompletion::openMaskImage()
{
	QString filename;
	filename = QFileDialog::getOpenFileName(this,
		tr("Open Image"),
		"",
		tr("Images ( *.jpg )"));
	if(filename.isEmpty())
	{
		return;
	}
	else
	{
		if(_imgMask != NULL)
		{
			delete _imgMask;
		}
		_imgMask = new QImage;

		if(!( _imgMask->load(filename))) //¼ÓÔØÍ¼Ïñ
		{
			QMessageBox::information(this,
				tr("Open Fail"),
				tr("Open Fail!"));
			delete _imgMask;
			_imgMask = NULL;
			return;
		}

		if(_imgMask->bitPlaneCount() != 8)
		{
			QMessageBox::information(this,
				tr("Open Fail"),
				tr("The Mask Image Depth should be 8!"));
			delete _imgMask;
			_imgMask = NULL;
			return;
		}

		if(_imgwidth == _imgMask->width() && _imgheight == _imgMask->height())
		{
			QImage *img = new QImage;
			img->load("E:\\pano\\3-2.jpg");

			if(_mask != NULL)
			{
				delete []_mask;
			}
			_mask = new int[_imgheight * _imgwidth];

			unsigned char *data = img->bits();
			unsigned char *data2 = _imgMask->bits();

			_ukpixels_count = 0;
			for(int i = 0; i < _imgheight; i++)
			{
				data = img->scanLine(i);
				data2 = _imgMask->scanLine(i); 

				for(int j = 0; j < _imgwidth; j++)
				{
					if(data[j] < 200)
					{
						if(data2[j] < 200)
						{
							_mask[i * _imgwidth + j] = 0; 
						}
						else
						{
							_mask[i * _imgwidth + j] = 1;
						}
						_ukpixels_count ++;
						if(_ukpixels_count == 1)
						{
							_upixels = new UKPixels[_ukpixels_count];
							_upixels[_ukpixels_count - 1].x = j;
							_upixels[_ukpixels_count - 1].y = i;
							if(data2[j] < 200)
							{
								_upixels[_ukpixels_count - 1].flag_boundary = false;
							}
							else
							{
								_upixels[_ukpixels_count - 1].flag_boundary = true;
							}
						}
						else
						{
							UKPixels *temp = new UKPixels[_ukpixels_count];
							for(int m = 0; m < _ukpixels_count - 1; m++)
							{
								temp[m] = _upixels[m]; 
							}
							temp[_ukpixels_count - 1].x = j;
							temp[_ukpixels_count - 1].y = i;
							
							if(data2[j] < 200)
							{
								temp[_ukpixels_count - 1].flag_boundary = false;
							}
							else
							{
								temp[_ukpixels_count - 1].flag_boundary = true;
							}

							delete []_upixels;
							_upixels = temp;
							temp = NULL;
						}
					}
					else
					{
						_mask[i * _imgwidth + j] = 1;
					}
				}
			}
			delete img;
		}
		else
		{
			QMessageBox::information(this,
				tr("Open Fail"),
				tr("Open Fail! The Mask Image should be as large as the RGB Image!"));
			delete _imgMask;
			_imgMask = NULL;
			return;
		}
	}
}

void ImageCompletion::toYCbCrChannel()
{
	unsigned char *data = _imgRGB->bits();
	unsigned char R;
	unsigned char G;
	unsigned char B;
	if(_imgYCbCr != NULL)
	{
		delete []_imgYCbCr;
	}
	_imgYCbCr = new double[_imgwidth * _imgheight * 3];

	for(int i = 0; i < _imgheight; i++)
	{
		data = _imgRGB->scanLine(i);
		for(int j = 0; j < _imgwidth; j++)
		{
			R = data[j * 4 + 2];
			G = data[j * 4 + 1];
			B = data[j * 4];

			_imgYCbCr[_imgwidth * _imgheight * 0 + i * _imgwidth + j] = 16 + (0.256789 * R + 0.504129 * G + 0.097906 * B);
			_imgYCbCr[_imgwidth * _imgheight * 1 + i * _imgwidth + j] = 128 + (-0.148223 * R - 0.290992 * G + 0.439215 * B);
			_imgYCbCr[_imgwidth * _imgheight * 2 + i * _imgwidth + j] = 128 + (0.439215 * R - 0.367789 * G - 0.071426 * B);
		}
	}
}

void ImageCompletion::findRectangle()
{
	int a1 = -10000;
	int a2 = 10000;
	int b1 = -10000;
	int b2 = 10000;

	for(int i = 0; i < _ukpixels_count; i++)
	{
		if(_upixels[i].x > a1)
			a1 = _upixels[i].x;
		if(_upixels[i].x < a2)
			a2 = _upixels[i].x;
		if(_upixels[i].y > b1)
			b1 = _upixels[i].y;
		if(_upixels[i].y < b2)
			b2 = _upixels[i].y;
	}
	_x0 = ((a2 - 1 * (a1 - a2 + 1)) >= 0) ? a2 - 1 * (a1 - a2 + 1) : 0;
	_x1 = ((a1 + 1 * (a1 - a2 + 1)) < _imgwidth) ? a1 + 1 * (a1 - a2 + 1) : _imgwidth - 1;
	_y0 = ((b2 - 1 * (b1 - b2 + 1)) >= 0) ? b2 - 1 * (b1 - b2 + 1) : 0;
	_y1 = ((b1 + 1 * (b1 - b2 + 1)) < _imgheight) ? b1 + 1 * (b1 - b2 + 1) : _imgheight - 1;
}

bool ImageCompletion::judge1(int x, int y, int size)
{
	bool flag = true;
	for(int i = y; i < y + _patchsize; i++)
	{
		for(int j = x; j < x + _patchsize; j++)
		{
			if(_mask[i * _imgwidth + j] == 0)
				flag = false;
		}
	}
	return flag;
}

double* ImageCompletion::computeWHTVector(int x, int y, int PatchSize)
{
	double *wht = new double[24];
	for(int m = 0; m < 16; m++)
	{
		double SUM = 0;
		for(int i = y; i < y + PatchSize; i++)
		{
			for(int j = x; j < x + PatchSize; j++)
			{
				SUM = SUM + _imgYCbCr[_imgwidth * _imgheight * 0 + i * _imgwidth + j] * 
					_WHT[PatchSize * PatchSize * m + (i - y) * PatchSize + j - x];
			}
		}
		wht[m] = SUM;
	}
	for(int m = 0; m < 4; m++)
	{
		double SUM = 0;
		for(int i = y; i < y + PatchSize; i++)
		{
			for(int j = x; j < x + PatchSize; j++)
			{
				SUM = SUM + _imgYCbCr[_imgwidth * _imgheight * 1 + i * _imgwidth + j] * 
					_WHT[PatchSize * PatchSize * m + (i - y) * PatchSize + j - x];
			}
		}
		wht[m + 16] = SUM;
	}
	for(int m = 0; m < 4; m++)
	{
		double SUM = 0;
		for(int i = y; i < y + PatchSize; i++)
		{
			for(int j = x; j < x + PatchSize; j++)
			{
				SUM = SUM + _imgYCbCr[_imgwidth * _imgheight * 2 + i * _imgwidth + j] * 
					_WHT[PatchSize * PatchSize * m + (i - y) * PatchSize + j - x];
			}
		}
		wht[m + 20] = SUM;
	}
	return wht;
}

void ImageCompletion::generateWHTVector()
{
	findRectangle();
	createWHTBases(_patchsize);

	for(int i = _y0; i <= _y1 - _patchsize + 1; i++)
	{
		for(int j = _x0; j <= _x1 - _patchsize + 1; j++)
		{
			if(judge1(j, i, _patchsize))
			{
				_patch_count ++;
				if(_patch_count == 1)
				{
					_imgPatch = new Patch[_patch_count];
					_imgPatch[_patch_count - 1]._x = j;
					_imgPatch[_patch_count - 1]._y = i;
					_imgPatch[_patch_count - 1]._patchWHT = computeWHTVector(j, i, _patchsize);
				}
				else
				{
					Patch *temp = new Patch[_patch_count];
					for(int m = 0; m < _patch_count - 1; m++)
					{
						temp[m]._x = _imgPatch[m]._x;
						temp[m]._y = _imgPatch[m]._y;
						temp[m]._patchWHT = _imgPatch[m]._patchWHT;
					}
					temp[_patch_count - 1]._x = j;
					temp[_patch_count - 1]._y = i;
					temp[_patch_count - 1]._patchWHT = computeWHTVector(j, i, _patchsize);

					delete []_imgPatch;
					_imgPatch = temp;
					temp = NULL;
				}
			}
		}
	}
}

bool ImageCompletion::judge2(Node *node, int count)
{
	bool flag = false;
	for(int i = 0; i < count; i++)
	{
		if(node[i]._patchArray_count > 8 && (!node[i]._flag))
		{
			flag = true;
		}
	}
	return flag;
}

int ImageCompletion::calHyperPlane(Node *node, long double &yuzhi, bool &flag)
{
	double min[24];
	double max[24];
	double dmax = -10000;
	int dimension_index = -1;

	for(int i = 0; i < 24; i++)
	{
		min[i] = 10000;
		max[i] = -10000;
	}

	for(int i = 0; i < node->_patchArray_count; i++)
	{
		for(int j = 0; j < 24; j++)
		{
			if(node->_patchArray[i]._patchWHT[j] > max[j])
				max[j] = node->_patchArray[i]._patchWHT[j];
			if(node->_patchArray[i]._patchWHT[j] < min[j])
				min[j] = node->_patchArray[i]._patchWHT[j];
		}
	}
	for(int i = 0; i < 24; i++)
	{
		if(dmax < max[i] - min[i])
		{
			dmax = max[i] - min[i];
			dimension_index = i;
		}
	}

	double *a = new double[node->_patchArray_count];
	for(int i = 0; i < node->_patchArray_count; i++)
	{
		a[i] = node->_patchArray[i]._patchWHT[dimension_index];
	}
	for(int i = 0; i < (node->_patchArray_count / 2) + 1; i++)
	{
		for(int j = i + 1; j < node->_patchArray_count; j++)
		{
			if(a[i] < a[j])
			{
				double temp;
				temp = a[j];
				a[j] = a[i];
				a[i] = temp;
			}
		}
	}
	if(node->_patchArray_count % 2 == 0)
	{
		yuzhi = (a[node->_patchArray_count / 2 - 1] + a[node->_patchArray_count / 2]) / 2;
	}
	else
	{
		yuzhi = a[(node->_patchArray_count - 1) / 2];
	}

	if(max[dimension_index] == yuzhi || min[dimension_index] == yuzhi)
		flag = true;

	delete []a;
	return dimension_index;
}

void ImageCompletion::generateKDTree()
{
	_nodes_count = 1;
	_nodeArray = new Node[_nodes_count];
	_nodeArray[_nodes_count - 1]._patchArray_count = _patch_count;
	_nodeArray[_nodes_count - 1]._patchArray = new Patch[_patch_count];
	for(int i = 0; i < _patch_count; i++)
	{
		_nodeArray[_nodes_count - 1]._patchArray[i] = _imgPatch[i];
	}
	vector<Node> tempnode;
	while(judge2(_nodeArray, _nodes_count))
	{
		if(tempnode.size() != 0)
			tempnode.clear();
		for(int i = 0; i < _nodes_count; i++)
		{
			Node node = _nodeArray[i];
			if(node._patchArray_count > 8)
			{
				Node node1, node2;
				long double yuzhi;
				int index = -1;
				index = calHyperPlane(&node, yuzhi, node._flag);
				if(!node._flag)
				{
					int a = tempnode.size();

					for(int j = 0; j < node._patchArray_count; j++)
					{
						if(node._patchArray[j]._patchWHT[index] <= yuzhi)
						{
							node._patchArray[j]._node_index = a;
							node1._patchArray_count ++;
							if(node1._patchArray_count == 1)
							{
								node1._patchArray = new Patch[node1._patchArray_count];
								node1._patchArray[0] = node._patchArray[j];
							}
							else
							{
								Patch *temppatch = new Patch[node1._patchArray_count];
								for(int m = 0; m < node1._patchArray_count - 1; m++)
								{
									temppatch[m] = node1._patchArray[m];
								}
								temppatch[node1._patchArray_count - 1] = node._patchArray[j];
								delete []node1._patchArray;
								node1._patchArray = temppatch;
								temppatch = NULL;
							}
						}
						else
						{
							node._patchArray[j]._node_index = a + 1;
							node2._patchArray_count ++;
							if(node2._patchArray_count == 1)
							{
								node2._patchArray = new Patch[node2._patchArray_count];
								node2._patchArray[0] = node._patchArray[j];
							}
							else
							{
								Patch *temppatch = new Patch[node2._patchArray_count];
								for(int m = 0; m < node2._patchArray_count - 1; m++)
								{
									temppatch[m] = node2._patchArray[m];
								}
								temppatch[node2._patchArray_count - 1] = node._patchArray[j];
								delete []node2._patchArray;
								node2._patchArray = temppatch;
							}
						}
					}
					tempnode.push_back(node1);
					tempnode.push_back(node2);
				}

				else
				{
					int a = tempnode.size();
					for(int j = 0; j < node._patchArray_count; j++)
					{
						node._patchArray[j]._node_index = a;
					}
					tempnode.push_back(node);
				}
			}
			else if(node._patchArray_count <= 8 && node._patchArray_count > 0)
			{
				int a = tempnode.size();
				for(int j = 0; j < node._patchArray_count; j++)
				{
					node._patchArray[j]._node_index = a;
				}
				tempnode.push_back(node);
			}
		}
		_nodes_count = tempnode.size();
		if(_nodeArray != NULL)
			delete []_nodeArray;
		_nodeArray = new Node[_nodes_count];
		for(int i = 0; i < _nodes_count; i++)
		{
			_nodeArray[i] = tempnode[i];
		}
	}
	if(_index_frompixeltopatch != NULL)
	{
		delete []_index_frompixeltopatch;
	}
	_index_frompixeltopatch = new int[_imgheight * _imgwidth * 2];
	for(int i = 0; i < _imgheight; i++)
	{
		for(int j = 0; j < _imgwidth; j++)
		{
			_index_frompixeltopatch[i * _imgwidth + j] = -1;
			_index_frompixeltopatch[_imgheight * _imgwidth + i * _imgwidth + j] = -1;
		}
	}
	for(int i = 0; i < _nodes_count; i++)
	{
		for(int j = 0; j < _nodeArray[i]._patchArray_count; j++)
		{
			int x = _nodeArray[i]._patchArray[j]._x;
			int y = _nodeArray[i]._patchArray[j]._y;
			_index_frompixeltopatch[y * _imgwidth + x] = i;
			_index_frompixeltopatch[_imgheight * _imgwidth + y * _imgwidth + x] = j;
		}
	}
}

double ImageCompletion::calL2Distance1(Patch patch1, Patch patch2)
{
	double *a1 = patch1._patchWHT;
	double *a2 = patch2._patchWHT;
	double Sum = 0;
	for(int i = 0; i < 24; i++)
	{
		Sum += (a1[i] - a2[i]) * (a1[i] - a2[i]);
	}
	return Sum;
}

double ImageCompletion::calL2Distance2(Patch patch1, Patch patch2)
{
	int x1 = patch1._x;
	int y1 = patch1._y;
	int x2 = patch2._x;
	int y2 = patch2._y;
	double Sum = 0;
	double *data1 = NULL;
	double *data2 = NULL;
	double *data3 = NULL;
	double *data4 = NULL;
	double *data5 = NULL;
	double *data6 = NULL;
	for(int i = 0; i < _patchsize; i++)
	{
		data1 = &_imgYCbCr[(i + y1) * _imgwidth];
		data2 = &_imgYCbCr[(i + y2) * _imgwidth];
		data3 = &_imgYCbCr[(_imgwidth * _imgheight) + (i + y1) * _imgwidth];
		data4 = &_imgYCbCr[(_imgwidth * _imgheight) + (i + y2) * _imgwidth];
		data5 = &_imgYCbCr[(_imgwidth * _imgheight) * 2 + (i + y1) * _imgwidth];
		data6 = &_imgYCbCr[(_imgwidth * _imgheight) * 2 + (i + y1) * _imgwidth];
		for(int j = 0; j < _patchsize; j++)
		{
			Sum += (data1[(j + x1)] - data2[(j + x2)]) * (data1[(j + x1)] - data2[(j + x2)]) + 
				(data3[(j + x1)] - data4[(j + x2)]) * (data3[(j + x1)] - data4[(j + x2)]) + 
				(data5[(j + x1)] - data6[(j + x2)]) * (data5[(j + x1)] - data6[(j + x2)]);
		}
	}
	return Sum;
}

int ImageCompletion::propagationProcess(int x, int y)
{
	double minvalue = 10000000;
	int a = _index_frompixeltopatch[y * _imgwidth + x];
	int b = _index_frompixeltopatch[_imgwidth * _imgheight + y * _imgwidth + x];
	Patch patch1 = _nodeArray[a]._patchArray[b];
	Patch patch2;
	int index_node = -1;

	if(x - 1 >= 0 && _index_frompixeltopatch[y * _imgwidth + x - 1] != -1 && _flagmatched[y * _imgwidth + x - 1] == true)
	{
		int deltx = _deltx[y * _imgwidth + x - 1];
		int delty = _delty[y * _imgwidth + x - 1];
		int c = _index_frompixeltopatch[(y + delty) * _imgwidth + x - 1 + deltx];
		int d = _index_frompixeltopatch[_imgwidth * _imgheight + (y + delty) * _imgwidth + x - 1 + deltx];
		patch2 = _nodeArray[c]._patchArray[d];
		if(calL2Distance1(patch1, patch2) < minvalue)
		{
			minvalue = calL2Distance1(patch1, patch2);
			index_node = c;
		}
	}

	if(x - 1 >= 0 && _index_frompixeltopatch[y * _imgwidth + x - 1] != -1 && _flagmatched[y * _imgwidth + x - 1] == true)
	{
		int deltx = _deltx2[y * _imgwidth + x - 1];
		int delty = _delty2[y * _imgwidth + x - 1];
		int c = _index_frompixeltopatch[(y + delty) * _imgwidth + x - 1 + deltx];
		int d = _index_frompixeltopatch[_imgwidth * _imgheight + (y + delty) * _imgwidth + x - 1 + deltx];
		patch2 = _nodeArray[c]._patchArray[d];
		if(calL2Distance1(patch1, patch2) < minvalue)
		{
			minvalue = calL2Distance1(patch1, patch2);
			index_node = c;
		}
	}

	if(y - 1 >= 0 && _index_frompixeltopatch[(y - 1) * _imgwidth + x] != -1 && _flagmatched[(y - 1) * _imgwidth + x] == true)
	{
		int deltx = _deltx[(y - 1) * _imgwidth + x];
		int delty = _delty[(y - 1) * _imgwidth + x];
		int c = _index_frompixeltopatch[(y - 1 + delty) * _imgwidth + x + deltx];
		int d = _index_frompixeltopatch[_imgwidth * _imgheight + (y - 1 + delty) * _imgwidth + x + deltx];
		patch2 = _nodeArray[c]._patchArray[d];
		if(calL2Distance1(patch1, patch2) < minvalue)
		{
			minvalue = calL2Distance1(patch1, patch2);
			index_node = c;
		}
	}

	if(y - 1 >= 0 && _index_frompixeltopatch[(y - 1) * _imgwidth + x] != -1 && _flagmatched[(y - 1) * _imgwidth + x] == true)
	{
		int deltx = _deltx2[(y - 1) * _imgwidth + x];
		int delty = _delty2[(y - 1) * _imgwidth + x];
		int c = _index_frompixeltopatch[(y - 1 + delty) * _imgwidth + x + deltx];
		int d = _index_frompixeltopatch[_imgwidth * _imgheight + (y - 1 + delty) * _imgwidth + x + deltx];
		patch2 = _nodeArray[c]._patchArray[d];
		if(calL2Distance1(patch1, patch2) < minvalue)
		{
			minvalue = calL2Distance1(patch1, patch2);
			index_node = c;
		}
	}

	return index_node;
}

void ImageCompletion::statisticOffset(int &min1, int &max1, int &min2, int &max2)
{
	int min = 10000000;
	int max = -10000000;
	for(int i = _y0; i < _y1 - _patchsize + 1; i++)
	{
		for(int j = _x0; j < _x1 - _patchsize + 1; j++)
		{
			if(_flagmatched[i * _imgwidth + j] == true)
			{
				if(_deltx[i * _imgwidth + j] < min)
					min = _deltx[i * _imgwidth + j];
				if(_deltx[i * _imgwidth + j] > max)
					max = _deltx[i * _imgwidth + j];
			}
		}
	}
	min1 = min;
	max1 = max;
	min = 10000000;
	max = -10000000;
	for(int i = _y0; i < _y1 - _patchsize + 1; i++)
	{
		for(int j = _x0; j < _x1 - _patchsize + 1; j++)
		{
			if(_flagmatched[i * _imgwidth + j] == true)
			{
				if(_delty[i * _imgwidth + j] < min)
					min = _delty[i * _imgwidth + j];
				if(_delty[i * _imgwidth + j] > max)
					max = _delty[i * _imgwidth + j];
			}
		}
	}
	min2 = min;
	max2 = max;
}

void ImageCompletion::gaussianFilter()
{
	int Gaussian[7][7] = {{1,4,7,10,7,4,1},
	{4,12,26,33,26,12,4},
	{7,26,55,71,55,26,7},
	{10,33,71,91,71,33,10},
	{7,26,55,71,55,26,7},
	{4,12,26,33,26,12,4},
	{1,4,7,10,7,4,1}};

	int Sum = 0;
	for(int i = 3; i < _sizey - 3; i++)
	{
		for(int j = 3; j < _sizex - 3; j++)
		{
			Sum = 0;
			int p = 0;
			for(int m = i - 3; m < i + 4; m++)
			{
				int q = 0;
				for(int n = j - 3; n < j + 4; n++)
				{
					Sum += _histogram[m * _sizex + n] * Gaussian[p][q];
					q++;
				}
				p++;
			}
			_histogram[i * _sizex + j] = (double)Sum / 1115;
		}
	}
}

void ImageCompletion::statisticsPatchOffset()
{
	if(_flagmatched != NULL)
	{
		delete []_flagmatched;
	}
	_flagmatched = new bool[_imgheight * _imgwidth];
	if(_deltx != NULL)
	{
		delete []_deltx;
	}
	_deltx = new int[_imgwidth * _imgheight];
	if(_delty != NULL)
	{
		delete []_delty;
	}
	_delty = new int[_imgwidth * _imgheight];

	if(_deltx2 != NULL)
	{
		delete []_deltx2;
	}
	_deltx2 = new int[_imgwidth * _imgheight];
	if(_delty2 != NULL)
	{
		delete []_delty2;
	}
	_delty2 = new int[_imgwidth * _imgheight];

	for(int i = 0; i < _imgheight * _imgwidth; i++)
	{
		_flagmatched[i] = false;
		_deltx[i] = 0;
		_delty[i] = 0;
		_deltx2[i] = 0;
		_delty2[i] = 0;
	}
	int node_index1;
	int node_index2;
	int patch_index;
	Patch patch1;
	Patch patch2;

	double minvalue = 10000000;
	double minvalue2 = 10000000;
	int x2 = 0;
	int y2 = 0;
	int x3 = 0;
	int y3 = 0;
	bool flag = false;
	bool flag2 = false;

	int t = 0;
	if(_y1 - _y0 > _x1 - _x0)
		t = _y1 - _y0 + 1;
	else
		t = _x1 - _x0 + 1;

	for(int j = _y0; j < _y1 - _patchsize + 1; j++)
	{
		for(int i = _x0; i < _x1 - _patchsize + 1; i++)
		{
			x2 = y2 = 0;
			x3 = y3 = 0;
			minvalue = 10000000;
			minvalue2 = 10000000;
			flag = false;
			flag2 = false;
			node_index1 = node_index2 = patch_index = -1;
			node_index1 = _index_frompixeltopatch[j * _imgwidth + i];
			patch_index = _index_frompixeltopatch[_imgwidth * _imgheight + j * _imgwidth + i];

			if(node_index1 != -1 && patch_index != -1 && _flagmatched[j * _imgwidth + i] == false)
			{
				patch1 = _nodeArray[node_index1]._patchArray[patch_index];
				for(int m = 0; m < _nodeArray[node_index1]._patchArray_count; m++)
				{
					if(m != patch_index)
					{
						patch2 = _nodeArray[node_index1]._patchArray[m];
						if(calL2Distance2(patch1, patch2) < minvalue)
						{
							if(sqrt((double)((patch2._x - patch1._x) * (patch2._x - patch1._x) + (patch2._y - patch1._y) * (patch2._y - patch1._y)))
								 > (double)t / 5)
							{
								if(minvalue < minvalue2)
								{
									minvalue2 = minvalue;
									x3 = x2;
									y3 = y2;
									flag2 = true;
								}
								minvalue = calL2Distance2(patch1, patch2);
								x2 = patch2._x;
								y2 = patch2._y;
								flag = true;
								continue;
							}
						}
					}
				}
				node_index2 = propagationProcess(i, j);
				if(node_index2 != -1)
				{
					for(int m = 0; m < _nodeArray[node_index2]._patchArray_count; m++)
					{
						patch2 = _nodeArray[node_index2]._patchArray[m];
						if(calL2Distance2(patch1, patch2) < minvalue)
						{
							if(sqrt((double)((patch2._x - patch1._x) * (patch2._x - patch1._x) + (patch2._y - patch1._y) * (patch2._y - patch1._y)))
								> (double)t / 5)
							{
								if(minvalue < minvalue2)
								{
									minvalue2 = minvalue;
									x3 = x2;
									y3 = y2;
									flag2 = true;
								}
								minvalue = calL2Distance2(patch1, patch2);
								x2 = patch2._x;
								y2 = patch2._y;
								flag = true;
								continue;
							}
						}
					}
				}
				if(flag)
				{
					_deltx[j * _imgwidth + i] = x2 - patch1._x;
					_delty[j * _imgwidth + i] = y2 - patch1._y;
					_flagmatched[j * _imgwidth + i] = true;
				}
				if(flag2)
				{
					_deltx2[j * _imgwidth + i] = x3 - patch1._x;
					_delty2[j * _imgwidth + i] = y3 - patch1._y;
				}
			}
		}
	}

	statisticOffset(_minoffsetx, _maxoffsetx, _minoffsety, _maxoffsety);
	if(_histogram != NULL)
	{
		delete []_histogram;
	}
	_sizex = _maxoffsetx - _minoffsetx + 1;
	_sizey = _maxoffsety - _minoffsety + 1; 
	_histogram = new double[_sizey * _sizex];
	for(int i = 0; i < _sizey; i++)
	{
		for(int j = 0; j < _sizex; j++)
		{
			_histogram[i * _sizex + j] = 0;
		}
	}
	for(int i = _y0; i < _y1 - _patchsize + 1; i++)
	{
		for(int j = _x0; j < _x1 - _patchsize + 1; j++)
		{
			if(_flagmatched[i * _imgwidth + j] == true)
			{
				int a = _deltx[i * _imgwidth + j];
				int b = _delty[i * _imgwidth + j];
				a = a - _minoffsetx;
				b = b - _minoffsety;
				_histogram[b * _sizex + a] += 1;
			}
		}
	}

	FILE *fp;
	fp = fopen("E:\\a.txt", "wt");
	if(fp != NULL)
	{
		for(int i = 0; i < _sizey; i++)
		{
			for(int j = 0; j < _sizex; j++)
			{
				fprintf(fp, "%lf   ", _histogram[i * _sizex + j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);

	gaussianFilter();

	fp = fopen("E:\\b.txt", "wt");
	if(fp != NULL)
	{
		for(int i = 0; i < _sizey; i++)
		{
			for(int j = 0; j < _sizex; j++)
			{
				fprintf(fp, "%lf   ", _histogram[i * _sizex + j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);

	fp = fopen("E:\\c.txt", "wt");
	if(fp != NULL)
	{
		for(int i = 0; i < _imgheight; i++)
		{
			for(int j = 0; j < _imgwidth; j++)
			{
				fprintf(fp, "%d   ", _deltx[i * _imgwidth + j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);

	fp = fopen("E:\\d.txt", "wt");
	if(fp != NULL)
	{
		for(int i = 0; i < _imgheight; i++)
		{
			for(int j = 0; j < _imgwidth; j++)
			{
				fprintf(fp, "%d   ", _delty[i * _imgwidth + j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);

	pickoutKCandidates();
}  

void ImageCompletion::pickoutKCandidates()
{
	searchPeaks();
	int temp1 = 0;
	int temp2 = 0;
	int temp3 = 0;
	for(int i = 0; i < _peaks_count - 1; i++)
	{
		for(int j = i + 1; j < _peaks_count; j++)
		{
			if(_peaksvalue[i] < _peaksvalue[j])
			{
				temp1 = _peaksvalue[j];
				_peaksvalue[j] = _peaksvalue[i];
				_peaksvalue[i] = temp1;

				temp2 = _peaksx[j];
				_peaksx[j] = _peaksx[i];
				_peaksx[i] = temp2;

				temp3 = _peaksy[j];
				_peaksy[j] = _peaksy[i];
				_peaksy[i] = temp3;
			}
		}
	}
	if(_peaks_count >= 60)
	{
		_candidates_count = 60;
		if(_offsetx_candidates != NULL)
		{
			delete []_offsetx_candidates;
		}
		if(_offsety_candidates != NULL)
		{
			delete []_offsety_candidates;
		}
		_offsetx_candidates = new int[_candidates_count];
		_offsety_candidates = new int[_candidates_count];


		for(int i = 0; i < _candidates_count; i++)
		{
			_offsetx_candidates[i] = _peaksx[i] + _minoffsetx;

			_offsety_candidates[i] = _peaksy[i] + _minoffsety;
		}
	}
	else
	{
		_candidates_count = _peaks_count;
		if(_offsetx_candidates != NULL)
		{
			delete []_offsetx_candidates;
		}
		if(_offsety_candidates != NULL)
		{
			delete []_offsety_candidates;
		}
		_offsetx_candidates = new int[_candidates_count];
		_offsety_candidates = new int[_candidates_count];


		for(int i = 0; i < _candidates_count; i++)
		{
			_offsetx_candidates[i] = _peaksx[i] + _minoffsetx;

			_offsety_candidates[i] = _peaksy[i] + _minoffsety;
		}
	}

	_candidates_count += 1;
	int *candidatetemp1 = new int[_candidates_count];
	int *candidatetemp2 = new int[_candidates_count];
	for(int i = 0; i < _candidates_count - 1; i++)
	{
		candidatetemp1[i] = _offsetx_candidates[i];
		candidatetemp2[i] = _offsety_candidates[i];
	}
	candidatetemp1[_candidates_count - 1] = 0;
	candidatetemp2[_candidates_count - 1] = 0;

	delete []_offsetx_candidates;
	delete []_offsety_candidates;
	_offsetx_candidates = candidatetemp1;
	_offsety_candidates = candidatetemp2;
	candidatetemp1 = NULL;
	candidatetemp2 = NULL;
}

void ImageCompletion::searchPeaks()
{
	int index1 = 0;
	int index2 = 0;
	double tempvalue = 0;
	for(int i = 4; i < _sizey - 4; i = i + 9)
	{
		for(int j = 4; j < _sizex - 4; j = j + 9)
		{
			index1 = 0;
			index2 = 0;
			tempvalue = -1;
			for(int m = i - 4; m < i + 5; m++)
			{
				for(int n = j - 4; n < j + 5; n++)
				{
					if(tempvalue < _histogram[m * _sizex + n])
					{
						tempvalue = _histogram[m * _sizex + n];
						index1 = n;
						index2 = m;
					}
				}
			}
			if(tempvalue > 0)
			{
				_peaks_count ++;
				if(_peaks_count == 1)
				{
					_peaksvalue = new double[_peaks_count];
					_peaksx = new int[_peaks_count];
					_peaksy = new int[_peaks_count];
					_peaksvalue[_peaks_count - 1] = tempvalue;
					_peaksx[_peaks_count - 1] = index1;
					_peaksy[_peaks_count - 1] = index2;
				}
				else
				{
					double *temp1 = new double[_peaks_count];
					int *temp2 = new int[_peaks_count];
					int *temp3 = new int[_peaks_count];
					for(int m = 0; m < _peaks_count - 1; m++)
					{
						temp1[m] = _peaksvalue[m];
						temp2[m] = _peaksx[m];
						temp3[m] = _peaksy[m];
					}
					temp1[_peaks_count - 1] = tempvalue;
					temp2[_peaks_count - 1] = index1;
					temp3[_peaks_count - 1] = index2;
					delete []_peaksvalue;
					delete []_peaksx;
					delete []_peaksy;
					_peaksvalue = temp1;
					_peaksx = temp2;
					_peaksy = temp3;
					temp1 = NULL;
					temp2 = NULL;
					temp3 = NULL;
				}
			}
		}
	}
}

int ImageCompletion::findVNeighbor(int index, int x, int y)
{
	int a = -1;
	for(int i = index; i < _ukpixels_count; i++)
	{
		if(_upixels[i].x == x && _upixels[i].y == y + 1)
		{
			a = i;
			break;
		}
	}
	return a;
}

void ImageCompletion::gcOptimization()
{
	if(_optimization_result != NULL)
	{
		delete []_optimization_result;
	}

	_optimization_result = new int[_ukpixels_count];
	int num_pixels = _ukpixels_count;
	int num_labels = _candidates_count;
	int *data = new int[num_pixels * num_labels];
			
	for(int i = 0; i < num_pixels; i++)
	{
		for(int l = 0; l < num_labels; l++)
		{
			if(_upixels[i].flag_boundary )
			{
				if(l != num_labels - 1)
				{
					data[i * num_labels + l] = INT_MAX;
					continue;
				}
				else
				{
					data[i * num_labels + l] = 0;
					continue;
				}
			}
			else
			{
				if( l == num_labels - 1)
				{
					data[i * num_labels + l] = INT_MAX;
				}
				else
				{
					int delt1 = _offsetx_candidates[l];
					int delt2 = _offsety_candidates[l];
					int x = _upixels[i].x;
					int y = _upixels[i].y;
					if(x + delt1 >= 0 && x + delt1 < _imgwidth && y + delt2 >= 0 && y + delt2 < _imgheight)
					{
						if(_mask[(y + delt2) * _imgwidth + x + delt1] == 1)
						{
							data[i * num_labels + l] = 0;
						}
						else
						{
							data[i * num_labels + l] = INT_MAX;
						}
					}
					else
					{
						data[i * num_labels + l] = INT_MAX;
					}
				}
			}
		}
	}

	int *smooth = new int[num_labels * num_labels];

	for(int l1 = 0; l1 < num_labels; l1++)
	{
		for (int l2 = 0; l2 < num_labels; l2++)
		{
			if(l1 == l2)
				smooth[l1 * num_labels + l2] = 1;
			else
				smooth[l1 * num_labels + l2] = 1;
		}
	}

	//	 try
	{
		int Success = 1;
		int a = 0;
		for(int i = 0; i < num_pixels; i++)
		{
			int l = rand() % num_labels;
			int x = _upixels[i].x;
			int y = _upixels[i].y;
			int deltx = _offsetx_candidates[l];
			int delty = _offsety_candidates[l];
			if((y + delty) >= 0 && (y + delty) < _imgheight && (x + deltx) >= 0 && (x + deltx) < _imgwidth && 
				_mask[(y + delty) * _imgwidth + (x + deltx)] == 1)
			{
				_optimization_result[i] = l;
			}
			else
			{
				if(a < 10)
				{
					a ++;
					i--;
				}
				else
				{
					_optimization_result[i] = l;
				}
			}
		}

		unsigned char* Imagedata = _imgRGB->bits();
		int64_t min_energy;
		GCoptimizationGeneralGraph *gc = NULL;

		gc = new GCoptimizationGeneralGraph(num_pixels, num_labels);

		gc->setDataCost(data);
		gc->setSmoothCost(smooth);

//		while(Success == 1)
		{
			Success = 0;

			for(int i = 0; i < num_pixels; i++)
			{
				gc->setLabel(i, _optimization_result[i]);
			}

			for(int i = 0; i < num_pixels - 1; i++)
			{
				int x = _upixels[i].x;
				int y = _upixels[i].y;
				if(_upixels[i + 1].y == y && _upixels[i + 1].x == x + 1)
				{
					int l1 = _optimization_result[i];
					int l2 = _optimization_result[i + 1];
					if(l1 == l2)
					{
						gc->setNeighbors(i, i + 1, 0);
					}
					else
					{
						int deltx1 = _offsetx_candidates[l1];
						int delty1 = _offsety_candidates[l1];
						int deltx2 = _offsetx_candidates[l2];
						int delty2 = _offsety_candidates[l2];

						int x1 = x + deltx1;
						int y1 = y + delty1;
						int x1_2 = x + deltx2;
						int y1_2 = y + delty2;

						int x2 = x + 1 + deltx1;
						int y2 = y + delty1;
						int x2_2 = x + 1 + deltx2;
						int y2_2 = y + delty2;

						int a[3], b[3], c[3], d[3], e[6], f[6], g[6], h[6];
						int Sum = 0;
						if(x1 >= 0 && x1 < _imgwidth - 1 && y1 >= 0 && y1 < _imgheight - 1 && x1_2 >= 0 && x1_2 < _imgwidth - 1 && y1_2 >= 0 && y1_2 < _imgheight - 1
							&& x2 >= 0 && x2 < _imgwidth - 1 && y2 >= 0 && y2 < _imgheight - 1 && x2_2 >= 0 && x2_2 < _imgwidth - 1 && y2_2 >= 0 && y2_2 < _imgheight - 1
							&& _mask[y1 * _imgwidth + x1] == 1 && _mask[y1_2 * _imgwidth + x1_2] == 1 && _mask[y2 * _imgwidth + x2] == 1 && _mask[y2_2 * _imgwidth + x2_2] == 1)
						{
							Imagedata = _imgRGB->scanLine(y1);
							a[0] = Imagedata[x1 * 4 + 2];
							a[1] = Imagedata[x1 * 4 + 1];
							a[2] = Imagedata[x1 * 4 + 0];

							e[0] = Imagedata[(x1 + 1) * 4 + 2] - a[0];
							e[1] = Imagedata[(x1 + 1) * 4 + 1] - a[1];
							e[2] = Imagedata[(x1 + 1) * 4 + 0] - a[2];

							Imagedata = _imgRGB->scanLine(y1 + 1);
							e[3] = Imagedata[x1 * 4 + 2] - a[0];
							e[4] = Imagedata[x1 * 4 + 1] - a[1];
							e[5] = Imagedata[x1 * 4 + 0] - a[2];

							Imagedata = _imgRGB->scanLine(y1_2);
							b[0] = Imagedata[x1_2 * 4 + 2];
							b[1] = Imagedata[x1_2 * 4 + 1];
							b[2] = Imagedata[x1_2 * 4 + 0];

							f[0] = Imagedata[(x1_2 + 1) * 4 + 2] - b[0];
							f[1] = Imagedata[(x1_2 + 1) * 4 + 1] - b[1];
							f[2] = Imagedata[(x1_2 + 1) * 4 + 0] - b[2];

							Imagedata = _imgRGB->scanLine(y1_2 + 1);
							f[3] = Imagedata[x1_2 * 4 + 2] - b[0];
							f[4] = Imagedata[x1_2 * 4 + 1] - b[1];
							f[5] = Imagedata[x1_2 * 4 + 0] - b[2];

							Imagedata = _imgRGB->scanLine(y2);
							c[0] = Imagedata[x2 * 4 + 2];
							c[1] = Imagedata[x2 * 4 + 1];
							c[2] = Imagedata[x2 * 4 + 0];

							g[0] = Imagedata[(x2 + 1) * 4 + 2] - c[0];
							g[1] = Imagedata[(x2 + 1) * 4 + 1] - c[1];
							g[2] = Imagedata[(x2 + 1) * 4 + 0] - c[2];
							   
							Imagedata = _imgRGB->scanLine(y2 + 1);
							g[3] = Imagedata[x2 * 4 + 2] - c[0];
							g[4] = Imagedata[x2 * 4 + 1] - c[1];
							g[5] = Imagedata[x2 * 4 + 0] - c[2];

							Imagedata = _imgRGB->scanLine(y2_2);
							d[0] = Imagedata[x2_2 * 4 + 2];
							d[1] = Imagedata[x2_2 * 4 + 1];
							d[2] = Imagedata[x2_2 * 4 + 0];

							h[0] = Imagedata[(x2_2 + 1) * 4 + 2] - d[0];
							h[1] = Imagedata[(x2_2 + 1) * 4 + 1] - d[1];
							h[2] = Imagedata[(x2_2 + 1) * 4 + 0] - d[2];

							Imagedata = _imgRGB->scanLine(y2_2 + 1);
							h[3] = Imagedata[x2_2 * 4 + 2] - d[0];
							h[4] = Imagedata[x2_2 * 4 + 1] - d[1];
							h[5] = Imagedata[x2_2 * 4 + 0] - d[2];

							Sum = (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) + 
								(a[2] - b[2])*(a[2] - b[2]) + (c[0] - d[0])*(c[0] - d[0]) + (c[1] - d[1])*(c[1] - d[1])
								+ (c[2] - d[2])*(c[2] - d[2]) + 2 * ((e[0] - f[0]) * (e[0] - f[0]) + (e[1] - f[1]) * (e[1] - f[1]) + 
								(e[2] - f[2]) * (e[2] - f[2]) + (g[0] - h[0]) * (g[0] - h[0]) + (g[1] - h[1]) * (g[1] - h[1]) + 
								(g[2] - h[2]) * (g[2] - h[2]) + (e[3] - f[3]) * (e[3] - f[3]) + (e[4] - f[4]) * (e[4] - f[4]) + 
								(e[5] - f[5]) * (e[5] - f[5]) + (g[3] - h[3]) * (g[3] - h[3]) + (g[4] - h[4]) * (g[4] - h[4]) + 
								(g[5] - h[5]) * (g[5] - h[5]));
							gc->setNeighbors(i, i + 1, Sum);
						}
					}
				}
			}

			for(int i = 0; i < num_pixels; i++)
			{
				int x = _upixels[i].x;
				int y = _upixels[i].y;
				int index = findVNeighbor(i, x, y);
				if(index != -1)
				{
					int l1 = _optimization_result[i];
					int l2 = _optimization_result[index];
					if(l1 == l2)
					{
						gc->setNeighbors(i, index, 0);
					}
					else
					{
						int deltx1 = _offsetx_candidates[l1];
						int delty1 = _offsety_candidates[l1];
						int deltx2 = _offsetx_candidates[l2];
						int delty2 = _offsety_candidates[l2];

						int x1 = x + deltx1;
						int y1 = y + delty1;
						int x1_2 = x + deltx2;
						int y1_2 = y + delty2;

						int x2 = x + deltx1;
						int y2 = y + 1 + delty1;
						int x2_2 = x + deltx2;
						int y2_2 = y + 1 + delty2;

						int a[3], b[3], c[3], d[3], e[6], f[6], g[6], h[6];
						int Sum = 0;
						if(x1 >= 0 && x1 < _imgwidth - 1 && y1 >= 0 && y1 < _imgheight - 1 && x1_2 >= 0 && x1_2 < _imgwidth - 1 && y1_2 >= 0 && y1_2 < _imgheight - 1
							&& x2 >= 0 && x2 < _imgwidth - 1 && y2 >= 0 && y2 < _imgheight - 1 && x2_2 >= 0 && x2_2 < _imgwidth - 1 && y2_2 >= 0 && y2_2 < _imgheight - 1
							&& _mask[y1 * _imgwidth + x1] == 1 && _mask[y1_2 * _imgwidth + x1_2] == 1 && _mask[y2 * _imgwidth + x2] == 1 && _mask[y2_2 * _imgwidth + x2_2] == 1)
						{
							Imagedata = _imgRGB->scanLine(y1);
							a[0] = Imagedata[x1 * 4 + 2];
							a[1] = Imagedata[x1 * 4 + 1];
							a[2] = Imagedata[x1 * 4 + 0];

							e[0] = Imagedata[(x1 + 1) * 4 + 2] - a[0];
							e[1] = Imagedata[(x1 + 1) * 4 + 1] - a[1];
							e[2] = Imagedata[(x1 + 1) * 4 + 0] - a[2];

							Imagedata = _imgRGB->scanLine(y1 + 1);
							e[3] = Imagedata[x1 * 4 + 2] - a[0];
							e[4] = Imagedata[x1 * 4 + 1] - a[1];
							e[5] = Imagedata[x1 * 4 + 0] - a[2];

							Imagedata = _imgRGB->scanLine(y1_2);
							b[0] = Imagedata[x1_2 * 4 + 2];
							b[1] = Imagedata[x1_2 * 4 + 1];
							b[2] = Imagedata[x1_2 * 4 + 0];

							f[0] = Imagedata[(x1_2 + 1) * 4 + 2] - b[0];
							f[1] = Imagedata[(x1_2 + 1) * 4 + 1] - b[1];
							f[2] = Imagedata[(x1_2 + 1) * 4 + 0] - b[2];

							Imagedata = _imgRGB->scanLine(y1_2 + 1);
							f[3] = Imagedata[x1_2 * 4 + 2] - b[0];
							f[4] = Imagedata[x1_2 * 4 + 1] - b[1];
							f[5] = Imagedata[x1_2 * 4 + 0] - b[2];

							Imagedata = _imgRGB->scanLine(y2);
							c[0] = Imagedata[x2 * 4 + 2];
							c[1] = Imagedata[x2 * 4 + 1];
							c[2] = Imagedata[x2 * 4 + 0];

							g[0] = Imagedata[(x2 + 1) * 4 + 2] - c[0];
							g[1] = Imagedata[(x2 + 1) * 4 + 1] - c[1];
							g[2] = Imagedata[(x2 + 1) * 4 + 0] - c[2];

							Imagedata = _imgRGB->scanLine(y2 + 1);
							g[3] = Imagedata[x2 * 4 + 2] - c[0];
							g[4] = Imagedata[x2 * 4 + 1] - c[1];
							g[5] = Imagedata[x2 * 4 + 0] - c[2];

							Imagedata = _imgRGB->scanLine(y2_2);
							d[0] = Imagedata[x2_2 * 4 + 2];
							d[1] = Imagedata[x2_2 * 4 + 1];
							d[2] = Imagedata[x2_2 * 4 + 0];

							h[0] = Imagedata[(x2_2 + 1) * 4 + 2] - d[0];
							h[1] = Imagedata[(x2_2 + 1) * 4 + 1] - d[1];
							h[2] = Imagedata[(x2_2 + 1) * 4 + 0] - d[2];

							Imagedata = _imgRGB->scanLine(y2_2 + 1);
							h[3] = Imagedata[x2_2 * 4 + 2] - d[0];
							h[4] = Imagedata[x2_2 * 4 + 1] - d[1];
							h[5] = Imagedata[x2_2 * 4 + 0] - d[2];

							Sum = (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) + 
								(a[2] - b[2])*(a[2] - b[2]) + (c[0] - d[0])*(c[0] - d[0]) + (c[1] - d[1])*(c[1] - d[1])
								+ (c[2] - d[2])*(c[2] - d[2]) + 2 * ((e[0] - f[0]) * (e[0] - f[0]) + (e[1] - f[1]) * (e[1] - f[1]) + 
								(e[2] - f[2]) * (e[2] - f[2]) + (g[0] - h[0]) * (g[0] - h[0]) + (g[1] - h[1]) * (g[1] - h[1]) + 
								(g[2] - h[2]) * (g[2] - h[2]) + (e[3] - f[3]) * (e[3] - f[3]) + (e[4] - f[4]) * (e[4] - f[4]) + 
								(e[5] - f[5]) * (e[5] - f[5]) + (g[3] - h[3]) * (g[3] - h[3]) + (g[4] - h[4]) * (g[4] - h[4]) + 
								(g[5] - h[5]) * (g[5] - h[5]));
							gc->setNeighbors(i, index, Sum);
						}
					}
				}
			}

			min_energy = gc->compute_energy();
			gc->swap(10);
			if(gc->compute_energy() < min_energy)
			{
				min_energy = gc->compute_energy();
				Success = 1;
				for(int i = 0; i < num_pixels; i++)
				{
					_optimization_result[i] = gc->whatLabel(i);
				}
			}
		}

		int win_width = _x1 - _x0 + 1;
		int win_height = _y1 - _y0 + 1;

 		int *result = new int[win_height * win_width];
		for(int i = 0; i < win_height; i++)
		{
			for(int j = 0; j < win_width; j++)
			{
				result[i * win_width + j] = -1;
			}
		}
		for(int i = 0; i < num_pixels; i++)
		{
			int x = _upixels[i].x;
			int y = _upixels[i].y;
			int l = _optimization_result[i];
			result[(y - _y0) * win_width + x - _x0] = l;
		}
		FILE *fp;
		fp = fopen("E:\\aaa.txt", "wt");
		if(fp != NULL)
		{
			for(int i = 0; i < win_height; i++)
			{
				for(int j = 0; j < win_width; j++)
				{
					fprintf(fp, "%d   ", result[i * win_width + j]);
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);
		delete []result;
	}
}

void ImageCompletion::completionImage()
{
	int R = 0;
	int G = 0;
	int B = 0;
	int offsetx = 0;
	int offsety = 0;
	unsigned char *data = _imgRGB->bits();
	for(int i = 0; i < _ukpixels_count; i++)
	{
		int l = _optimization_result[i];
		int x = _upixels[i].x;
		int y = _upixels[i].y;

		offsetx = _offsetx_candidates[l];
		offsety = _offsety_candidates[l];

		data = _imgRGB->scanLine(y + offsety);
		R = data[(x + offsetx) * 4 + 2];
		G = data[(x + offsetx) * 4 + 1];
		B = data[(x + offsetx) * 4 + 0];

		data = _imgRGB->scanLine(y);
		data[x * 4 + 2] = R;
		data[x * 4 + 1] = G;
		data[x * 4 + 0] = B;
	}
	ui.label->setPixmap(QPixmap::fromImage(*_imgRGB));
	_imgRGB->save("E:\\a.jpg");
}

void ImageCompletion::createWHTBases(int PatchSize)
{
	if(_WHT != NULL)
		delete []_WHT;
	_WHT = new int[PatchSize * PatchSize * 16];
	int i = 0;
	int j = 0;
	for(i = 0; i < PatchSize * PatchSize * 16; i++)
	{
		_WHT[i] = 1;
	}
	//base2
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 1 + i * PatchSize + j] = -1;
		}
	}
	//base3
	for(i = 0; i < PatchSize; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 2 + i * PatchSize + j] = -1;
		}
	}
	//base4
	for(i = 0; i < PatchSize / 2; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 3 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 3 + i * PatchSize + j] = -1;
		}
	}
	//base5
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 4 + i * PatchSize + j] = -1;
		}
	}
	//base6
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 5 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 5 + i * PatchSize + j] = -1;
		}
	}
	//base7
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 6 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 6 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 6 + i * PatchSize + j] = -1;
		}
	}
	//base8
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	//base9
	for(i = 0; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 8 + i * PatchSize + j] = -1;
		}
	}
	//base10
	for(i = 0; i < PatchSize / 2; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	//base11
	for(i = 0; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize/2; j++)
		{
			_WHT[PatchSize * PatchSize * 10 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 10 + i * PatchSize + j] = -1;
		}
	}
	//base12
	for(i = 0; i < PatchSize / 2; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
	}
	//base13
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 12 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 12 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 12 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 12 + i * PatchSize + j] = -1;
		}
	}
	//base14
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	//base15
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
		for(j =PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	//base 16
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}


/*



	//base 2
	for(i = 0; i < PatchSize; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 1 + i * PatchSize + j] = -1;
		}
	}
	//base 3
	for(i = 0; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 2 + i * PatchSize + j] = -1;
		}
	}
	//base 4
	for(i = 0; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize/2; j++)
		{
			_WHT[PatchSize * PatchSize * 3 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 3 + i * PatchSize + j] = -1;
		}
	}
	//base 5
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 4 + i * PatchSize + j] = -1;
		}
	}
	//base 6
	for(i = 0; i < PatchSize / 2; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 5 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 5 + i * PatchSize + j] = -1;
		}
	}
	//base 7
	for(i = 0; i < PatchSize / 2; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 6 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 6 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 6 + i * PatchSize + j] = -1;
		}
	}
	//base 8
	for(i = 0; i < PatchSize / 2; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 7 + i * PatchSize + j] = -1;
		}
	}
	//base 9
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 8 + i * PatchSize + j] = -1;
		}
	}
	//base 10
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 9 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 9 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 9 + i * PatchSize + j] = -1;
		}
	}
	//base 11
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 10 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 10 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 10 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 10 + i * PatchSize + j] = -1;
		}
	}
	//base 12
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
		for(j =PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 11 + i * PatchSize + j] = -1;
		}
	}
	//base 13
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 12 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 12 + i * PatchSize + j] = -1;
		}
	}
	//base 14
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = PatchSize / 2; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 13 + i * PatchSize + j] = -1;
		}
	}
	//base 15
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 14 + i * PatchSize + j] = -1;
		}
	}
	//base 16
	for(i = 0; i < PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 4; i < PatchSize / 2; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize / 2; i < PatchSize - PatchSize / 4; i++)
	{
		for(j = PatchSize / 4; j < PatchSize / 2; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize - PatchSize / 4; j < PatchSize; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	for(i = PatchSize - PatchSize / 4; i < PatchSize; i++)
	{
		for(j = 0; j < PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
		for(j = PatchSize / 2; j < PatchSize - PatchSize / 4; j++)
		{
			_WHT[PatchSize * PatchSize * 15 + i * PatchSize + j] = -1;
		}
	}
	*/
}