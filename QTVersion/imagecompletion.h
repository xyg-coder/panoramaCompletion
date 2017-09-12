#ifndef IMAGECOMPLETION_H
#define IMAGECOMPLETION_H

#include <QtGui/QMainWindow>
#include "ui_imagecompletion.h"

#include <omp.h>

#include <QAction>
#include <QMessageBox>
#include <QImage>
#include <QFileDialog>
#include <math.h>
#include <stdio.h>
#include "Patch.h"
#include "Node.h"
#include "GCoptimization.h"
#include <stdint.h>

#include <iostream>
using namespace std;

struct UKPixels
{
	int x;
	int y;
	bool flag_boundary;
};

class ImageCompletion : public QMainWindow
{
	Q_OBJECT

public:
	ImageCompletion(QWidget *parent = 0, Qt::WFlags flags = 0);
	~ImageCompletion();

public:
	QImage *_imgRGB;       //输入的原始待修补彩色图像
	double *_imgYCbCr;     //转到YCbCr颜色空间后的灰度数据
	int *_WHT;             //WHT变换的16个变换基
	QImage *_imgMask;      //定位图像缺损部位的8位掩膜图像
	int *_mask;            //掩膜的二进制表示 1为已知区域 0为未知区域
	int _imgwidth;         //图像的宽度
	int _imgheight;        //图像的高度 原始图像与掩膜图像尺寸相同
	int _ukpixels_count;   //图像中待修补的像素个数
	int _patchsize;        //patch的尺寸
	Patch *_imgPatch;      //指向在统计窗口内所有的Patch
	int _patch_count;      //patch的总数目

	int _x0;               //搜索窗口的左上角x坐标
	int _y0;               //搜索窗口的左上角y坐标
	int _x1;               //搜索窗口的右下角x坐标
	int _y1;               //搜索窗口的右下角y坐标

	UKPixels *_upixels;    //所有待修补像素的图像坐标
	Node *_nodeArray;      // KD-Tree下所有叶子节点对应的指针
	int _nodes_count;      //KD-Tree下叶子节点的数目
	int *_index_frompixeltopatch;
	bool *_flagmatched;   //判断该位置的patch是否完成了match
	int *_deltx;
	int *_delty;
	int *_deltx2;
	int *_delty2;

	int _minoffsetx;
	int _minoffsety;
	int _maxoffsetx;
	int _maxoffsety;

	double *_histogram;
	int _sizex;
	int _sizey;

	double *_peaksvalue;
	int *_peaksx;
	int *_peaksy;
	int _peaks_count;

	int *_offsetx_candidates;
	int *_offsety_candidates;
	int _candidates_count;

	int *_optimization_result;

public:
	void createWHTBases(int PatchSize);
	void findRectangle();
	bool judge1(int x, int y, int size); //判断以(x,y)为左上角的size x size的
										 //patch中是否包含待修补的像素
	bool judge2(Node *node, int count);  //生成K-D树时  用来判断一个节点是否需要进行分解
	double* computeWHTVector(int x, int y, int PatchSize);
	int calHyperPlane(Node *node, long double &yuzhi, bool &flag);
	double calL2Distance1(Patch patch1, Patch patch2);
	int propagationProcess(int x, int y);
	double calL2Distance2(Patch patch1, Patch patch2);
	void statisticOffset(int &min1, int &max1, int &min2, int &max2);
	void gaussianFilter();
	void pickoutKCandidates();
	void searchPeaks();
	int findVNeighbor(int index, int x, int y);

private slots:
	void openRGBImage();
	void openMaskImage();
	void toYCbCrChannel();
	void generateWHTVector();
	void generateKDTree();
	void statisticsPatchOffset();
	void gcOptimization();
	void completionImage();

private:
	Ui::ImageCompletionClass ui;
};

#endif // IMAGECOMPLETION_H
