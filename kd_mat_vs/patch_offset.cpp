#include "mex.h"
#include "kd_tree.h"
#define Image_feature 0
#define Image 1
#define Mask 2

/*
输入参数;
1. patch的feature数组;(feature传入是不用倒置)
2. image的feature数组;
3. patch;
4. image;(pixel数组传入时需要倒置);
*/
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{
	unsigned char *image_8 = NULL;
	unsigned int **image = NULL;
	image_8 = (unsigned char*)mxGetData(in[Image]);
	int *mask_matrix = (int *)mxGetData(in[Mask]);
	int mask_col = mxGetM(in[Mask]);
	struct feature*image_feature = NULL;
	int n = mxGetM(in[Image_feature]);//特征向量维数;
	int nImage_feature = mxGetN(in[Image_feature]);
	int image_row = mxGetN(in[Image]);
	int image_col = mxGetM(in[Image]);

	//feature_init(&patch_feature, in[Patch_feature], n, nPatch_feature, patch_col - 7,patch_row-7);
	feature_init(&image_feature, in[Image_feature], n, nImage_feature, image_col - 7, image_row - 7);
	//struct feature*f = image_feature + (image_col - 7) * 9 + 64;
	//image_init(patch_image_8, &patch_image, patch_row, patch_col);
	image_init(image_8, &image, image_row, image_col);
	struct mask* mask = (struct mask*)malloc(sizeof(struct mask));
	mask_init(mask_matrix, &mask, image_row, mask_col);
	int n_feature=patch_match(image_feature, image, image_col, image_row, nImage_feature, mask,out);
	
	//以下为输出部分;
	for (int i = 0; i < image_row; i++)
	{
		free(image[i]);
	}
	free(image);
	for (int i = 0; i <n_feature; i++)
	{
		free((image_feature + i)->fwd_match);
	}
	free(image_feature);
	mask_release(&mask, image_row);
}

/*
feat为传进来的feature指针（未构造）;
feat_array为matlab中传进来的数组;
dim为特征子的维数;
n为特征子的总数;
mode为初始化的模式,1为patch,2为image;*/
void feature_init(struct feature**feat, const mxArray*feat_array, int dim, int n, int col, int row)
{
	int *descr_array = (int *)mxGetData(feat_array);
	*feat = (struct feature*)malloc(n*sizeof(struct feature));
	for (int i = 0; i < n; i++)
	{
		struct feature*tmp = *feat + i;
		tmp->y = i / col;
		tmp->x = i%col;
		tmp->d = dim;
		tmp->feature_data = NULL;
		tmp->descr = descr_array + (i*dim);
		tmp->valid = true;
		tmp->p_m_leaf = tmp->p_leaf = NULL;
	}
}

/* 将传进来的unsigned char数据转化为unsigned int二维数组;*/
void image_init(unsigned char*image, unsigned int***new_image, int row, int col)
{
	*new_image = (unsigned int**)malloc(row*sizeof(unsigned int*));
	for (int i = 0; i < row; i++)
	{
		(*new_image)[i] = (unsigned int*)malloc(col*sizeof(unsigned int));
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			(*new_image)[i][j] = (*(image + col*i + j));
		}
	}
}

void mask_init(int *mask_matrix, struct mask**mask, int row, int col)
{
	(*mask)->row = row;
	(*mask)->lines = (struct line_mask**)malloc(row*sizeof(struct line_mask*));
	(*mask)->n = (int *)malloc(row*sizeof(int));
	(*mask)->min = (int *)malloc(row*sizeof(int));
	(*mask)->max = (int *)malloc(row*sizeof(int));
	for (int i = 0; i < row; i++)
	{
		int num = (*(mask_matrix + i*col));
		*((*mask)->n + i) = num;
		*((*mask)->lines + i) = (struct line_mask*)malloc(num*sizeof(struct line_mask));
		for (int j = 0; j < num; j++)
		{
			int left = (*(mask_matrix + i*col + j * 2 + 1)) - 1;
			int right = (*(mask_matrix + i*col + j * 2 + 2)) - 1;
			(*((*((*mask)->lines + i)) + j)).left = left;
			(*((*((*mask)->lines + i)) + j)).right = right;
			if (j == 0)
			{
				*((*mask)->min + i) = left;
			}
			if (j == num - 1)
			{
				*((*mask)->max + i) = right;
			}
		}
	}
}

void mask_release(struct mask**mask, int row)
{
	free((*mask)->min);
	free((*mask)->max);
	for (int i = 0; i < row; i++)
	{
		free(*((*mask)->lines + i));
	}
	free((*mask)->n);
	free((*mask)->lines);
	free((*mask));
}
