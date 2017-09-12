
#include "GCoptimization.h"
#include <iostream>
#include <string>
#include <fstream>

struct offset_num
{
	int offset_x;
	int offset_y;
	int num;
};

struct grid_mask
{
	int first_row;
	int last_row;
	int *left_x;//对应每一行mask像素的最小x值;
	int *right_x;//对应每一行mask像素的最大x值;
};

struct coordinate
{
	int x;
	int y;
};//坐标值;

struct grid_extra_data
{
	unsigned int***pimage_red;
	unsigned int***pimage_green;
	unsigned int***pimage_blue;
	struct grid_mask**pmask;
	struct grid_mask**pbroaden_mask;
	struct offset_num***poffsets;
	int row;
	int col;
	int rect_width;
	int rect_height;
	int rect_left;
	int rect_up;
	int num_labels;//offset_数量;
	double **pimage_grid_x;
	double **pimage_grid_y;
};


void grid_label_optimize(unsigned int**image_red, unsigned int**image_green, unsigned int**image_blue, 
    struct grid_mask*mask,struct grid_mask*broaden_mask,
    struct offset_num**offset_num, int k, int **pNewResult,int **pOldResult,
	int rect_first_row, int rect_last_row, int rect_left, int rect_right, int *edge, int edge_num,
	int image_col, int image_row,double *image_grid_x,double *image_grid_y);
void outputdata(std::string file_name, int *result, int k, struct offset_num**offsets, int edge_num, int pixel_num);
struct offset_num**offset_init(int *data_num, int *data_x, int *data_y, int col, int row);
struct offset_num**sort(struct offset_num**offset_num, int k, int n);
int dichotomy_find_place(struct offset_num*offset, struct offset_num**k_offset, int first_one, int last_one);
void image_init(unsigned int*image, unsigned int***new_image, int row, int col);
void inputdata(std::string file_name, unsigned int**image_8_red, unsigned int**image_8_green, unsigned int**image_8_blue,
	int *row_bin_num, int *col_bin_num, int **bin_num, int **bin_x, int **bin_y, int *image_row, int *image_col, int **mask_matrix,
	int **broaden_mask_matrix,int *edge_num, int **edge_data, int *first_row, int *last_row, int *left, int *right,
	double **grid_image_x,double **grid_image_y);
void deletedata(unsigned int**image_8_red, unsigned int**image_8_green, unsigned int**image_8_blue,
	int *row_bin_num, int *col_bin_num, int **bin_num, int **bin_x, int **bin_y, int *image_row, int *image_col, int **mask_matrix,
	int *edge_num, int **edge_data);
void outputdata(std::string file_name, int *NewResult, int *OldResult, int k, struct offset_num**offset, int edge_num, int pixel_num);
bool grid_pixel_valid(int x, int y, struct grid_mask*true_maskint, int image_col, int image_row);
struct coordinate grid_get_coor(int no, int mask_width, int mask_height, int left_col, int first_row);
void grid_mask_init(int *mask_matrix, struct grid_mask**pmask, int row);
int GirdSmoothCostFn(int siteID1, int siteID2, int labelID1, int labelID2, void *extra_data);
int NewDataCost(int siteID, int labelID, void *extra_data);
int OldGirdSmoothCostFn(int siteID1, int siteID2, int labelID1, int labelID2, void *extra_data);
int mask_num(struct grid_mask*mask);//返回该mask中mask像素个数;