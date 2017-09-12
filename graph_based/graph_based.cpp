#define MaxDataCost 1000000
#define MaxSmoothCost 10000
#include "graph_based.h"
#define NOT_EXIST -1
#define PI 3.1415926535898

int main()
{
	unsigned int*image_8_red, *image_8_green, *image_8_blue;
	int *edge_data, *mask_matrix, image_row, image_col;
	int *NewResult, *OldResult;
	int rect_first_row, rect_last_row, rect_left, rect_right;
	int row_bin_num, col_bin_num, *offset_data_num, *offset_x, *offset_y;
	int edge_num;
	double *image_grid_x, *image_grid_y;//储存图像的梯度;
	int *broaden_mask_matrix;
	std::string input_path_name = "E:\\c++\\graph_based grid\\graph_based grid\\data.txt";
	std::string output_path_name = "E:\\c++\\graph_based grid\\graph_based grid\\result.txt";
	inputdata(input_path_name,&image_8_red,&image_8_green,&image_8_blue,&row_bin_num,&col_bin_num,&offset_data_num,&offset_x,
		&offset_y,&image_row,&image_col,&mask_matrix,&broaden_mask_matrix,&edge_num,&edge_data,&rect_first_row,&rect_last_row,&rect_left,&rect_right,
		&image_grid_x,&image_grid_y);
	unsigned int**image_red = NULL;
	unsigned int**image_green = NULL;
	unsigned int**image_blue = NULL;
	image_init(image_8_red, &image_red, image_row, image_col);
	image_init(image_8_green, &image_green, image_row, image_col);
	image_init(image_8_blue, &image_blue, image_row, image_col);
	struct grid_mask* mask = (struct grid_mask*)malloc(sizeof(struct grid_mask));
	struct grid_mask* broaden_mask = (struct grid_mask*)malloc(sizeof(struct grid_mask));
	grid_mask_init(mask_matrix, &mask,image_row);
	grid_mask_init(broaden_mask_matrix, &broaden_mask, image_row);
	struct offset_num**offsets = offset_init(offset_data_num, offset_x, offset_y, col_bin_num, row_bin_num);
	int k = 60;
	int offsets_num = row_bin_num*col_bin_num;
	struct offset_num**k_offsets = sort(offsets, k, offsets_num);
	int new_k = k;
	for (int i = 0; i < k;i++)
	{
		struct offset_num*offset = (*(k_offsets + i));
		if (offset->num<20)
		{
			new_k = i;
			break;
		}
	}
	//new_k = k;
	(*(k_offsets + new_k)) = (*(k_offsets + k));
	grid_label_optimize(image_red, image_green, image_blue, mask, broaden_mask,
		k_offsets, new_k, &NewResult, &OldResult,rect_first_row, rect_last_row, rect_left,
		rect_right, edge_data, edge_num, image_col, image_row,image_grid_x,image_grid_y);
	int num_pixels = (rect_last_row - rect_first_row + 1)*(rect_right - rect_left + 1);
	outputdata(output_path_name, NewResult,OldResult, k, k_offsets, edge_num, num_pixels);
	
	deletedata(&image_8_red, &image_8_green, &image_8_blue, &row_bin_num, &col_bin_num, &offset_data_num,
		&offset_x, &offset_y, &image_row, &image_col, &mask_matrix, &edge_num, &edge_data);
	free(mask->left_x);
	free(mask->right_x);
	free(mask);
	free(NewResult);
	free(OldResult);
	for (int i = 0; i < image_row; i++)
	{
		free(image_red[i]);
	}
	free(image_red);
	for (int i = 0; i < image_row; i++)
	{
		free(image_green[i]);
	}
	free(image_green);
	for (int i = 0; i < image_row; i++)
	{
		free(image_blue[i]);
	}
	free(image_blue);
	for (int i = 0; i < offsets_num; i++)
	{
		free((*(offsets + i)));
	}
	free(offsets);
	free(k_offsets);
	free(image_grid_x);
	free(image_grid_y);
}

void outputdata(std::string file_name, int *NewResult, int *OldResult,int k, struct offset_num**offsets, int edge_num, int pixel_num)
{
	std::ofstream out(file_name);
	out << k << "\n";
	for (int i = 0; i < k; i++)
	{
		out << (*(offsets + i))->offset_x << " " << (*(offsets + i))->offset_y << "\n";
	}
	out << pixel_num << "\n";
	for (int i = 0; i < pixel_num; i++)
	{
		out << OldResult[i] << "\n";
	}
	for (int i = 0; i < pixel_num;i++)
	{
		out << NewResult[i] << "\n";
	}
	out.close();
}

void image_init(unsigned int*image, unsigned int***new_image, int row, int col)
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

struct offset_num**offset_init(int *data_num, int *data_x, int *data_y, int col, int row)
{
	struct offset_num**offsets = (struct offset_num**)malloc(col*row*sizeof(struct offset*));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			struct offset_num*offset = (struct offset_num*)malloc(sizeof(offset_num));
			offset->num = (*(data_num + j*row + i));
			offset->offset_x = (*(data_x + j*row + i));
			offset->offset_y = (*(data_y + j*row + i));
			*(offsets + i*col + j) = offset;
		}
	}
	return offsets;
}

bool grid_pixel_valid(int x, int y, struct grid_mask*true_mask, int image_col, int image_row)
{
	if (x < 0 || x >= image_col)
	{
		return false;
	}
	if (y < 0 || y >= image_row)
	{
		return false;
	}
	if (y < true_mask->first_row)
	{
		return true;
	}
	if (y > true_mask->last_row)
	{
		return true;
	}
	int left = (*(true_mask->left_x + y));
	int right = (*(true_mask->right_x + y));
	if (x >= left&&x <= right)
	{
		return false;
	}
	else
	{
		return true;
	}
}

struct coordinate grid_get_coor(int no, int mask_width, int mask_height, int left_col, int first_row)
{
	struct coordinate coor;
	coor.y = no / mask_width;
	coor.x = no - (coor.y*mask_width);
	coor.y = coor.y + first_row;
	coor.x = coor.x + left_col;
	return coor;
}

void grid_mask_init(int *mask_matrix, struct grid_mask**pmask, int row)
{
	(*pmask)->left_x = (int *)malloc(row*sizeof(int));
	(*pmask)->right_x = (int *)malloc(row*sizeof(int));
	bool come_up = false;//如果第一个还没有出现则是false;
	bool the_last_one = false;
	for (int i = 0; i < row; i++)
	{
		int n = mask_matrix[i * 3];
		if (n == 1 && !come_up)
		{
			(*pmask)->first_row = i;
			come_up = true;
		}
		if (n == 1)
		{
			int left = mask_matrix[i * 3 + 1];
			int right = mask_matrix[i * 3 + 2];
			(*(pmask))->left_x[i] = left;
			(*(pmask))->right_x[i] = right;
		}
		else
		{
			if (come_up&&(!the_last_one))
			{
				(*pmask)->last_row = i - 1;
				the_last_one = true;
			}
			(*(pmask))->left_x[i] = -1;
			(*(pmask))->right_x[i] = -1;
		}
		if ((i==row-1)&&(!the_last_one))
		{
			(*pmask)->last_row = row - 1;
		}
	}
}

struct offset_num**sort(struct offset_num**offset_num, int k, int n)
{
	struct offset_num**k_offsets = (struct offset_num**)malloc((k + 1)*sizeof(struct offset_num*));
	int num = 0;//k_offsets中现有的数量;
	for (int i = 0; i < n; i++)
	{
		if (num == k)
		{
			if ((*(offset_num + i))->num <= (*(k_offsets + k - 1))->num)
			{
				continue;
			}
			else
			{
				int place = dichotomy_find_place(*(offset_num + i), k_offsets, 0, k - 1);
				for (int j = k - 1; j > place; j--)
				{
					*(k_offsets + j) = (*(k_offsets + j - 1));
				}
				*(k_offsets + place) = (*(offset_num + i));
			}
		}
		else
		{
			if (num == 0)
			{
				*(k_offsets) = (*(offset_num + i));
				num = 1;
			}
			else if ((*(offset_num + i))->num <= (*(k_offsets + num - 1))->num)
			{
				*(k_offsets + num) = (*(offset_num + i));
				num++;
			}
			else
			{
				int place = dichotomy_find_place(*(offset_num + i), k_offsets, 0, num - 1);
				for (int j = num; j > place; j--)
				{
					(*(k_offsets + j)) = (*(k_offsets + j - 1));
				}
				*(k_offsets + place) = (*(offset_num + i));
				num++;
			}
		}
	}
	struct offset_num*offset = (struct offset_num*)malloc(sizeof(struct offset_num));
	offset->offset_x = 0;
	offset->offset_y = 0;
	offset->num = 0;
	(*(k_offsets + k)) = offset;
	return k_offsets;
}

int dichotomy_find_place(struct offset_num*offset, struct offset_num**k_offset, int first_one, int last_one)
{
	int num = offset->num;
	while (first_one<last_one)
	{
		if (last_one - first_one == 1)
		{
			return num>(*(k_offset + first_one))->num ? first_one : last_one;
		}
		int mid = (first_one + last_one) / 2;
		int mid_num = (*(k_offset + mid))->num;
		if (num > mid_num)
		{
			last_one = mid;
		}
		else
		{
			first_one = mid;
		}
	}
	return first_one;
}

void grid_label_optimize(unsigned int**image_red, unsigned int**image_green, unsigned int**image_blue,
struct grid_mask*mask, struct grid_mask*broaden_mask,
struct offset_num**offset_num, int k, int **pNewResult,int **pOldResult,
	int rect_first_row, int rect_last_row, int rect_left, int rect_right, int *edge, int edge_num,
	int image_col, int image_row, double *image_grid_x, double *image_grid_y)
{
	rect_first_row--;
	rect_last_row--;
	rect_left--;
	rect_right--;
	int rect_width = (rect_right - rect_left + 1);
	int rect_height = (rect_last_row - rect_first_row + 1);
	int num_pixels = rect_width*rect_height;
	int current_no = 0;//现在循环到了的mask编号;
	//int num_labels = k + 1;
	int num_labels = k + 1;
	double *grad_image = NULL;
	struct grid_extra_data*extra_data = (struct grid_extra_data*)malloc(sizeof(struct grid_extra_data));
	extra_data->pimage_red = &image_red;
	extra_data->pimage_green = &image_green;
	extra_data->pimage_blue = &image_blue;
	extra_data->pmask = &mask;
	extra_data->poffsets = &offset_num;
	extra_data->row = image_row;
	extra_data->col = image_col;
	extra_data->num_labels = num_labels;
	extra_data->rect_height = rect_height;
	extra_data->rect_width = rect_width;
	extra_data->rect_left = rect_left;
	extra_data->rect_up = rect_first_row;
	extra_data->pimage_grid_x = &image_grid_x;
	extra_data->pimage_grid_y = &image_grid_y;
	extra_data->pbroaden_mask = &broaden_mask;
	//OldResult:
	GCoptimizationGridGraph *gc = new GCoptimizationGridGraph((rect_right - rect_left + 1), (rect_last_row - rect_first_row + 1), num_labels);
	int *data_cost = new int[num_labels*num_pixels];
/*
	for (int i = 0; i < num_pixels; i++)
	{
		bool not_initialized = true;
		if (is_edge[i])
		{
			for (int j = 0; j < num_labels; j++)
			{
				if (j == num_labels - 1)
				{
					gc->setLabel(i, j);
					data_cost[i*num_labels + j] = 0;
				}
				else
				{
					data_cost[i*num_labels + j] = MaxDataCost;
				}
			}
		}
		else
		{
			struct coordinate coor = grid_get_coor(i, rect_width, rect_height, rect_left, rect_first_row);
			for (int j = 0; j < num_labels; j++)
			{
				int new_x = coor.x + (*(offset_num + j))->offset_x;
				int new_y = coor.y + (*(offset_num + j))->offset_y;
				if (grid_pixel_valid(new_x, new_y, mask, image_col, image_row))
				{
					data_cost[i*num_labels + j] = 0;
					if (not_initialized)
					{
						gc->setLabel(i, j);
						not_initialized = false;
					}
				}
				else
				{
					data_cost[i*num_labels + j] = MaxDataCost;
				}
			}
		}
	}*/
	for (int i = 0; i < num_pixels; i++)
	{
		bool not_initialized = true;
		struct coordinate coor = grid_get_coor(i, rect_width, rect_height, rect_left, rect_first_row);
		if (grid_pixel_valid(coor.x,coor.y,mask,image_col,image_row))
		{
			for (int j = 0; j < num_labels-1;j++)
			{
				data_cost[i*num_labels + j] = MaxDataCost;
			}
			data_cost[i*num_labels + num_labels] = 0;
			gc->setLabel(i, num_labels - 1);
		}
		else
		{
			for (int j = 0; j < num_labels; j++)
			{
				int new_x = coor.x + (*(offset_num + j))->offset_x;
				int new_y = coor.y + (*(offset_num + j))->offset_y;
				if (grid_pixel_valid(new_x, new_y, mask, image_col, image_row))
				{
					data_cost[i*num_labels + j] = 0;
					if (not_initialized)
					{
						gc->setLabel(i, j);
						not_initialized = false;
					}
				}
				else
				{
					data_cost[i*num_labels + j] = MaxDataCost;
				}
			}
		}
	}
	gc->setDataCost(data_cost);
	gc->setSmoothCost(OldGirdSmoothCostFn, extra_data);
	//gc->setLabelOrder(true);
	long int sum_cost = gc->compute_energy();
	std::cout << "old method cost before optimize is " << sum_cost;
	gc->swap(-1);
	sum_cost = gc->compute_energy();
	std::cout << "\n old method after optimizing the energy is " << sum_cost;
	(*pOldResult) = (int *)malloc(num_pixels*sizeof(int));
	for (int i = 0; i < num_pixels; i++)
	{
		(*pOldResult)[i] = gc->whatLabel(i);
	}
	delete gc;
	gc = NULL;
	//delete[]data_cost;
	//finish OldResult;

	//newResult:
	GCoptimizationGridGraph *new_gc = new GCoptimizationGridGraph((rect_right - rect_left + 1), (rect_last_row - rect_first_row + 1), num_labels);
	for (int i = 0; i < num_pixels; i++)
	{
		struct coordinate coor = grid_get_coor(i, rect_width, rect_height, rect_left, rect_first_row);
		if (grid_pixel_valid(coor.x,coor.y,broaden_mask,image_col,image_row))
		{
			new_gc->setLabel(i,((*pOldResult)[i]));
		}		else
		{
			for (int j = 0; j < num_labels; j++)
			{
				int new_x = coor.x + (*(offset_num + j))->offset_x;
				int new_y = coor.y + (*(offset_num + j))->offset_y;
				if (grid_pixel_valid(new_x, new_y, mask, image_col, image_row))
				{
					new_gc->setLabel(i, j);
					break;
				}
			}
		}
	}
	new_gc->setDataCost(NewDataCost, extra_data);
	new_gc->setSmoothCost(GirdSmoothCostFn, extra_data);
	long int new_sum_cost = new_gc->compute_energy();
	std::cout << "\n new method cost before optimize is " << new_sum_cost;
	new_gc->swap(-1);
	new_sum_cost = new_gc->compute_energy();
	std::cout << "\n new method after optimizing the energy is " << new_sum_cost;
	(*pNewResult) = (int *)malloc(num_pixels*sizeof(int));
	for (int i = 0; i < num_pixels; i++)
	{
		(*pNewResult)[i] = new_gc->whatLabel(i);
	}
	//delete[]data_cost;
	delete new_gc;
	new_gc = NULL;
	//finish newResult;
}

int GirdSmoothCostFn(int siteID1, int siteID2, int labelID1, int labelID2, void *extra_data)
{
	if (labelID1==labelID2)
	{
		return 0;
	}
	struct grid_extra_data*ex = (struct grid_extra_data*)extra_data;
	unsigned int **image_red = (*(ex->pimage_red));
	unsigned int **image_green = (*(ex->pimage_green));
	unsigned int **image_blue = (*(ex->pimage_blue));
	struct grid_mask*mask = (*(ex->pmask));
	struct offset_num**offsets = (*(ex->poffsets));
	int image_col = ex->col;
	int image_row = ex->row;
	int rect_width = ex->rect_width;
	int rect_height = ex->rect_height;
	int rect_left = ex->rect_left;
	int rect_up = ex->rect_up;
	//double *grid_image_x = (*(ex->pimage_grid_x));
	double *grid_image_y = (*(ex->pimage_grid_y));
	struct coordinate coor1 = grid_get_coor(siteID1, rect_width, rect_height, rect_left, rect_up);
	struct coordinate coor2 = grid_get_coor(siteID2, rect_width, rect_height, rect_left, rect_up);
	int offset_x1 = (*(offsets + labelID1))->offset_x;
	int offset_y1 = (*(offsets + labelID1))->offset_y;
	int offset_x2 = (*(offsets + labelID2))->offset_x;
	int offset_y2 = (*(offsets + labelID2))->offset_y;
	int new_s1_l1_x = coor1.x + offset_x1;
	int new_s1_l1_y = coor1.y + offset_y1;
	int new_s1_l2_x = coor1.x + offset_x2;
	int new_s1_l2_y = coor1.y + offset_y2;
	int new_s2_l1_x = coor2.x + offset_x1;
	int new_s2_l1_y = coor2.y + offset_y1;
	int new_s2_l2_x = coor2.x + offset_x2;
	int new_s2_l2_y = coor2.y + offset_y2;
	bool valid_s1_l1 = grid_pixel_valid(new_s1_l1_x, new_s1_l1_y, mask, image_col, image_row);
	bool valid_s1_l2 = grid_pixel_valid(new_s1_l2_x, new_s1_l2_y, mask, image_col, image_row);
	bool valid_s2_l1 = grid_pixel_valid(new_s2_l1_x, new_s2_l1_y, mask, image_col, image_row);
	bool valid_s2_l2 = grid_pixel_valid(new_s2_l2_x, new_s2_l2_y, mask, image_col, image_row);
	if (valid_s1_l1&&valid_s1_l2&&valid_s2_l1&&valid_s2_l2)
	{
		int cost = 0;
		cost += ((image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x])*(image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x])*(image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x])*(image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x])*(image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x])*(image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x])*(image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x]));
		/*double line_s1_l1, line_s1_l2, line_s2_l1, line_s2_l2;
		line_s1_l1 = (*(line_data + new_s1_l1_y*image_col + new_s1_l1_x));
		line_s1_l2 = (*(line_data + new_s1_l2_y*image_col + new_s1_l2_x));
		line_s2_l1 = (*(line_data + new_s2_l1_y*image_col + new_s2_l1_x));
		line_s2_l2 = (*(line_data + new_s2_l2_y*image_col + new_s2_l2_x));
		double d_line = (line_s1_l1 - line_s1_l2)*(line_s1_l1 - line_s1_l2) + (line_s2_l1-line_s2_l2)*(line_s2_l1-line_s2_l2);	
		cost += d_line;*/
		double grid_s1_l1, grid_s1_l2,grid_s2_l1,grid_s2_l2, d_grid;
		grid_s1_l1 = (*(grid_image_y + new_s1_l1_y*image_col + new_s1_l1_x));
		grid_s1_l2 = (*(grid_image_y + new_s1_l2_y*image_col + new_s1_l2_x));
		grid_s2_l1 = (*(grid_image_y + new_s2_l1_y*image_col + new_s2_l1_x));
		grid_s2_l2 = (*(grid_image_y + new_s2_l2_y*image_col + new_s2_l2_x));
		d_grid = (grid_s1_l1 - grid_s1_l2)*(grid_s1_l1 - grid_s1_l2) + (grid_s2_l1 - grid_s2_l2)*(grid_s2_l1 - grid_s2_l2);
		cost += d_grid;
		return cost;
	}
	else if (valid_s1_l1&&valid_s2_l2&&valid_s1_l2)
	{
		int cost = 0;
		cost += ((image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x])*(image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x])*(image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x])*(image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x]));
		//cost /= 2;
		double grid_s1_l1, grid_s1_l2, d_grid;
		grid_s1_l1 = (*(grid_image_y + new_s1_l1_y*image_col + new_s1_l1_x));
		grid_s1_l2 = (*(grid_image_y + new_s1_l2_y*image_col + new_s1_l2_x));
		d_grid = (grid_s1_l1 - grid_s1_l2)*(grid_s1_l1 - grid_s1_l2)*2;
		cost += d_grid;
		return cost;
	}
	else if (valid_s1_l1&&valid_s2_l2&&valid_s2_l1)
	{
		int cost = 0;
		cost += ((image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x])*(image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x])*(image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x])*(image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x]));
		//cost /= 2;
		double grid_s2_l1, grid_s2_l2, d_grid;
		grid_s2_l1 = (*(grid_image_y + new_s2_l1_y*image_col + new_s2_l1_x));
		grid_s2_l2 = (*(grid_image_y + new_s2_l2_y*image_col + new_s2_l2_x));
		d_grid = (grid_s2_l1 - grid_s2_l2)*(grid_s2_l1 - grid_s2_l2)*2;
		cost += d_grid;
		return cost;
	}
	else
	{
		return MaxSmoothCost;
	}
}

void inputdata(std::string file_name, unsigned int**image_8_red, unsigned int**image_8_green, unsigned int**image_8_blue,
	int *row_bin_num, int *col_bin_num, int **bin_num, int **bin_x, int **bin_y, int *image_row, int *image_col, int **mask_matrix,
	int **broaden_mask_matrix, int *edge_num, int **edge_data, int *first_row, int *last_row, int *left, int *right,
	double **grid_image_x, double **grid_image_y)
{
	std::ifstream in(file_name);
	int data;
	in >> (*row_bin_num) >> (*col_bin_num);
	int row_bin = (*row_bin_num);
	int col_bin = (*col_bin_num);
	(*bin_num) = (int *)malloc(row_bin*col_bin*sizeof(int));
	(*bin_x) = (int *)malloc(row_bin*col_bin*sizeof(int));
	(*bin_y) = (int *)malloc(row_bin*col_bin*sizeof(int));
	for (int i = 0; i < row_bin; i++)
	{
		for (int j = 0; j < col_bin; j++)
		{
			in >> (*((*bin_num) + i*col_bin + j));
		}
	}
	for (int i = 0; i < row_bin; i++)
	{
		for (int j = 0; j < col_bin; j++)
		{
			in >> (*((*bin_x) + i*col_bin + j));
		}
	}
	for (int i = 0; i < row_bin; i++)
	{
		for (int j = 0; j < col_bin; j++)
		{
			in >> (*((*bin_y) + i*col_bin + j));
		}
	}
	in >> (*image_row) >> (*image_col);
	int image_height = (*image_row);
	int image_width = (*image_col);
	(*(image_8_red)) = (unsigned int*)malloc(image_height*image_width*sizeof(int));
	(*(image_8_green)) = (unsigned int*)malloc(image_height*image_width*sizeof(int));
	(*(image_8_blue)) = (unsigned int*)malloc(image_height*image_width*sizeof(int));
	for (int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < image_width; j++)
		{
			in >> (*((*image_8_red) + i*image_width + j));
		}
	}
	int mask_col = 3;
	(*mask_matrix) = (int *)malloc(image_height * 3 * sizeof(int));
	for (int i = 0; i < image_height; i++)
	{
		in >> (*((*mask_matrix) + i * 3)) >> (*((*mask_matrix) + i * 3 + 1)) >> (*((*mask_matrix) + i * 3 + 2));
	}
	in >> (*first_row) >> (*last_row) >> (*left) >> (*right);
	in >> data;
	in >> (*edge_num);
	(*edge_data) = (int *)malloc((*edge_num)*sizeof(int));
	for (int i = 0; i < (*edge_num); i++)
	{
		in >> (*((*edge_data) + i));
	}
	for (int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < image_width; j++)
		{
			in >> (*((*image_8_green) + i*image_width + j));
		}
	}
	for (int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < image_width; j++)
		{
			in >> (*((*image_8_blue) + i*image_width + j));
		}
	}
	*grid_image_x = (double *)malloc(image_width*image_height*sizeof(double));
	*grid_image_y = (double *)malloc(image_width*image_height*sizeof(double));
	for (int i = 0; i < image_height;i++)
	{
		for (int j = 0; j < image_width;j++)
		{
			in >> (*((*grid_image_x) + i*image_width + j));
		}
	}
	for (int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < image_width; j++)
		{
			in >> (*((*grid_image_y) + i*image_width + j));
		}
	}
	(*broaden_mask_matrix) = (int *)malloc(image_height * 3 * sizeof(int));
	for (int i = 0; i < image_height; i++)
	{
		in >> (*((*broaden_mask_matrix) + i * 3)) >> (*((*broaden_mask_matrix) + i * 3 + 1)) >> (*((*broaden_mask_matrix) + i * 3 + 2));
	}
	in.close();
}

void deletedata(unsigned int**image_8_red, unsigned int**image_8_green, unsigned int**image_8_blue,
	int *row_bin_num, int *col_bin_num, int **bin_num, int **bin_x, int **bin_y, int *image_row, int *image_col, int **mask_matrix,
	int *edge_num, int **edge_data)
{
	free(*image_8_red);
	free(*image_8_green);
	free(*image_8_blue);
	free(*bin_num);
	free(*bin_x);
	free(*bin_y);
	free(*mask_matrix);
	free(*edge_data);
}

int NewDataCost(int siteID, int labelID, void *extra_data)
{
	struct grid_extra_data*ex = (struct grid_extra_data*)extra_data;
	unsigned int **image_red = (*(ex->pimage_red));
	unsigned int **image_green = (*(ex->pimage_green));
	unsigned int **image_blue = (*(ex->pimage_blue));
	struct grid_mask*mask = (*(ex->pmask));
	struct grid_mask*broaden_mask = (*(ex->pbroaden_mask));
	struct offset_num**offsets = (*(ex->poffsets));
	int image_width = ex->col;
	int image_height = ex->row;
	int rect_width = ex->rect_width;
	int rect_height = ex->rect_height;
	int rect_left = ex->rect_left;
	int rect_up = ex->rect_up;
	double *grid_data = (*(ex->pimage_grid_y));
	int offset_x = (*(offsets + labelID))->offset_x;
	int offset_y = (*(offsets + labelID))->offset_y;
	struct coordinate coor = grid_get_coor(siteID, rect_width, rect_height, rect_left, rect_up);
	if (grid_pixel_valid(coor.x,coor.y,broaden_mask,image_width,image_height))
	{
		if (offset_x==0&&offset_y==0)
		{
			return 0;
		}
		else
		{
			return MaxDataCost;
		}
	}
	else if (grid_pixel_valid(coor.x + offset_x, coor.y + offset_y, mask, image_width, image_height))
	{
		if (grid_pixel_valid(coor.x, coor.y, mask, image_width, image_height))
		{
			if (offset_x==0&&offset_y==0)
			{
				return MaxDataCost;
			}
			int cost = 0;
			int new_x = coor.x + offset_x;
			int new_y = coor.y + offset_y;
			cost += ((image_red[new_y][new_x] - image_red[coor.y][coor.x])*(image_red[new_y][new_x] - image_red[coor.y][coor.x]));
			cost += ((image_green[new_y][new_x] - image_green[coor.y][coor.x])*(image_green[new_y][new_x] - image_green[coor.y][coor.x]));
			cost += ((image_blue[new_y][new_x] - image_blue[coor.y][coor.x])*(image_blue[new_y][new_x] - image_blue[coor.y][coor.x]));
			double new_grid, old_grid;
			new_grid = (*(grid_data + new_y*image_width + new_x));
			old_grid = (*(grid_data + coor.y*image_width + coor.x));
			cost += ((new_grid - old_grid)*(new_grid - old_grid));
			return cost;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return MaxDataCost;
	}
}

int OldGirdSmoothCostFn(int siteID1, int siteID2, int labelID1, int labelID2, void *extra_data)
{
	if (labelID1 == labelID2)
	{
		return 0;
	}
	struct grid_extra_data*ex = (struct grid_extra_data*)extra_data;
	unsigned int **image_red = (*(ex->pimage_red));
	unsigned int **image_green = (*(ex->pimage_green));
	unsigned int **image_blue = (*(ex->pimage_blue));
	struct grid_mask*mask = (*(ex->pmask));
	struct offset_num**offsets = (*(ex->poffsets));
	int image_col = ex->col;
	int image_row = ex->row;
	int rect_width = ex->rect_width;
	int rect_height = ex->rect_height;
	int rect_left = ex->rect_left;
	int rect_up = ex->rect_up;
	double *grid_image_y = (*(ex->pimage_grid_y));
	struct coordinate coor1 = grid_get_coor(siteID1, rect_width, rect_height, rect_left, rect_up);
	struct coordinate coor2 = grid_get_coor(siteID2, rect_width, rect_height, rect_left, rect_up);
	int offset_x1 = (*(offsets + labelID1))->offset_x;
	int offset_y1 = (*(offsets + labelID1))->offset_y;
	int offset_x2 = (*(offsets + labelID2))->offset_x;
	int offset_y2 = (*(offsets + labelID2))->offset_y;
	int new_s1_l1_x = coor1.x + offset_x1;
	int new_s1_l1_y = coor1.y + offset_y1;
	int new_s1_l2_x = coor1.x + offset_x2;
	int new_s1_l2_y = coor1.y + offset_y2;
	int new_s2_l1_x = coor2.x + offset_x1;
	int new_s2_l1_y = coor2.y + offset_y1;
	int new_s2_l2_x = coor2.x + offset_x2;
	int new_s2_l2_y = coor2.y + offset_y2;
	bool valid_s1_l1 = grid_pixel_valid(new_s1_l1_x, new_s1_l1_y, mask, image_col, image_row);
	bool valid_s1_l2 = grid_pixel_valid(new_s1_l2_x, new_s1_l2_y, mask, image_col, image_row);
	bool valid_s2_l1 = grid_pixel_valid(new_s2_l1_x, new_s2_l1_y, mask, image_col, image_row);
	bool valid_s2_l2 = grid_pixel_valid(new_s2_l2_x, new_s2_l2_y, mask, image_col, image_row);
	if (valid_s1_l1&&valid_s1_l2&&valid_s2_l1&&valid_s2_l2)
	{
		int cost = 0;
		cost += ((image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x])*(image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x])*(image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x])*(image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x])*(image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x])*(image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x])*(image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x]));
		return cost;
	}
	else if (valid_s1_l1&&valid_s2_l2&&valid_s1_l2)
	{
		int cost = 0;
		cost += ((image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x])*(image_red[new_s1_l1_y][new_s1_l1_x] - image_red[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x])*(image_green[new_s1_l1_y][new_s1_l1_x] - image_green[new_s1_l2_y][new_s1_l2_x]));
		cost += ((image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x])*(image_blue[new_s1_l1_y][new_s1_l1_x] - image_blue[new_s1_l2_y][new_s1_l2_x]));
		return cost;
	}
	else if (valid_s1_l1&&valid_s2_l2&&valid_s2_l1)
	{
		int cost = 0;
		cost += ((image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x])*(image_red[new_s2_l1_y][new_s2_l1_x] - image_red[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x])*(image_green[new_s2_l1_y][new_s2_l1_x] - image_green[new_s2_l2_y][new_s2_l2_x]));
		cost += ((image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x])*(image_blue[new_s2_l1_y][new_s2_l1_x] - image_blue[new_s2_l2_y][new_s2_l2_x]));
		return cost;
	}
	else
	{
		return MaxSmoothCost;
	}
}

int mask_num(struct grid_mask*mask)
{
	int num = 0;
	for (int i = mask->first_row; i <= mask->last_row;i++)
	{
		num += (mask->right_x[i] - mask->left_x[i] + 1);
	}
	return num;
}
