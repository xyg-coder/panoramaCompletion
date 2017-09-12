//propogate时可能要考虑超过范围问题，暂时没考虑;

#include "kd_tree.h"
#include <stdio.h>
#include <malloc.h>
#include <memory.h>


#define MINPQ_INIT_NALLOCD 512
#define INT_MAX 100000

#ifndef ABS
#define ABS(x) ( ( (x) < 0 )? -(x) : (x) )
#endif


static inline int parent(int i)
{
	return (i - 1) / 2;
}

static inline int left(int i)
{
	return 2 * i + 1;
}

static inline int right(int i)
{
	return 2 * i + 2;
}

//默认是patch大小是8*8;
//n为image 中features的number;
//返回值是不受空洞影响的特征数;
int patch_match(struct feature* image_features,unsigned int**image,
	int image_cols, int image_rows, int n,struct mask*mask,mxArray*out[])
{
	struct feature**nbrs=NULL;
	struct kd_node*kd_root = NULL;
	struct feature**feature_coor = (struct feature**)calloc(n,sizeof(struct feature*));
	kd_root = kdtree_build(image_features, &n,feature_coor,image_cols-7,mask);
	int max_ = image_rows >= image_cols ? image_rows : image_cols;
	int thresh_2 = max_*max_ / (15.0*15.0);

	//调试;
/*
	struct feature*feat2 = patch_features;
	int k = kdtree_bbf_knn_for_patch1(kd_root, feat2, 2, &nbrs, 200);
	if (k == 2)
	{
		int d0 = image_pixel_dist(patch_image, image, 0, 0, nbrs[0]->y, nbrs[0]->x, 8, 8);
		int d1 = image_pixel_dist(patch_image, image, 0, 0, nbrs[1]->y, nbrs[1]->x, 8, 8);
		if (d0 > d1)
		{
			struct feature*tmp = nbrs[0];
			nbrs[0] = nbrs[1];
			nbrs[1] = tmp;
		}
		patch_features[0].fwd_match = *nbrs;
	}
	free(nbrs);
	int j = 1;
	struct feature*feat = patch_features + patch_image_cols-7;
	struct kd_node*leaf1 = propogate(image_features,feat2,NULL, patch_image,
		image, image_rows, image_cols-7, patch_image_rows, patch_image_cols,feature_coor);
	int k2 = kdtree_bbf_knn(kd_root, feat, leaf1, 2, &nbrs, 80);
	if (k2 == 2)
	{
		int d0 = image_pixel_dist(patch_image, image, j, 0, nbrs[0]->y, nbrs[0]->x, 8, 8);
		int d1 = image_pixel_dist(patch_image, image, j, 0, nbrs[1]->y, nbrs[1]->x, 8, 8);
		if (d0 > d1)
		{
			struct feature*tmp = nbrs[0];
			nbrs[0] = nbrs[1];
			nbrs[1] = tmp;
		}
		patch_features[j].fwd_match = *nbrs;
	}
	free(nbrs);*/
	//调试;
//	struct feature*feat = feature_coor[107 * (image_cols - 7) + 79];
	for (int i = 0; i < image_rows - 7; i++)
	{
		for (int j = 0; j < image_cols - 7; j++)
		{
			struct feature*feat = feature_coor[(image_cols - 7)*i + j];
			if (feat==NULL)
			{
				continue;
			}
			if (i == 0 && j == 0)
			{
				int k = kdtree_bbf_knn_for_patch1(kd_root, feat, 2, &nbrs, 200,thresh_2);
				if (k == 2)
				{
					int d0 = image_pixel_dist(image, image, 0, 0, nbrs[0]->y, nbrs[0]->x, 8, 8);
					int d1 = image_pixel_dist(image, image, 0, 0, nbrs[1]->y, nbrs[1]->x, 8, 8);
					if (d0 > d1)
					{
						struct feature*tmp = nbrs[0];
						nbrs[0] = nbrs[1];
						nbrs[1] = tmp;
					}
				}
				feat->fwd_match_num = k;
				feat->fwd_match = (struct feature**)malloc(k*sizeof(struct feature*));
				for (int iter = 0; iter < k;iter++)
				{
					*(feat->fwd_match + iter) = (*(nbrs + iter));
				}
				//merge_features(kd_root);
			}
			else if (i != 0 && j != 0)
			{
				struct kd_node*node = kd_root;
				while (node->m_leaf != 1)
				{
					int ki = node->ki;
					int kv = node->kv;
					if (feat->descr[ki] <= kv)
					{
						if (!node->kd_left&&node->kd_left->m_leaf!=1)
						{
							int a = 10;
						}
						node = node->kd_left;
					}
					else
					{
						if (!node->kd_right&&node->kd_right->m_leaf!=1)
						{
							int a = 10;
						}
						node = node->kd_right;
					}
				}
				struct kd_node*leaf1 = propogate(node,feat, feature_coor[(image_cols - 7)*(i - 1) + j], feature_coor[(image_cols - 7)*i + j-1], image,
					image, image_rows, image_cols-7,feature_coor,mask);
				int k = kdtree_bbf_knn(node, feat, leaf1, 2, &nbrs, 150,thresh_2);
				if (k == 2)
				{
					int d0 = image_pixel_dist(image, image, i, j, nbrs[0]->y, nbrs[0]->x, 8, 8);
					int d1 = image_pixel_dist(image, image, i, j, nbrs[1]->y, nbrs[1]->x, 8, 8);
					if (d0 > d1)
					{
						struct feature*tmp = nbrs[0];
						nbrs[0] = nbrs[1];
						nbrs[1] = tmp;
					}
				}
					feat->fwd_match_num = k;
					feat->fwd_match = (struct feature**)malloc(k*sizeof(struct feature*));
					for (int iter = 0; iter < k; iter++)
					{
						*(feat->fwd_match + iter) = (*(nbrs + iter));
					}
			}
			else if (i == 0 && j != 0)
			{
				struct kd_node*node = kd_root;
				while (node->m_leaf != 1)
				{
					int ki = node->ki;
					int kv = node->kv;
					if (feat->descr[ki] <= kv)
					{
						node = node->kd_left;
					}
					else
					{
						node = node->kd_right;
					}
				}
				struct kd_node*leaf1 = propogate(node,feat, NULL, feature_coor[(image_cols - 7)*i + j - 1], image,
					image, image_rows, image_cols-7, feature_coor,mask);
				int k = kdtree_bbf_knn(node, feat, leaf1, 2, &nbrs, 150,thresh_2);
				if (k == 2)
				{
					int d0 = image_pixel_dist(image, image, 0, j, nbrs[0]->y, nbrs[0]->x, 8, 8);
					int d1 = image_pixel_dist(image, image, 0, j, nbrs[1]->y, nbrs[1]->x, 8, 8);
					if (d0 > d1)
					{
						struct feature*tmp = nbrs[0];
						nbrs[0] = nbrs[1];
						nbrs[1] = tmp;
					}
				}
				feat->fwd_match_num = k;
				feat->fwd_match = (struct feature**)malloc(k*sizeof(struct feature*));
				for (int iter = 0; iter < k; iter++)
				{
					*(feat->fwd_match + iter) = (*(nbrs + iter));
				}
			}
			else
			{
				struct kd_node*node = kd_root;
				while (node->m_leaf != 1)
				{
					int ki = node->ki;
					int kv = node->kv;
					if (feat->descr[ki] <= kv)
					{
						node = node->kd_left;
					}
					else
					{
						node = node->kd_right;
					}
				}
				struct kd_node*leaf1 = propogate(node,feat, feature_coor[(image_cols - 7)*(i - 1) + j], NULL, image,
					image, image_rows, image_cols-7,feature_coor,mask);
				int k = kdtree_bbf_knn(node, feat, leaf1, 2, &nbrs, 150,thresh_2);
				if (k == 2)
				{
					int d0 = image_pixel_dist(image, image, i, 0, nbrs[0]->y, nbrs[0]->x, 8, 8);
					int d1 = image_pixel_dist(image, image, i, 0, nbrs[1]->y, nbrs[1]->x, 8, 8);
					if (d0 > d1)
					{
						struct feature*tmp = nbrs[0];
						nbrs[0] = nbrs[1];
						nbrs[1] = tmp;
					}
				}
				feat->fwd_match_num = k;
				feat->fwd_match = (struct feature**)malloc(k*sizeof(struct feature*));
				for (int iter = 0; iter < k; iter++)
				{
					*(feat->fwd_match + iter) = (*(nbrs + iter));
				}
			}
			free(nbrs);
		}
	}

	struct offset_num **diagram = dominant_offset_form(image_cols, image_rows, feature_coor, image_cols - 7, image_rows - 7);
	//struct offset_num**k_offsets = sort(diagram, 60, (image_rows * 2 + 1)*(image_cols * 2 + 1));
	mxArray*out1 = output_diagram(diagram, image_cols, image_rows);
	out[0] = out1;
	mxArray*out2 = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *data_out2 = (double *)mxGetData(out2);
	*data_out2 = image_cols;
	out[1] = out2;
	mxArray*out3 = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *data_out3 = (double *)mxGetData(out3);
	*data_out3 = image_rows;
	out[2] = out3;
	kdtree_release(kd_root);
	
	for (int i = 0; i < (image_rows / 2 * 2 + 1)*(image_cols / 2 * 2 + 1); i++)
	{
		free(*(diagram + i));
	}
	//free(k_offsets);
	free(diagram);
	free(feature_coor);
	return n;
}

//n: the number of the features
//mode 为1时，是第一个patch的kd_tree的build过程;
struct kd_node* kdtree_build(struct feature* features, int*nn,struct feature**feature_coor,int image_col,struct mask*mask)
{
	struct kd_node* kd_root;
	int n = (*nn) - 1;
	if (!features || n <= 0)
	{
		fprintf(stderr, "Warning: kdtree_build(): no features, %s, line %d\n",
			__FILE__, __LINE__);
		return NULL;
	}
	int i1 = 0, j1 = n - 1;
	struct feature*front = NULL, *last = features + j1;
	int num_mask = 0;
	for (int i = 0; i < mask->row - 7; i++)
	{
		for (int j = 0; j < image_col; j++)
		{
			if (!valid(j, i, mask))
			{
				if (front == NULL)
				{
					front = features + i*image_col + j;
				}
				feature_coor[i*image_col + j] = NULL;
				(features + i*image_col + j)->valid = false;
				num_mask++;
			}
		}
	}
	if (front)
	{
		while (front != last&&front)
		{
 			while (front->valid&&front!=last)
			{
				front++;
			}
			n--;
			while (!last->valid&&last!=front)
			{
				n--;
				last--;
			}
			*front = *last;
			if (last-front==1||last==front)
			{
				break;
			}
			front++;
			last--;
		}
	}

	kd_root = kd_node_init(features, n);
	expand_kd_node_subtree(kd_root,feature_coor,image_col,1);
	*nn = n;
	/*for (int i = 0; i < image_col;i++)
	{
		for (int j = 0; j < 431;j++)
		{
			struct feature* f1 = feature_coor[353 * j + i];
			if (!f1)
			{
				int uu = 10;
			}
		}
	}*/
	//struct feature* f2 = feature_coor[353 * 9 + 65];
	return kd_root;
}

/*
Initializes a kd tree node with a set of features.  The node is not
expanded, and no ordering is imposed on the features.

@param features an array of image features
@param n number of features

@return Returns an unexpanded kd-tree node.
*/
static struct kd_node* kd_node_init(struct feature* features, int n)
{
	struct kd_node* kd_node;

	kd_node = (struct kd_node*)malloc(sizeof(struct kd_node));
	memset(kd_node, 0, sizeof(struct kd_node));
	kd_node->ki = -1;
	kd_node->features = features;
	kd_node->n = n;

	return kd_node;
}

/*
Recursively expands a specified kd tree node into a tree whose leaves
contain m entries.

@param kd_node an unexpanded node in a kd tree
@param mode 1则是还在寻找m_leaf，0则是已经找到;
*/
static void expand_kd_node_subtree(struct kd_node* kd_node,struct feature**feature_coor,int image_col,int mode)
{
	/*if (kd_node->n==3)
	{
		struct feature*feat = kd_node->features;
		struct feature*feat2 = feat + 1;
		struct feature*feat3 = feat2 + 1;
		if (feat->x == 64 && feat->y == 9||feat2->x==64&&feat2->y==9||feat3->x==64&&feat3->y==9)
		{
			int m = 10;
		}
	}*/
	/* base case: leaf node */
	if (kd_node->n ==1||kd_node->n==0)
	{
		kd_node->leaf = 1;
		for (int i = 0; i < kd_node->n;i++)
		{
			kd_node->features[i].p_leaf = kd_node;
		}
		if (kd_node->n==1)
		{
			int x = kd_node->features[0].x;
			int y = kd_node->features[0].y;
			feature_coor[y*image_col + x] = &(kd_node->features[0]);
		}
		return;
	}
	else if (/*kd_node->n>=8&&*/kd_node->n<=64&&mode)
	{
		kd_node->m_leaf = 1;
		for (int i = 0; i < kd_node->n;i++)
		{
			kd_node->features[i].p_m_leaf = kd_node;
		}
		mode = 0;
	}

	assign_part_key(kd_node);
	partition_features(kd_node,feature_coor,image_col);

	if (kd_node->kd_left)
		expand_kd_node_subtree(kd_node->kd_left,feature_coor,image_col,mode);
	if (kd_node->kd_right)
		expand_kd_node_subtree(kd_node->kd_right,feature_coor,image_col,mode);

	
}


/*
Determines the descriptor index at which and the value with which to
partition a kd tree node's features.

@param kd_node a kd tree node
*/
static void assign_part_key(struct kd_node* kd_node)
{
	struct feature* features;
	int kv, x, mean, var, var_max = 0;
	int* tmp;
	int d, n, i, j, ki = 0;

	features = kd_node->features;
	n = kd_node->n;
	d = features[0].d;

	/* partition key index is that along which descriptors have most variance */
	for (j = 0; j < d; j++)
	{
		mean = var = 0;
		for (i = 0; i < n; i++)
			mean += features[i].descr[j];
		mean /= n;
		for (i = 0; i < n; i++)
		{
			x = features[i].descr[j] - mean;
			var += x * x;
		}
		var /= n;

		if (var > var_max)
		{
			ki = j;
			var_max = var;
		}
	}

	/* partition key value is median of descriptor values at ki */
	tmp = (int *)calloc(n, sizeof(double));
	for (i = 0; i < n; i++)
		tmp[i] = features[i].descr[ki];
	kv = median_select(tmp, n);
	free(tmp);

	kd_node->ki = ki;
	kd_node->kv = kv;
}

/*
Finds the median value of an array.  The array's elements are re-ordered
by this function.

@param array an array; the order of its elelemts is reordered
@param n number of elements in array

@return Returns the median value of array.
*/
static double median_select(int* array, int n)
{
	return rank_select(array, n, (n - 1) / 2);
}

/*
Finds the element of a specified rank in an array using the linear time
median-of-medians algorithm by Blum, Floyd, Pratt, Rivest, and Tarjan.
The elements of the array are re-ordered by this function.

@param array an array; the order of its elelemts is reordered
@param n number of elements in array
@param r the zero-based rank of the element to be selected

@return Returns the element from array with zero-based rank r.
*/
static double rank_select(int* array, int n, int r)
{
	int* tmp, med;
	int gr_5, gr_tot, rem_elts, i, j, *index;

	/* base case */
	if (n == 1)
		return array[0];

	/* divide array into groups of 5 and sort them */
	gr_5 = n / 5;
	gr_tot = gr_5 + 1;
	rem_elts = n % 5;
	tmp = array;
	for (i = 0; i < gr_5; i++)
	{
		insertion_sort(tmp, 5);
		tmp += 5;
	}
	insertion_sort(tmp, rem_elts);

	/* recursively find the median of the medians of the groups of 5 */
	tmp = (int *)calloc(gr_tot, sizeof(double));
	for (i = 0, j = 2; i < gr_5; i++, j += 5)
		tmp[i] = array[j];
	if (rem_elts)
		tmp[i++] = array[n - 1 - rem_elts / 2];
	med = rank_select(tmp, i, (i - 1) / 2);
	free(tmp);

	/* partition around median of medians and recursively select if necessary */
	index = partition_array(array, n, med);
	if (r >= index[0]&&r<=index[1])
		return med;
	else if (r < index[0])
		return rank_select(array, index[0], r);
	else
	{
		array += (index[1] + 1);
		return rank_select(array, (n - index[1] - 1), (r - index[1] - 1));
	}
	free (index);
}

/*
Sorts an array in place into increasing order using insertion sort.

@param array an array
@param n number of elements
*/
static void insertion_sort(int* array, int n)
{
	int k;
	int i, j;

	for (i = 1; i < n; i++)
	{
 		k = array[i];
		j = i - 1;
		while (j >= 0 && array[j] > k)
		{
			array[j + 1] = array[j];
			j -= 1;
		}
		array[j + 1] = k;
	}
}

/*
Partitions an array around a specified value.

@param array an array
@param n number of elements
@param pivot value around which to partition

@return Returns index of the pivot after partitioning
*/
static int* partition_array(int* array, int n, int pivot)
{
	int tmp;
	int p, i, j;
	int n_med = 0;//等于piyot的个数;
	int *index = (int *)malloc(2 * sizeof(int));
	i = -1;
	for (j = 0; j < n; j++)
		if (array[j] <= pivot)
		{
		tmp = array[++i];
		array[i] = array[j];
		array[j] = tmp;
		if (array[i] == pivot)
		{
			if (n_med == 0)
			{
				index[0] = i;
			}
			p = i;
			n_med++;
		}
		}
	array[p] = array[i];
	array[i] = pivot;
	index[1] = i;
	return index;
}

/*
Partitions the features at a specified kd tree node to create its two
children.

@param kd_node a kd tree node whose partition key is set
*/
static void partition_features(struct kd_node* kd_node,struct feature**feature_coor,int image_col)
{
	struct feature* features, tmp;
	int kv;
	int n, ki, p, i, j = -1;

	features = kd_node->features;
	n = kd_node->n;
	ki = kd_node->ki;
	kv = kd_node->kv;
	for (i = 0; i < n; i++)
		if (features[i].descr[ki] <= kv)
		{
		tmp = features[++j];
		features[j] = features[i];
		features[i] = tmp;
		if (features[j].descr[ki] == kv)
			p = j;
		}
	tmp = features[p];
	features[p] = features[j];
	features[j] = tmp;

	/* if all records fall on same side of partition, make node a leaf */
	if (j == n - 1||j == n - 2)
	{
		kd_node->leaf = 1;
		for (int i = 0; i < n;i++)
		{
			if ((features+i)->p_m_leaf==NULL)
			{
				kd_node->m_leaf = 1;
				(features + i)->p_m_leaf = kd_node;
			}
			(features+i)->p_leaf = kd_node;
			int x = (features + i)->x;
			int y = (features + i)->y;
			feature_coor[y*image_col + x] = features + i;
		}
		return;
	}

	kd_node->kd_left = kd_node_init(features, j + 1);
	kd_node->kd_right = kd_node_init(features + (j + 1), (n - j - 1));
}

//用于在第一次kd_tree构造之后将不必要的kd_tree进行合并;
/*
static void merge_features(struct kd_node*kd_node)
{
	if (kd_node->n>64)
	{
		merge_features(kd_node->kd_left);
		merge_features(kd_node->kd_right);
	}
	else
	{
		for (int i = 0; i < kd_node->n;i++)
		{
			kd_node->features[i].feature_leaf = kd_node;
			kdtree_release(kd_node->kd_left);
			kdtree_release(kd_node->kd_right);
		}
	}
}*/

/*
Finds an image feature's approximate k nearest neighbors in a kd tree using
Best Bin First search.

@param kd_root root of an image feature kd tree
@param feat image feature for whose neighbors to search
@param k number of neighbors to find
@param nbrs pointer to an array in which to store pointers to neighbors
in order of increasing descriptor distance
@param max_nn_chks search is cut off after examining this many tree entries

@return Returns the number of neighbors found and stored in nbrs, or
-1 on error.
*/
//此为第一个patch采用的方法;
int kdtree_bbf_knn_for_patch1(struct kd_node* kd_root, struct feature* feat, int k,
struct feature*** nbrs, int max_nn_chks,int distant_thresh_2)
{
	struct kd_node* expl;
	struct min_pq* min_pq;
	struct feature* tree_feat, ** _nbrs;
	struct bbf_data* bbf_data;
	int i, t = 0, n = 0;

	if (!nbrs || !feat || !kd_root)
	{
		fprintf(stderr, "Warning: NULL pointer error, %s, line %d\n",
			__FILE__, __LINE__);
		return -1;
	}

	_nbrs = (struct feature**)calloc(k, sizeof(struct feature*));
	//以下为回溯做准备;
	min_pq = minpq_init();
	minpq_insert(min_pq, kd_root, 0);//将kd_root插入minpq队列;
	while (min_pq->n > 0 && t < max_nn_chks)
	{
		expl = (struct kd_node*)minpq_extract_min(min_pq);//提取优先级最高的节点,复制给当前节点;
		if (!expl)
		{
			fprintf(stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
				__FILE__, __LINE__);
			goto fail;
		}

		expl = explore_to_leaf(expl, feat, min_pq);//从当前搜索节点expl一直搜索到叶子节点，
		//搜索过程中将未搜索的节点根据优先级放入队列，返回值为叶子节点;

		if (!expl)
		{
			fprintf(stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
				__FILE__, __LINE__);
			goto fail;
		}

		//比较查找最近邻;
		for (i = 0; i < expl->n; i++)
		{
			tree_feat = &expl->features[i];
			bbf_data =(struct bbf_data*) malloc(sizeof(struct bbf_data));
			if (!bbf_data)
			{
				fprintf(stderr, "Warning: unable to allocate memory,"
					" %s line %d\n", __FILE__, __LINE__);
				goto fail;
			}
			bbf_data->old_data = tree_feat->feature_data;//保存第i个特征点的feature_data域以前的值;
			if ((feat->x-tree_feat->x)*(feat->x-tree_feat->x)+(feat->y-tree_feat->y)*(feat->y-tree_feat->y)<=distant_thresh_2)
			{
				continue;
			}
			bbf_data->d = descr_dist_sq(feat, tree_feat);//当前搜索点和目标点之间的欧氏距离;
			tree_feat->feature_data = bbf_data;
			n += insert_into_nbr_array(tree_feat, _nbrs, n, k);//当最近邻数组中元素已经达到k时，继续插入元素个数不会增加，但是会更新元素值;
		}
		t++;
	}

	minpq_release(&min_pq,feat);
	for (i = 0; i < n; i++)
	{
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		free(bbf_data);
	}
	*nbrs = _nbrs;
	return n;

fail:
	minpq_release(&min_pq,feat);
	for (i = 0; i < n; i++)
	{
		bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
		_nbrs[i]->feature_data = bbf_data->old_data;
		free(bbf_data);
	}
	free(_nbrs);
	*nbrs = NULL;
	return -1;
}

/*
// * /针对除了第一个patch所采用的方法;
int kdtree_bbf_knn(struct kd_node* kd_root,struct kd_node*leaf1, struct feature* feat, int k,
struct feature***nbrs)
{
	struct kd_node*expl;
	struct feature* tree_feat, **_nbrs;
	struct bbf_data*bbf_data;
}* /*/

/*
De-allocates memory held by a kd tree

@param kd_root pointer to the root of a kd tree
*/
void kdtree_release(struct kd_node* kd_root)
{
	if (!kd_root)
		return;
	kdtree_release(kd_root->kd_left);
	kdtree_release(kd_root->kd_right);
	free(kd_root);
}

/*
Creates a new minimizing priority queue.
*/
struct min_pq* minpq_init()
{
	struct min_pq* min_pq;

	min_pq = (struct min_pq*)malloc(sizeof(struct min_pq));
	min_pq->pq_array =(struct pq_node*) calloc(MINPQ_INIT_NALLOCD, sizeof(struct pq_node));
	min_pq->nallocd = MINPQ_INIT_NALLOCD;
	min_pq->n = 0;

	return min_pq;
}


/**
Inserts an element into a minimizing priority queue.

@param min_pq a minimizing priority queue
@param data the data to be inserted
@param key the key to be associated with \a data

@return Returns 0 on success or 1 on failure.
*/
int minpq_insert(struct min_pq* min_pq, void* data, int key)
{
	int n = min_pq->n;

	/* double array allocation if necessary */
	if (min_pq->nallocd == n)
	{
		min_pq->nallocd = array_double((void**)&min_pq->pq_array,
			min_pq->nallocd,
			sizeof(struct pq_node));
		if (!min_pq->nallocd)
		{
			fprintf(stderr, "Warning: unable to allocate memory, %s, line %d\n",
				__FILE__, __LINE__);
			return 1;
		}
	}


	min_pq->pq_array[n].data = data;
	min_pq->pq_array[n].key = INT_MAX;
	decrease_pq_node_key(min_pq->pq_array, min_pq->n, key);
	min_pq->n++;

	return 0;
}

/*
Decrease a minimizing pq element's key, rearranging the pq if necessary

@param pq_array minimizing priority queue array
@param i index of the element whose key is to be decreased
@param key new value of element <EM>i</EM>'s key; if greater than current
key, no action is taken
*/
static void decrease_pq_node_key(struct pq_node* pq_array, int i, int key)
{
	struct pq_node tmp;

	if (key > pq_array[i].key)
		return;

	pq_array[i].key = key;
	while (i > 0 && pq_array[i].key < pq_array[parent(i)].key)
	{
		tmp = pq_array[parent(i)];
		pq_array[parent(i)] = pq_array[i];
		pq_array[i] = tmp;
		i = parent(i);
	}
}


/*
Removes and returns the element of a minimizing priority queue with the
smallest key.

@param min_pq a minimizing priority queue

@return Returns the element of \a min_pq with the smallest key of NULL
if \a min_pq is empty
*/
void* minpq_extract_min(struct min_pq* min_pq)
{
	void* data;

	if (min_pq->n < 1)
	{
		fprintf(stderr, "Warning: PQ empty, %s line %d\n", __FILE__, __LINE__);
		return NULL;
	}
	data = min_pq->pq_array[0].data;
	min_pq->n--;
	min_pq->pq_array[0] = min_pq->pq_array[min_pq->n];
	restore_minpq_order(min_pq->pq_array, 0, min_pq->n);

	return data;
}

/*
Recursively restores correct priority queue order to a minimizing pq array

@param pq_array a minimizing priority queue array
@param i index at which to start reordering
@param n number of elements in \a pq_array
*/
static void restore_minpq_order(struct pq_node* pq_array, int i, int n)
{
	struct pq_node tmp;
	int l, r, min = i;

	l = left(i);
	r = right(i);
	if (l < n)
		if (pq_array[l].key < pq_array[i].key)
			min = l;
	if (r < n)
		if (pq_array[r].key < pq_array[min].key)
			min = r;

	if (min != i)
	{
		tmp = pq_array[min];
		pq_array[min] = pq_array[i];
		pq_array[i] = tmp;
		restore_minpq_order(pq_array, min, n);
	}
}

/*
Explores a kd tree from a given node to a leaf.  Branching decisions are
made at each node based on the descriptor of a given feature.  Each node
examined but not explored is put into a priority queue to be explored
later, keyed based on the distance from its partition key value to the
given feature's desctiptor.
从给定的kd tree追溯到叶子，检查但是没有explore的节点放入优先队列当中优先级的排列顺序为空间距离;，?
@param kd_node root of the subtree to be explored
@param feat feature upon which branching decisions are based
@param min_pq a minimizing priority queue into which tree nodes are placed
as described above

@return Returns a pointer to the leaf node at which exploration ends or
NULL on error.
*/
static struct kd_node* explore_to_leaf(struct kd_node* kd_node,
struct feature* feat,
struct min_pq* min_pq)
{
	struct kd_node* unexpl, *expl = kd_node;
	int kv;
	int ki;

	while (expl  &&  !expl->leaf)
	{
		ki = expl->ki;
		kv = expl->kv;

		if (ki >= feat->d)
		{
			fprintf(stderr, "Warning: comparing imcompatible descriptors, %s" \
				" line %d\n", __FILE__, __LINE__);
			return NULL;
		}
		if (feat->descr[ki] <= kv)
		{
			unexpl = expl->kd_right;
			expl = expl->kd_left;
		}
		else
		{
			unexpl = expl->kd_left;
			expl = expl->kd_right;
		}

		if (minpq_insert(min_pq, unexpl, ABS(kv - feat->descr[ki])))
		{
			fprintf(stderr, "Warning: unable to insert into PQ, %s, line %d\n",
				__FILE__, __LINE__);
			return NULL;
		}
	}

	return expl;
}

/*
Calculates the squared Euclidian distance between two feature descriptors.

@param f1 first feature
@param f2 second feature

@return Returns the squared Euclidian distance between the descriptors of
f1 and f2.
*/
int descr_dist_sq(struct feature* f1, struct feature* f2)
{
	int diff, dsq = 0;
	int* descr1, *descr2;
	int i, d;

	
	descr1 = f1->descr;
	descr2 = f2->descr;
	d = f2->d;

	for (i = 0; i < d; i++)
	{
		diff = descr1[i] - descr2[i];
		dsq += diff*diff;
	}
	return dsq;
}

/*
Inserts a feature into the nearest-neighbor array so that the array remains
in order of increasing descriptor distance from the search feature.

@param feat feature to be inderted into the array; it's feature_data field
should be a pointer to a bbf_data with d equal to the squared descriptor
distance between feat and the search feature
@param nbrs array of nearest neighbors neighbors
@param n number of elements already in nbrs and
@param k maximum number of elements in nbrs

@return If feat was successfully inserted into nbrs, returns 1; otherwise
returns 0.
*/
static int insert_into_nbr_array(struct feature* feat, struct feature** nbrs,
	int n, int k)
{

	struct bbf_data* fdata, *ndata;
	double dn, df;
	int i, ret = 0;

	if (n == 0)
	{
		nbrs[0] = feat;
		return 1;
	}

	/* check at end of array */
	fdata = (struct bbf_data*)feat->feature_data;
	df = fdata->d;
	ndata = (struct bbf_data*)nbrs[n - 1]->feature_data;
	dn = ndata->d;
	if (df >= dn)
	{
		if (n == k)
		{
			feat->feature_data = fdata->old_data;
			free(fdata);
			return 0;
		}
		nbrs[n] = feat;
		return 1;
	}

	/* find the right place in the array */
	if (n < k)
	{
		nbrs[n] = nbrs[n - 1];
		ret = 1;
	}
	else
	{
		nbrs[n - 1]->feature_data = ndata->old_data;
		free(ndata);
	}
	i = n - 2;
	while (i >= 0)
	{
		ndata = (struct bbf_data*)nbrs[i]->feature_data;
		dn = ndata->d;
		if (dn <= df)
			break;
		nbrs[i + 1] = nbrs[i];
		i--;
	}
	i++;
	nbrs[i] = feat;

	return ret;
}

/*
De-allocates the memory held by a minimizing priorioty queue

@param min_pq pointer to a minimizing priority queue
*/
void minpq_release(struct min_pq** min_pq,struct feature*feat)
{
	if (!min_pq)
	{
		fprintf(stderr, "Warning: NULL pointer error, %s line %d\n", __FILE__,
			__LINE__);
		return;
	}
	if ((*min_pq) && (*min_pq)->pq_array)
	{
		free((*min_pq)->pq_array);
		free(*min_pq);
		*min_pq = NULL;
	}
}

//欧几里得距离;
//patch_y,patch_x:the location of the first element of the patch
//image_x.image_y:the location of the first element of the image
int image_pixel_dist(unsigned int**patch_pixel, unsigned int**image_pixel,
	int patch_y, int patch_x, int image_y, int image_x,int width,int height)
{
	int sum = 0;
	for (int i = 0; i < height;i++)
	{
		for (int j = 0; j < width;j++)
		{
			int dist = patch_pixel[patch_y + i][patch_x + j] - image_pixel[image_y + i][image_x + j];
			sum += (dist*dist);
		}
	}
	return sum;
}

/*
//根据邻近已知good match确定下一个leaf;
//patch_width:补丁的像素宽度;
struct kd_node*kdroot_build(struct kd_node*kd_root,struct feature*patch_features,struct feature *up_feature,
int patch_width,struct feature*left_feature, unsigned int**patch_pixel, unsigned int**image_pixel)
{
	struct kd_node* leaf0, *leaf1,*new_root;
	if (!up_feature&&!left_feature)
	{
		return NULL;
	}
	else if (up_feature&&left_feature)
	{
		struct feature*left_feature1 = left_feature->fwd_match[0];
		struct feature*left_feature2 = left_feature->fwd_match[1];
		struct feature*up_feature1 = up_feature->fwd_match[0];
		struct feature*up_feature2 = up_feature->fwd_match[1];
		int left_x1 = left_feature1->x;
		int left_x2 = left_feature2->x;
		int left_y1 = left_feature1->y;
		int left_y2 = left_feature2->y;
		int up_x1 = up_feature1->x;
		int up_x2 = up_feature2->x;
		int up_y1 = up_feature1->y;
		int up_y2 = up_feature2->y; 
		int x = left_feature->x + 1;
		int y = left_feature->y;
		int dist[4];
		dist[0] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, left_y1, left_x1 + 1,8,8);
		dist[1] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, left_y2, left_x2 + 1, 8, 8);
		dist[2] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, up_y1 + 1, up_x1, 8, 8);
		dist[3] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, up_y2 + 1, up_x2, 8, 8);
		if (dist[0]<=dist[1]&&dist[0]<=dist[2]&&dist[0]<=dist[3])
		{
			leaf1 = left_feature1->p_m_leaf;
		}
		else if (dist[1] <= dist[0] && dist[1] <= dist[2]&&dist[1]<=dist[3])
		{
			leaf1 = left_feature2->p_m_leaf;
		}
		else if (dist[2] <= dist[1] && dist[2] <= dist[0] && dist[2] <= dist[3])
		{
			leaf1 = up_feature1->p_m_leaf;
		}
		else
		{
			leaf1 = up_feature2->p_m_leaf;
		}
		leaf0 = kd_root;
		while (leaf0->m_leaf!=1)
		{
			int ki = leaf0->ki;
			int kv = leaf0->kv;
			struct feature*feat = patch_features + patch_width * y + x;
			if (feat->descr[ki]<=kv)
			{
				leaf0 = leaf0->kd_left;
			}
			else
			{
				leaf0 = leaf0->kd_right;
			}
		}
	}
	else if (up_feature&&(!left_feature))
	{
		struct feature*up_feature1 = up_feature->fwd_match[0];
		struct feature*up_feature2 = up_feature->fwd_match[1];
		int up_x1 = up_feature1->x;
		int up_x2 = up_feature2->x;
		int up_y1 = up_feature1->y;
		int up_y2 = up_feature2->y;
		int x = up_feature->x;
		int y = up_feature->y + 1;
		int dist[2];
		dist[0] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, up_y1 + 1, up_x1, 8, 8);
		dist[1] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, up_y2 + 1, up_x2, 8, 8);
		if (dist[0]<dist[1])
		{
			leaf1 = up_feature1->p_m_leaf;
		}
		else
		{
			leaf1 = up_feature2->p_m_leaf;
		}
		leaf0 = kd_root;
		while (leaf0->m_leaf != 1)
		{
			int ki = leaf0->ki;
			int kv = leaf0->kv;
			struct feature*feat = patch_features + patch_width * y + x;
			if (feat->descr[ki] <= kv)
			{
				leaf0 = leaf0->kd_left;
			}
			else
			{
				leaf0 = leaf0->kd_right;
			}
		}
	}
	else
	{
		struct feature*left_feature1 = left_feature->fwd_match[0];
		struct feature*left_feature2 = left_feature->fwd_match[1];
		int left_x1 = left_feature1->x;
		int left_x2 = left_feature2->x;
		int left_y1 = left_feature1->y;
		int left_y2 = left_feature2->y;
		int x = left_feature->x + 1;
		int y = up_feature->y;
		int dist[2];
		dist[0] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, left_y1, left_x1+1, 8, 8);
		dist[1] = image_pixel_dist(patch_pixel, image_pixel, y,
			x, left_y2, left_x2 + 1, 8, 8);
		if (dist[0] < dist[1])
		{
			leaf1 = left_feature1->p_m_leaf;
		}
		else
		{
			leaf1 = left_feature2->p_m_leaf;
		}
		leaf0 = kd_root;
		while (leaf0->m_leaf != 1)
		{
			int ki = leaf0->ki;
			int kv = leaf0->kv;
			struct feature*feat = patch_features + patch_width * y + x;
			if (feat->descr[ki] <= kv)
			{
				leaf0 = leaf0->kd_left;
			}
			else
			{
				leaf0 = leaf0->kd_right;
			}
		}
	}
}*/

//以下为将其他邻近good-patch propogate到新的patch的过程;
//以距离进行排序，选择一个合适的顺序添加到kd_node**中；
/*
参数：features:image的feature数组;
up_feature和left_feature为patch已经匹配过的feature;
patch_image和image分别为patch的图像和image的图像;
image_row和image_col, patch_row和patch_col都为对应的长宽;*/
struct kd_node*propogate(struct kd_node*node,struct feature*features,struct feature*up_feature, struct feature*left_feature,
	unsigned int**patch_image, unsigned int** image,
	int image_row,int image_col/*,int patch_row,int patch_col*/,struct feature**feature_coor,struct mask*mask)
{
	if (!up_feature&&!left_feature)
	{
		fprintf(stderr, "Warning: NULL pointer error, %s, line %d\n",
			__FILE__, __LINE__);
		return NULL;
	}
	else if (up_feature&&left_feature)
	{
		int patch_x = up_feature->x;
		int patch_y = up_feature->y+1;
		/*if (patch_x==140&&patch_y==221)
		{
			int m = 20;
		}*/
		int patch_num = up_feature->fwd_match_num + left_feature->fwd_match_num;
		int *new_x = (int *)malloc(patch_num*sizeof(int)), *new_y = (int *)malloc(patch_num*sizeof(int));
		for (int i = 0; i < left_feature->fwd_match_num;i++)
		{
			new_x[i] = (*(left_feature->fwd_match + i))->x + 1;
			new_y[i] = (*(left_feature->fwd_match + i))->y;
		}
		for (int i = 0; i < up_feature->fwd_match_num;i++)
		{
			new_x[i + left_feature->fwd_match_num] = (*(up_feature->fwd_match + i))->x;
			new_y[i + left_feature->fwd_match_num] = (*(up_feature->fwd_match + i))->y + 1;
		}
		int *dist = (int *)malloc(patch_num*sizeof(int));
		struct feature**feature_array = (struct feature**)malloc(patch_num*sizeof(struct feature*));
		int leaf_num = 0;//与node不同的m_leaf的数目;
		for (int i = 0; i < patch_num;i++)
		{
			int y = new_y[i];
			int x = new_x[i];
			int num = new_y[i] * image_col + new_x[i];
			if (new_x[i]<image_col&&new_y[i]<image_row - 7 && feature_coor[num] && feature_coor[num]->p_m_leaf != node)
			{
				dist[leaf_num] = image_pixel_dist(patch_image, image, patch_y, patch_x, new_y[i], new_x[i], 8, 8);
				feature_array[leaf_num] = feature_coor[num];
				leaf_num++;
			}
		}
		if (leaf_num)
		{
			int index=0, min_dist=dist[0];
			for (int i = 1; i < leaf_num;i++)
			{
				if (dist[i]<=min_dist)
				{
					min_dist = dist[i];
					index = i;
				}
			}
			return feature_array[index]->p_m_leaf;
		}
		else
		{
			return NULL;
		}
		free(new_x);
		free(new_y);
		free(dist);
		free(feature_array);
	/*	if (left_new_match1_x<image_col&&left_new_match2_x<image_col
			&&up_new_match1_y<image_row-7&&up_new_match2_y<image_row-7)
		{
			int dist[4];
			dist[0] = image_pixel_dist(patch_image, image,patch_y,patch_x,left_new_match1_y,left_new_match1_x, 8, 8);
			dist[1] = image_pixel_dist(patch_image, image, patch_y, patch_x, left_new_match2_y, left_new_match2_x, 8, 8);
			dist[2] = image_pixel_dist(patch_image, image, patch_y, patch_x, up_new_match1_y, up_new_match1_x, 8, 8);
			dist[3] = image_pixel_dist(patch_image, image, patch_y, patch_x, up_new_match2_y, up_new_match2_x, 8, 8);
			int new_x, new_y;
			if (dist[0]<=dist[1]&&dist[0]<=dist[2]&&dist[0]<=dist[3])
			{
				new_x = left_new_match1_x;
				new_y = left_new_match1_y;
			}
			else if (dist[1]<=dist[0]&&dist[1]<=dist[2]&&dist[0]<=dist[3])
			{
				new_x = left_new_match2_x;
				new_y = left_new_match2_y;
			}
			else if (dist[2]<=dist[0]&&dist[2]<=dist[1]&&dist[2]<=dist[3])
			{
				new_x = up_new_match1_x;
				new_y = up_new_match1_y;
			}
			else
			{
				new_x = up_new_match2_x;
				new_y = up_new_match2_y;
			}
			return feature_coor[image_col*new_y + new_x]->p_m_leaf;
		}*/
	}
	else if (!up_feature&&left_feature)
	{
		/*int patch_x = left_feature->x + 1;
		int patch_y = left_feature->y;
		int left_new_match1_x = left_feature->fwd_match[0].x + 1;
		int left_new_match1_y = left_feature->fwd_match[0].y;
		int left_new_match2_x = left_feature->fwd_match[1].x + 1;
		int left_new_match2_y = left_feature->fwd_match[1].y;
		if (left_new_match1_x < image_col && left_new_match2_x < image_col)
		{
			int dist[2];
			dist[0] = image_pixel_dist(patch_image, image, patch_y, patch_x, left_new_match1_y, left_new_match1_x, 8, 8);
			dist[1] = image_pixel_dist(patch_image, image, patch_y, patch_x, left_new_match2_y, left_new_match2_x, 8, 8);
			int new_x, new_y;
			if (dist[0]<dist[1])
			{
				new_x = left_new_match1_x;
				new_y = left_new_match1_y;
			}
			else
			{
				new_x = left_new_match2_x;
				new_y = left_new_match2_y;
			}
			//return feature_coor[image_col*new_y + new_x]->p_m_leaf;
			struct feature*f = feature_coor[image_col*new_y + new_x];
			return f->p_m_leaf;
		}
		else
		{
			return NULL;
		}*/
		int patch_x = left_feature->x + 1;
		int patch_y = left_feature->y;
		int patch_num =left_feature->fwd_match_num;
		int *new_x = (int *)malloc(patch_num*sizeof(int)), *new_y = (int *)malloc(patch_num*sizeof(int));
		for (int i = 0; i < left_feature->fwd_match_num; i++)
		{
			new_x[i] = (*(left_feature->fwd_match+i))->x + 1;
			new_y[i] = (*(left_feature->fwd_match + i))->y;
		}
		int *dist = (int *)malloc(patch_num*sizeof(int));
		struct feature**feature_array = (struct feature**)malloc(patch_num*sizeof(struct feature*));
		int leaf_num = 0;//与node不同的m_leaf的数目;
		for (int i = 0; i < patch_num; i++)
		{
			int num = new_y[i] * image_col + new_x[i];
			if (new_x[i] < image_col&&new_y[i] < image_row - 7 && feature_coor[num] && feature_coor[num]->p_m_leaf != node)
			{
				dist[leaf_num] = image_pixel_dist(patch_image, image, patch_y, patch_x, new_y[i], new_x[i], 8, 8);
				feature_array[leaf_num] = feature_coor[num];
				leaf_num++;
			}
		}
		if (leaf_num)
		{
			int index = 0, min_dist = dist[0];
			for (int i = 1; i < leaf_num; i++)
			{
				if (dist[i] <= min_dist)
				{
					min_dist = dist[i];
					index = i;
				}
			}
			return feature_array[index]->p_m_leaf;
		}
		else
		{
			return NULL;
		}
		free(new_x);
		free(new_y);
		free(dist);
		free(feature_array);
	}
	else
	{
		int patch_x = up_feature->x;
		int patch_y = up_feature->y + 1;
		int patch_num = up_feature->fwd_match_num;
		int *new_x = (int *)malloc(patch_num*sizeof(int)), *new_y = (int *)malloc(patch_num*sizeof(int));
		for (int i = 0; i < up_feature->fwd_match_num; i++)
		{
			new_x[i] = (*(up_feature->fwd_match+i))->x;
			new_y[i] = (*(up_feature->fwd_match + i))->y + 1;
		}
		int *dist = (int *)malloc(patch_num*sizeof(int));
		struct feature**feature_array = (struct feature**)malloc(patch_num*sizeof(struct feature*));
		int leaf_num = 0;//与node不同的m_leaf的数目;
		for (int i = 0; i < patch_num; i++)
		{
			int num = new_y[i] * image_col + new_x[i];
			if (new_x[i] < image_col&&new_y[i] < image_row - 7 && feature_coor[num] && feature_coor[num]->p_m_leaf != node)
			{
				dist[leaf_num] = image_pixel_dist(patch_image, image, patch_y, patch_x, new_y[i], new_x[i], 8, 8);
				feature_array[leaf_num] = feature_coor[num];
				leaf_num++;
			}
		}
		if (leaf_num)
		{
			int index = 0, min_dist = dist[0];
			for (int i = 1; i < leaf_num; i++)
			{
				if (dist[i] <= min_dist)
				{
					min_dist = dist[i];
					index = i;
				}
			}
			return feature_array[index]->p_m_leaf;
		}
		else
		{
			return NULL;
		}
		free(new_x);
		free(new_y);
		free(dist);
		free(feature_array);
		/*if (up_new_match1_y < image_row - 7 && up_new_match2_y < image_row - 7)
		{
			int dist[2];
			dist[0] = image_pixel_dist(patch_image, image, patch_y, patch_x, up_new_match1_y, up_new_match1_x, 8, 8);
			dist[1] = image_pixel_dist(patch_image, image, patch_y, patch_x, up_new_match2_y, up_new_match2_x, 8, 8);
			int new_x, new_y;
			if (dist[0] <= dist[1])
			{
				new_x = up_new_match1_x;
				new_y = up_new_match1_y;
			}
			else
			{
				new_x = up_new_match2_x;
				new_y = up_new_match2_y;
			}
			return feature_coor[image_col*new_y + new_x]->p_m_leaf;
		}
		else
		{
			return NULL;
		}*/
	}
}


/*
find an image feature's 最邻近值在kdtree中使用BBF算法
@param kd_root root of an image feature kd tree
@param feat image feature for whose neighbors to search
@param leaf1 the pointer to the kd_node which is propogated from nearby matches
@param k number of neighbors to find
@param nbrs pointer to an array in which to store pointers to neighbors
in order of increasing descriptor distance
@param max_nn_chks search is cut off after examining this many tree entries*/
int kdtree_bbf_knn(struct kd_node*node, struct feature*feat, struct kd_node*leaf1, int k, struct feature***nbrs, int max_nn_chks,
	int distant_thresh_2)
{
	
	struct min_pq*min_pq = minpq_init();
	int i, t = 0, n = 0;
	struct kd_node*expl;
	struct feature* tree_feat = NULL, ** _nbrs = (struct feature**)calloc(k, sizeof(struct feature*));
	struct bbf_data* bbf_data;
	if (leaf1 == NULL)
	{
		minpq_insert(min_pq, node, 0);
		while (min_pq->n > 0 && t < max_nn_chks)
		{
			expl = (struct kd_node*)minpq_extract_min(min_pq);//提取优先级最高的节点,复制给当前节点;
			if (!expl)
			{
				fprintf(stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
					__FILE__, __LINE__);
				goto fail;
			}

			expl = explore_to_leaf(expl, feat, min_pq);//从当前搜索节点expl一直搜索到叶子节点，
			//搜索过程中将未搜索的节点根据优先级放入队列，返回值为叶子节点  
			if (!expl)
			{
				fprintf(stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
					__FILE__, __LINE__);
				goto fail;
			}

			//比较查找最近邻;
			for (i = 0; i < expl->n; i++)
			{
				tree_feat = &expl->features[i];
				if ((feat->x-tree_feat->x)*(feat->x-tree_feat->x)+(feat->y-tree_feat->y)*(feat->y-tree_feat->y)<=distant_thresh_2)
				{
					continue;
				}
				bbf_data = (struct bbf_data*)malloc(sizeof(struct bbf_data));
				if (!bbf_data)
				{
					fprintf(stderr, "Warning: unable to allocate memory,"
						" %s line %d\n", __FILE__, __LINE__);
					goto fail;
				}
				bbf_data->old_data = tree_feat->feature_data;//保存第i个特征点的feature_data域以前的值;
				bbf_data->d = descr_dist_sq(feat, tree_feat);//当前搜索点和目标点之间的欧氏距离;
				tree_feat->feature_data = bbf_data;
				n += insert_into_nbr_array(tree_feat, _nbrs, n, k);//当最近邻数组中元素已经达到k时，继续插入元素个数不会增加，但是会更新元素值;
			}
			t++;
		}

		minpq_release(&min_pq,feat);
		for (i = 0; i < n; i++)
		{
			bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
			_nbrs[i]->feature_data = bbf_data->old_data;
			free(bbf_data);
		}
		*nbrs = _nbrs;
		return n;
	}
	else
	{
		minpq_insert(min_pq, node, 0);
		int ki = leaf1->ki;
		int kv = leaf1->kv;
		minpq_insert(min_pq, leaf1, ABS(feat->descr[ki] - kv));
		while (min_pq->n > 0 && t < max_nn_chks)
		{
			expl = (struct kd_node*)minpq_extract_min(min_pq);//提取优先级最高的节点,复制给当前节点;
			if (!expl)
			{
				fprintf(stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
					__FILE__, __LINE__);
				goto fail;
			}

			expl = explore_to_leaf(expl, feat, min_pq);//从当前搜索节点expl一直搜索到叶子节点，
			//搜索过程中将未搜索的节点根据优先级放入队列，返回值为叶子节点  
			if (!expl)
			{
				fprintf(stderr, "Warning: PQ unexpectedly empty, %s line %d\n",
					__FILE__, __LINE__);
				goto fail;
			}

			//比较查找最近邻;
			for (i = 0; i < expl->n; i++)
			{
				tree_feat = &expl->features[i];
				if ((feat->x - tree_feat->x)*(feat->x - tree_feat->x) + (feat->y - tree_feat->y)*(feat->y - tree_feat->y) <= distant_thresh_2)
				{
					continue;
				}
				bbf_data = (struct bbf_data*)malloc(sizeof(struct bbf_data));
				if (!bbf_data)
				{
					fprintf(stderr, "Warning: unable to allocate memory,"
						" %s line %d\n", __FILE__, __LINE__);
					goto fail;
				}
				bbf_data->old_data = tree_feat->feature_data;//保存第i个特征点的feature_data域以前的值;
				bbf_data->d = descr_dist_sq(feat, tree_feat);//当前搜索点和目标点之间的欧氏距离;
				tree_feat->feature_data = bbf_data;
				n += insert_into_nbr_array(tree_feat, _nbrs, n, k);//当最近邻数组中元素已经达到k时，继续插入元素个数不会增加，但是会更新元素值;
			}
			t++;
		}

		minpq_release(&min_pq,feat);
		for (i = 0; i < n; i++)
		{
			bbf_data = (struct bbf_data*)_nbrs[i]->feature_data;
			_nbrs[i]->feature_data = bbf_data->old_data;
			free(bbf_data);
		}
		*nbrs = _nbrs;
		return n;

	fail:
		minpq_release(&min_pq,feat);
		for (i = 0; i < n; i++)
		{
			bbf_data = (struct bbf_data*) _nbrs[i]->feature_data;
			_nbrs[i]->feature_data = bbf_data->old_data;
			free(bbf_data);
		}
		free(_nbrs);
		*nbrs = NULL;
		return -1;
	}
}

bool valid(int left, int up, struct mask*mask)
{
	int right = left + 7;
	int down = up + 7;
	for (int i = up; i <= down;i++)
	{
		if (*(mask->lines+i)==NULL)
		{
			continue;
		}
		if (*(mask->min+i)>right||*(mask->max+i)<left)
		{
			continue;
		}
		for (int j = 0; j < *(mask->n + i);j++)
		{
			if (left >= (*(mask->lines + i)+j)->left&&left <= (*(mask->lines + i)+j)->right
				|| (left <= (*(mask->lines + i) + j)->left&&right >= (*(mask->lines + i) + j)->left))
			{
				return false;
			}
		}
	}
	return true;
}

//在-width to width和-height to height内的offset均被统计在内;
struct offset_num** dominant_offset_form(int width, int height,struct feature**feature_coor,int col,int row)
{
	/*int **digram = (int **)malloc((2*height+1)*sizeof(int *));*/
	struct offset_num**digram = (struct offset_num**)malloc((2 * height + 1)*(2 * width + 1)*sizeof(struct offset_num*));

	for (int i = 0; i < (2 * width + 1)*(2 * height + 1); i++)
	{
		*(digram + i) = (struct offset_num *)malloc(sizeof(struct offset_num));
	}
	for (int i = 0; i < 2 * height + 1;i++)
	{
		for (int j = 0; j < 2 * width + 1;j++)
		{
			(*(digram + i*(2 * width + 1) + j))->offset_x = j - width;
			(*(digram + i*(2 * width + 1) + j))->offset_y = i - height;
			(*(digram + i*(2 * width + 1) + j))->num = 0;
		}
	}
	for (int i = 0; i < row;i++)
	{
		for (int j = 0; j < col; j++)
		{
			int num = i*col + j;
			if (feature_coor[num]==NULL||feature_coor[num]->fwd_match_num==0)
			{
				continue;
			}
			else
			{
				int offset_x = (*(feature_coor[num]->fwd_match))->x - feature_coor[num]->x;
				int offset_y = (*(feature_coor[num]->fwd_match))->y - feature_coor[num]->y;
				if (offset_x>=-width&&offset_x<=width&&offset_y>=-height&&offset_y<=height)
				{
					(*(digram + (offset_y + height)*(2 * width + 1) + (offset_x + width)))->num++;
				}
			}
		}
	}
	return digram;
}



//n为offset_num中的总数量;
struct offset_num**sort(struct offset_num**offset_num,int k,int n)
{
	struct offset_num**k_offsets = (struct offset_num**)malloc(k*sizeof(struct offset_num*));
	int num = 0;//k_offsets中现有的数量;
	for (int i = 0; i < n;i++)
	{
		if (num==k)
		{
			if ((*(offset_num + i))->num <= (*(k_offsets + k-1))->num)
			{
				continue;
			}
			else
			{
				int place = dichotomy_find_place(*(offset_num + i), k_offsets, 0, k - 1);
				for (int j = k - 1; j > place;j--)
				{
					*(k_offsets + j) = (*(k_offsets + j - 1));
				}
				*(k_offsets + place) = (*(offset_num + i));
			}
		}
		else
		{
			if (num==0)
			{
				*(k_offsets) = (*(offset_num + i));
				num = 1;
			}
			else if ((*(offset_num+i))->num<=(*(k_offsets+num-1))->num)
			{
				*(k_offsets + num) = (*(offset_num + i));
				num++;
			}
			else
			{
				int place = dichotomy_find_place(*(offset_num + i), k_offsets, 0, num - 1);
				for (int j = num; j > place;j--)
				{
					(*(k_offsets + j)) = (*(k_offsets + j - 1));
				}
				*(k_offsets + place) = (*(offset_num + i));
				num++;
			}
		}
	}
	return k_offsets;
}

int dichotomy_find_place(struct offset_num*offset, struct offset_num**k_offset, int first_one, int last_one)
{
	int num = offset->num;
	while (first_one<last_one)
	{
		if (last_one-first_one==1)
		{
			return num>(*(k_offset + first_one))->num ? first_one : last_one;
		}
		int mid = (first_one + last_one) / 2;
		int mid_num = (*(k_offset + mid))->num;
		if (num>mid_num)
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

mxArray*output_diagram(struct offset_num**digram, int width, int height)
{
	int col = 2 * width + 1;
	int row = 2 * height + 1;
	mxArray* out = mxCreateDoubleMatrix(row, col, mxREAL);
	double *data = (double *)mxGetData(out);
	for (int j = 0; j < col; j++)
	{
		for (int i = 0; i < row; i++)
		{
			*(data + j*row + i) = (*(digram + i*col + j))->num;
		}
	}
	return out;
}

int array_double(void** array, int n, int size)
{
	void* tmp;

	tmp = realloc(*array, 2 * n * size);
	if (!tmp)
	{
		fprintf(stderr, "Warning: unable to allocate memory in array_double(),"
			" %s line %d\n", __FILE__, __LINE__);
		if (*array)
			free(*array);
		*array = NULL;
		return 0;
	}
	*array = tmp;
	return n * 2;
}