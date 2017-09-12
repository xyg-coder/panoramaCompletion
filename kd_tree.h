#include "mex.h"

struct feature
{
	int x;
	int y;
	int d;//descriptor length;
	int *descr;
	struct feature**fwd_match;//matching feature from forward image;
	struct kd_node* p_leaf;//此feature最终对应的leaf的指针;
	struct kd_node* p_m_leaf;//此feature对应的m_leaf;
	void *feature_data;
	bool valid;
	int fwd_match_num;//匹配出来的fwd_match数目;
};

struct kd_node
{
	int ki;//partition key index;
	int kv;//partition key value;
	int leaf;//1 if node is a leaf,0 otherwise;
	int m_leaf;//是含有features m(8~64)个的node;
	struct feature* features;
	int n;//number of features;
	struct kd_node* kd_left;
	struct kd_node* kd_right;
};

/** an element in a minimizing priority queue */
struct pq_node
{
	void* data;
	int key;
};


/** a minimizing priority queue */
struct min_pq
{
	struct pq_node* pq_array;    /* array containing priority queue */
	int nallocd;                 /* number of elements allocated */
	int n;                       /**< number of elements in pq */
};

struct bbf_data
{
	int d;
	void *old_data;
};

struct line_mask
{
	int left;
	int right;
};

struct mask
{
	struct line_mask**lines;//对应每一行line_mask数组;
	int *n;//对应每一行的mask段的数目;
	int *min;//每一行line_mask最小值;
	int *max;//每一行line_mask最大值;
	int row;//应该等于图像大小;
};

struct offset_num
{
	int offset_x;
	int offset_y;
	int num;
};


int patch_match(struct feature* image_features,unsigned int**image,
	int image_cols, int image_rows, int n,struct mask*mask,mxArray*out[]);
struct kd_node* kdtree_build(struct feature* features, int*n, struct feature**feature_coor, int image_col, struct mask*mask);
static struct kd_node* kd_node_init(struct feature* features, int n);
static void expand_kd_node_subtree(struct kd_node* kd_node,struct feature**feature_coor,int image_col,int mode);
static void assign_part_key(struct kd_node* kd_node);
static double median_select(int* array, int n);
static double rank_select(int* array, int n, int r);
static void insertion_sort(int* array, int n);
static int* partition_array(int* array, int n, int pivot);
static void partition_features(struct kd_node* kd_node,struct feature**feature_coor,int image_col);
int kdtree_bbf_knn_for_patch1(struct kd_node* kd_root, struct feature* feat, int k,
struct feature*** nbrs, int max_nn_chks, int distant_thresh_2);
void kdtree_release(struct kd_node* kd_root);
struct min_pq* minpq_init();
int minpq_insert(struct min_pq* min_pq, void* data, int key);
static void decrease_pq_node_key(struct pq_node* pq_array, int i, int key);
void* minpq_extract_min(struct min_pq* min_pq);
static void restore_minpq_order(struct pq_node* pq_array, int i, int n);
static struct kd_node* explore_to_leaf(struct kd_node* kd_node,struct feature* feat,
struct min_pq* min_pq);
int descr_dist_sq(struct feature* f1, struct feature* f2);
static int insert_into_nbr_array(struct feature* feat, struct feature** nbrs,
	int n, int k);
void minpq_release(struct min_pq** min_pq,struct feature*feat);
int image_pixel_dist(unsigned int**patch_pixel, unsigned int**image_pixel,
	int patch_y, int patch_x, int image_y, int image_x, int width, int height);
struct kd_node*propogate(struct kd_node*node,struct feature*features, struct feature*up_feature, struct feature*left_feature,
	unsigned int**patch_image, unsigned int** image,
	int image_row, int image_col, struct feature**feature_coor,struct mask*mask);
int kdtree_bbf_knn(struct kd_node*kd_root, struct feature*feat, struct kd_node*leaf1, int k, struct feature***nbrs, int max_nn_chks,
	int distant_thresh_2);
void feature_init(struct feature**feat, const mxArray*feat_array, int dim, int n, int col, int row);
void image_init(unsigned char*image, unsigned int***new_image, int row, int col);
bool valid(int left, int up, struct mask*mask);
struct offset_num** dominant_offset_form(int width, int height, struct feature**feature_coor, int col, int row);
struct offset_num** sort(struct offset_num**offset_num, int k,int n);
int dichotomy_find_place(struct offset_num*offset, struct offset_num**k_offset, int first_one, int last_one);
void mask_init(int *mask_matrix, struct mask**mask, int row, int col);
void mask_release(struct mask**mask, int row);
mxArray*output_diagram(struct offset_num**digram, int width, int height);
int array_double(void** array, int n, int size);