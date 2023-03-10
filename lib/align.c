#include "align.h"

unsigned char seq_comp_tb[256] = { 
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, 
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31, 
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47, 
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63, 
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95, 
     96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255 
};

void malloc_vd(int **V, uint8_t ***D, uint64_t max_mem_d){
	*V = (int *) malloc( max_mem_d * 2 * sizeof(int));
	int i;
	*D = NULL;
	for (i = 0; i < 3 && ! *D; i++){
		*D = (uint8_t **) malloc( max_mem_d * sizeof(uint8_t *) + max_mem_d * (max_mem_d + 1)/2 * sizeof(uint8_t));
	}
	if (*D) {
		uint64_t d;
		uint8_t * const _ = (uint8_t *) (*D + max_mem_d);
		for (d = 0; d < max_mem_d; d ++ ) {
			(*D)[d] = _ + d * (d + 1)/2;
		}
	}
}

void clean_V(int *V, int max_mem_d){
	memset(V, 0, max_mem_d * 2 * sizeof(int));
}

void destory_vd(int *V, uint8_t **D){
	free(V);
	free(D);
}

void reverse_str(char *str, int len)
{
	char tmp;
	char *p1 = str;
	char *p2 = str + len - 1;
	while (p1 < p2) {
		tmp = *p1;
		*p1++ = *p2;
		*p2-- = tmp;
	}
}

void revcomp_bseq(char *str, int len)
{
    char tmp;
    char *p1 = str;
    char *p2 = str + len - 1;
    while (p1 < p2) {
        tmp = *p1;
        *p1++ = seq_comp_tb[(uint8_t) *p2];
        *p2-- = seq_comp_tb[(uint8_t) tmp];
    }
    if (p1 == p2) *p1 = seq_comp_tb[(uint8_t) *p1];
}

void str_tolower(char *p) {
	for ( ; *p; ++p) *p |=0x20;
}

void str_toupper(char *p) {
	for ( ; *p; ++p) *p &= 0xdf;//~0x20;
}

void ide(const char *query_seq, int q_len, const char *target_seq, int t_len,
	int *V, uint8_t **D, int max_d, int band_size, int *mlen, int *blen) {

	int x, new_min_k, new_max_k, k, k2, d;
	int y = x = 0;
	int kk = 0;
	int min_k = 0;
	int max_k = 0;
	int best_m = -1;
	int k_offset = max_d;

	for (d = 0; d < max_d && max_k - min_k <= band_size; d ++ ) {
		for (k = min_k; k <= max_k;  k += 2) {
			kk = k < 0 ? -1 * k - 1 : k;
			if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
				x = V[ k + 1 + k_offset];
				D[d][kk] = 0;
			} else {
				x = V[ k - 1 + k_offset] + 1;
				D[d][kk] = 1;
			}


			y = x - k;
			while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
				x++;
				y++;
			}

			V[ k + k_offset] = x;

			if ( x + y > best_m) {
				best_m = x + y;
			}

			if ( x >= q_len || y >= t_len) {//for prefix aln
				*mlen = x - (k + d)/2;
				*blen = y + (k + d)/2;
				return;
			// printf("x:%d y:%d d:%d k:%d m:%d l:%d\n",x, y, d, k, x - (k + d)/2,y + (k + d)/2 );
			}
		}
		// For banding
		new_min_k = max_k;
		new_max_k = min_k;

		k2 = min_k;
		while (k2 < new_min_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_min_k = k2;
			k2 += 2;
		}
		
		k2 = max_k;
		while (k2 > new_max_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_max_k = k2;
			k2 -= 2;
		}

		max_k = new_max_k + 1;
		min_k = new_min_k - 1;
	}
}

//For trimed align
//atcgat, atcatcgatcatg
//---atcgat----, atcatcgatcatg, ide=100%
void alnpos(const char *query_seq, int q_len, const char *target_seq, int t_len,
    int *V, uint8_t **D, int max_d, int band_size, alignpos *align_rtn){

	int x, new_min_k, new_max_k, k, k2, d;
	int y = x = 0;
	int kk = 0;
	int min_k = 0;
	int max_k = 0;
	int best_m = -1;
	int k_offset = max_d;
	bool aligned = false;

	for (d = 0; d < max_d && max_k - min_k <= band_size; d ++ ) {
		for (k = min_k; k <= max_k;  k += 2) {
			kk = k < 0 ? -1 * k - 1 : k;
			if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
				x = V[ k + 1 + k_offset];
				D[d][kk] = 0;
			} else {
				x = V[ k - 1 + k_offset] + 1;
				D[d][kk] = 1;
			}


			y = x - k;
			while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
				x++;
				y++;
			}

			V[ k + k_offset] = x;

			if ( x + y > best_m) {
				best_m = x + y;
			}

			if ( x >= q_len || y >= t_len) {//for prefix aln
				aligned = true;
				break;
				// *mlen = x - (k + d)/2;
				// *blen = y + (k + d)/2;
				// return;
			// printf("x:%d y:%d d:%d k:%d m:%d l:%d\n",x, y, d, k, x - (k + d)/2,y + (k + d)/2 );
			}
		}

		if (!aligned){
			// For banding
			new_min_k = max_k;
			new_max_k = min_k;

			k2 = min_k;
			while (k2 < new_min_k ){
				if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_min_k = k2;
				k2 += 2;
			}
			
			k2 = max_k;
			while (k2 > new_max_k ){
				if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_max_k = k2;
				k2 -= 2;
			}

			max_k = new_max_k + 1;
			min_k = new_min_k - 1;
		}else{
			align_rtn->aln_q_e = x;
			align_rtn->aln_t_e = y;
			x --;//0-based

			int gap, aln_pos, pre_d, pre_k, pre_kk, pre_x;
			gap = aln_pos = 0;

			while (true){
				while (x >= 0 && x >= k && query_seq[x] == target_seq[x - k]){
					x --;
					aln_pos ++;
				}

				pre_d = d - 1;
				if (x < 0 || x - k < 0) break;
				if (x < k || (x >= 0 && D[d][kk])){
				// if (D[d][kk]){
					pre_k = k - 1;
					pre_x = x - 1;
				}else{
					pre_k = k + 1;
					pre_x = x;
				}

				pre_kk = pre_k < 0 ? -1 * pre_k - 1 : pre_k;
				aln_pos ++;
				gap ++;
				d = pre_d;
				k = pre_k;
				kk = pre_kk;
				x = pre_x;
			}
			align_rtn->aln_len = aln_pos;
			align_rtn->aln_mlen = aln_pos - gap;
			align_rtn->aln_q_s = x + 1;
			align_rtn->aln_t_s = x + 1 - k;
			break;
		}
	}
}

//from 5->3
void extend_fwd(const char *query_seq, int q_len, const char *target_seq, int t_len,
	int *V, uint8_t **D, int max_d, int band_size, float d_factor, int *bstx, int *bsty) {

	int x, new_min_k, new_max_k, k, k2, d;
	int y = x = 0;
	int kk = 0;
	int min_k = 0;
	int max_k = 0;
	int best_m = -1;
	int aligned = 0;
	int k_offset = max_d;
	int best_x, best_y, best_d, global_best_m = -1;
	float score, peak_score = 0;
	*bstx = *bsty = 0;
	
	for (d = 0; d < max_d && max_k - min_k <= band_size && !aligned; d ++ ) {
		for (k = min_k; k <= max_k;  k += 2) {
			kk = k < 0 ? -1 * k - 1 : k;
			if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
				x = V[ k + 1 + k_offset];
				D[d][kk] = 0;
			} else {
				x = V[ k - 1 + k_offset] + 1;
				D[d][kk] = 1;
			}

			y = x - k;
			while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
				x++;
				y++;
			}

			V[ k + k_offset] = x;

			if ( x + y > best_m) {
				best_m = x + y;
				if (best_m > global_best_m){
					best_x = x;
					best_y = y;
					best_d = d;
					// best_kk = kk;
					global_best_m = best_m;
					score = (best_x + best_y) * d_factor - best_d;
					if (score > peak_score){//find the global peak score
						peak_score = score;
						*bstx = best_x;
						*bsty = best_y;
					}else if (score < peak_score - 30) return;//TODO check
				}
			}
			

			if ( x >= q_len || y >= t_len) {//for prefix aln
				aligned = 1;
				best_x = x;
				best_y = y;
				best_d = d;
				score = (best_x + best_y) * d_factor - best_d;
				if (score > 0){//peak_score){
					peak_score = score;
					*bstx = best_x;
					*bsty = best_y;
				}
				return;
			}
		}
		// For banding
		new_min_k = max_k;
		new_max_k = min_k;

		k2 = min_k;
		while (k2 < new_min_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_min_k = k2;
			k2 += 2;
		}
		
		k2 = max_k;
		while (k2 > new_max_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_max_k = k2;
			k2 -= 2;
		}
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;
	}
}

//from 3->5
void extend_rev(const char *query_seq, int q_len, const char *target_seq, int t_len,
	int *V, uint8_t **D, int max_d, int band_size, float d_factor, int *bstx, int *bsty) {

	int x, new_min_k, new_max_k, k, k2, d;
	int y = x = 0;
	int kk = 0;
	int min_k = 0;
	int max_k = 0;
	int best_m = -1;
	int aligned = 0;
	int k_offset = max_d;
	int best_x, best_y, best_d, global_best_m = -1;
	float score, peak_score = 0;
	*bstx = *bsty = 0;
	
	for (d = 0; d < max_d && max_k - min_k <= band_size && !aligned; d ++ ) {
		for (k = min_k; k <= max_k;  k += 2) {
			kk = k < 0 ? -1 * k - 1 : k;
			if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
				x = V[ k + 1 + k_offset];
				D[d][kk] = 0;
			} else {
				x = V[ k - 1 + k_offset] + 1;
				D[d][kk] = 1;
			}

			y = x - k;
			while ( x < q_len && y < t_len && query_seq[q_len - x - 1] == target_seq[t_len - y - 1] ){
				x++;
				y++;
			}

			V[ k + k_offset] = x;

			if ( x + y > best_m) {
				best_m = x + y;
				if (best_m > global_best_m){
					best_x = x;
					best_y = y;
					best_d = d;
					// best_kk = kk;
					global_best_m = best_m;
					score = (best_x + best_y) * d_factor - best_d;
					if (score > peak_score){//find the global peak score
						peak_score = score;
						*bstx = best_x;
						*bsty = best_y;
					}else if (score < peak_score - 30) return;//TODO check
				}
			}
			

			if ( x >= q_len || y >= t_len) {//for prefix aln
				aligned = 1;
				best_x = x;
				best_y = y;
				best_d = d;
				// best_kk = kk;
				score = (best_x + best_y) * d_factor - best_d;
				if (score > 0){//peak_score){
					peak_score = score;
					*bstx = best_x;
					*bsty = best_y;
				}
				return;
			}
		}
		// For banding
		new_min_k = max_k;
		new_max_k = min_k;

		k2 = min_k;
		while (k2 < new_min_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_min_k = k2;
			k2 += 2;
		}
		
		k2 = max_k;
		while (k2 > new_max_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_max_k = k2;
			k2 -= 2;
		}
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;
	}
}

void core(char *query_seq, int q_len, char *target_seq, int t_len,
	alignment *align_rtn, int *V, uint8_t **D, int max_d, int band_size) {

	int x, new_min_k, new_max_k, k, k2, d;
	int y = x = 0;
	int kk = 0;
	int min_k = 0;
	int max_k = 0;
	int best_m = -1;
	bool aligned = false;
	int k_offset = max_d;

	for (d = 0; d < max_d && max_k - min_k <= band_size; d ++ ) {
		for (k = min_k; k <= max_k;  k += 2) {
			kk = k < 0 ? -1 * k - 1 : k;
			if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
				x = V[ k + 1 + k_offset];
				D[d][kk] = 0;
			} else {
				x = V[ k - 1 + k_offset] + 1;
				D[d][kk] = 1;
			}

			y = x - k;
			while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
				x++;
				y++;
			}

			V[ k + k_offset] = x;

			if ( x + y > best_m) {
				best_m = x + y;
			}

			// if ( x >= q_len || y >= t_len) {//for prefix aln
			// 	aligned = true;
			// 	break;
			// }
			if ( x >= q_len && y >= t_len) { //for global aln
				aligned = true;
				break;
			}
		}
		// For banding
		new_min_k = max_k;
		new_max_k = min_k;

		k2 = min_k;
		while (k2 < new_min_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_min_k = k2;
			k2 += 2;
		}
		
		k2 = max_k;
		while (k2 > new_max_k ){
			if (V[k2 + k_offset] * 2 - k2 >= best_m - 150) new_max_k = k2;
			k2 -= 2;
		}

		max_k = new_max_k + 1;
		min_k = new_min_k - 1;

		if (aligned == true){
			x--; // 0-based
			max_d = d;
			align_rtn->aln_t_e = align_rtn->aln_t_s + y - 1;// fix bug for tail with clip
			align_rtn->aln_t_len = y;
			align_rtn->aln_q_len = x + 1;

			int gap, aln_pos, pre_d, pre_k, pre_kk, pre_x, pre_y;
			gap = aln_pos = 0;

			while (true){
				while (x >= 0 && x >= k && query_seq[x] == target_seq[x - k]){
					align_rtn->t_aln_str[aln_pos] = align_rtn->q_aln_str[aln_pos] = query_seq[x];
					x --;
					aln_pos ++;
					gap = 0;
				}

				pre_d = d - 1;

				if (x < 0 && x - k < 0) break;
				if (x < k || (x >= 0 && D[d][kk])){
				// if (D[d][kk]){
					pre_k = k - 1;
					pre_x = x - 1;
				}else{
					pre_k = k + 1;
					pre_x = x;
				}

				pre_y = pre_x - pre_k;
				pre_kk = pre_k < 0 ? -1 * pre_k - 1 : pre_k;

				if (pre_x == x && pre_y !=  x - k){ //advance in y

					if (x - k < 0) gap = 260;
					else{
						align_rtn->q_aln_str[aln_pos] = '-';
						align_rtn->t_aln_str[aln_pos++] = target_seq[x - k];
					}
					

				}else{ //advance in x
					
					if (x < 0) gap = 260;
					else{
						align_rtn->q_aln_str[aln_pos] = query_seq[x];
						align_rtn->t_aln_str[aln_pos ++] = '-';
					}
				}
				
				if (gap ++ > 250){//only allow the max length of a gap = 250
					aln_pos = 2;
					break;
				}

				d = pre_d;
				k = pre_k;
				kk = pre_kk;
				x = pre_x;
			}
			align_rtn->aln_len = aln_pos;
			reverse_str(align_rtn->t_aln_str, align_rtn->aln_len);
			reverse_str(align_rtn->q_aln_str, align_rtn->aln_len);
			align_rtn->t_aln_str[align_rtn->aln_len] = '\0';
			align_rtn->q_aln_str[align_rtn->aln_len] = '\0';
			break;
		}

	}
}

void align_hq(char *query_seq, int q_len,
				char *target_seq, int t_len,
				alignment *align_rtn, int *V, uint8_t **D) {
	
	int max_d = (int) ((q_len + t_len > 1000 ? 0.1 : 0.5) * (q_len + t_len));
	int band_size = (int) ((q_len + t_len > 1000 ? 0.03 : 0.3) * (q_len + t_len));
	core(query_seq, q_len, target_seq, t_len, align_rtn, V, D, max_d, band_size);
}

void align(char *query_seq, int q_len,
				char *target_seq, int t_len,
				alignment *align_rtn, int *V, uint8_t **D) {
	int max_d = (int) (0.4 * (q_len + t_len));
	int band_size =  (int) ((q_len + t_len > 5000 ? 0.1 : 1) * (q_len + t_len));
	core(query_seq, q_len, target_seq, t_len, align_rtn, V, D, max_d, band_size);
}

/////////////////////NM ALIGN/////////////////
static uint8_t MMH[] = {64, 16, 4, 1};
static uint8_t INS[] = {128, 32, 8, 2};
static uint8_t DEL[] = {192, 48, 12, 3};

typedef struct {
	int mas; //2
	int mis; //-4
	int gos; //-4
	int ges; //-2
} sco;

void align_nd_tb(const char *s1, const char *s2, uint32_t s1_l, uint32_t s2_l, uint8_t**d, alignment *aln){
	uint8_t i, j;
	aln->aln_len = 0;
	while(s1_l != 0 || s2_l != 0){
		i = s2_l&3;
		j = d[s1_l][s2_l>>2] & DEL[i];

		// printf("%d %d %d %d\n", s1_l, s2_l, i, aln->aln_len);
		if (j == MMH[i]){
				aln->t_aln_str[aln->aln_len] = s1[--s1_l];
				aln->q_aln_str[aln->aln_len] = s2[--s2_l];
				aln->aln_len ++;
		}else if (j == DEL[i]){
			aln->q_aln_str[aln->aln_len] = s2[--s2_l];
			aln->t_aln_str[aln->aln_len] = '-';
			aln->aln_len ++;
		}else{
			aln->t_aln_str[aln->aln_len] = s1[--s1_l];
			aln->q_aln_str[aln->aln_len] = '-';
			aln->aln_len ++;
		}
	}
	aln->t_aln_str[aln->aln_len] = aln->q_aln_str[aln->aln_len] = '\0';
	reverse_str(aln->t_aln_str, aln->aln_len);
	reverse_str(aln->q_aln_str, aln->aln_len);
}

void align_nd_dp(const char *s1, const char *s2, const uint32_t s1_l, const uint32_t s2_l, uint8_t**d){
	uint32_t i, j;
	sco score = {
		.mas = 2,
		.mis = -4,
		.gos = -4,
		.ges = -2
	};

	int32_t cs, ms, is, ds, *s = (int32_t *) malloc((s2_l + 1) * sizeof(int32_t));
	s[0] = cs = 0;
	d[0][0] |= MMH[0];
	for (j = 1; j < s2_l + 1; j++){
		s[j] = s[j - 1] + ((d[0][(j - 1)>>2] & DEL[(j - 1)&3]) == DEL[(j - 1)&3] ? score.ges : score.gos);
		d[0][j>>2] |= DEL[j&3];
	}
	
	register int k;
	for (i = 1; i < s1_l + 1; i++){
		for (j = 0; j < s2_l + 1; j++){
			if (j == 0){
				cs = s[j] + ((d[i - 1][0]& DEL[0]) == INS[0] ? score.ges : score.gos);
				d[i][j] = INS[0];
			}else{
				k = j&3;
				ms = s[j - 1] + (s1[i - 1] == s2[j - 1] ? score.mas : score.mis);
				is = s[j] + ((d[i - 1][j>>2] & DEL[k]) == INS[k] ? score.ges : score.gos);
				ds = cs + ((d[i][(j - 1)>>2] & DEL[(j - 1)&3]) == DEL[(j - 1)&3] ? score.ges : score.gos);
				s[j - 1] = cs;

				d[i][j>>2] |= ms > is ? (ms > ds ? (cs = ms, MMH[k]) : (cs = ds, DEL[k])) : \
					(is > ds ? (cs = is, INS[k]) : (cs = ds, DEL[k]));
				if (j == s2_l) s[j] = cs;
			}
			// printf("%d %d %d %d\n",i,j,cs, d[i][j>>2] & DEL[j&3]);
		}
	}
	// printf("cs:%d\n",cs );
	free (s);
}

void align_nd(const char *s2, const uint32_t s2_l, const char *s1, const uint32_t s1_l, alignment *aln){
	uint64_t s2_size = (s2_l >> 2) + 1;
	uint64_t size = s2_size * (s1_l + 1) * sizeof(uint8_t) + (s1_l + 1) * sizeof(uint8_t *);

	uint8_t **d = (uint8_t **) calloc(1, size);
	// assert (d);
	d[0] = (uint8_t *)(d + s1_l + 1);
	uint32_t i;
	for (i = 1; i < s1_l + 1; i++ ) d[i] = d[i - 1] + s2_size;

	// alignment aln;
	// aln.t_aln_str = malloc(s1_l + s2_l);
	// aln.q_aln_str = malloc(s1_l + s2_l);
	align_nd_dp(s1, s2, s1_l, s2_l, d);
	align_nd_tb(s1, s2, s1_l, s2_l, d, aln);
	// printf("%s\n%s\n", aln.t_aln_str, aln.q_aln_str);
	free (d);
	// free (aln.t_aln_str);
	// free (aln.q_aln_str);
}
