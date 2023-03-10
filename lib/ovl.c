#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include "ovl.h"

bytes_t *init_encode_table(uint32_t n)
{
	int i, j, k, m;
	uint8_t tmp;
	bytes_t *ret = (bytes_t *)malloc(sizeof(bytes_t) * n);
	for(i = 0; i < n; i++) {
		k = i, j = 28, m = 0;
		if (k <= 127) ret[i].bytes[m++] = (k & 127);
		else {
			while (j >= 0) {
				tmp = ((k >> j) & 127);
				if (tmp > 0 || m > 0) ret[i].bytes[m++] = (tmp | 128);
				j -= 7;
			}
			ret[i].bytes[m - 1] &= 127;
		}
		ret[i].n = m;
	}
	return ret;
}

uint32_t *init_decode_table(uint64_t n)
{
	uint32_t *ret = (uint32_t *)malloc(sizeof(uint32_t) * n);
	uint32_t raw_n;
	uint64_t i, j;
	uint8_t m;
	for(i = 0; i < n; i++) {
		j = i, raw_n = 0, m = 0;
		while (j) {
			raw_n |= ((j & 127) << m);
			m += 7;
			j >>= 8;
		}
		ret[i] = raw_n;
	}
	return ret;
}

void destroy_table(void *arr)
{
	free(arr);
}

buffer_t *init_buffer(uint32_t buffersize)
{
	buffer_t *ret = (buffer_t *)calloc(1, sizeof(buffer_t));
	ret->buffer = (uint8_t *)malloc(sizeof(uint8_t) * buffersize);
	assert (ret->buffer);
	return ret;
}

void flush_buffer(FILE *fp, buffer_t *buf)
{
	if (fp && buf->buffer_i)
		fwrite(buf->buffer, sizeof(uint8_t), buf->buffer_i, fp);
	free(buf->buffer);
	free(buf);
}

void init_ovl_mode(FILE *fp, int mode){
	if (mode == 10){
		uint8_t buf[2] = {0, 255};
		fwrite(buf, sizeof(uint8_t), 2, fp);
	}
}

int find_ovlt_mode(FILE *fp){
	int mode = 10;
	overlap_i ovl;
	memset(&ovl, 0, sizeof(overlap_i));
	char *line = NULL;
	size_t len = 0;
	if (getline(&line, &len, fp)!=-1){
		sscanf(line, "%u\t%hhu\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", &ovl.qname, &ovl.rev, &ovl.qs, \
			&ovl.qe, &ovl.tname, &ovl.ts, &ovl.te, &ovl.qlen, &ovl.tlen, &ovl.identity);
		if (ovl.tlen == ovl.identity && ovl.tlen == 0){
			mode = 8;
		}
		rewind (fp);
	}
	if (line) free(line);
	return mode;
}

int find_ovlb_mode(FILE *fp){
	int mode = 8;
	uint8_t buf[2];
	int bread = fread(buf, sizeof(uint8_t), 2, fp);
	if (bread == 2){
		if (buf[0] == 0 && buf[1] == 255 ) mode = 10;
		else rewind (fp);
	}else{
		fprintf(stderr, "[ERROR] failed to find ovl file type or empty file\n");
		exit(1);
	}
	return mode;
}

void encode_ovl(FILE *fp, bytes_t *arr, prev_t *id, overlap *ovl, buffer_t *buf)
{
	int i, j, m, tmp;
	uint32_t result[8];
	uint32_t tlen;
	result[3] = ovl->qe - ovl->qs, tlen = ovl->te - ovl->ts;

	if (ovl->qname >= id->prev_qname) result[0] = ovl->qname - id->prev_qname;
	else ovl->rev |= 0b10, result[0] = id->prev_qname - ovl->qname;
	id->prev_qname = ovl->qname;

	if (ovl->tname >= id->prev_tname)  result[4] = ovl->tname - id->prev_tname;
	else ovl->rev |= 0b100, result[4] = id->prev_tname - ovl->tname;
	id->prev_tname = ovl->tname;

	if (result[3] >= tlen) result[6] = result[3] - tlen;
	else ovl->rev |= 0b1000, result[6] = tlen - result[3];
	
	result[1] = ovl->rev; result[2] = ovl->qs; result[5] = ovl->ts; result[7] = ovl->match;

	for(i = 0; i < 8; i++) {
		if (result[i] < NUM_LIMIT) {
			memcpy(buf->buffer + buf->buffer_i, arr[result[i]].bytes, arr[result[i]].n);
			buf->buffer_i += arr[result[i]].n;
		} else {
			j = 28, m = 0;
			uint8_t buffer1[5];
			while (j >= 0) {
				tmp = ((result[i] >> j) & 127);
				if (tmp > 0 || m > 0) buffer1[m++] = (tmp | 128);
				j -= 7;
			}
			buffer1[m - 1] &= 127;
			memcpy(buf->buffer + buf->buffer_i, buffer1, m);
			buf->buffer_i += m;
		}
	}
	if (buf->buffer_i >= 100000000) {
		fwrite(buf->buffer, sizeof(uint8_t), buf->buffer_i, fp);
		buf->buffer_i = 0;
	}
}

int decode_ovl(FILE *fp, uint32_t *arr, prev_t *id, uint32_t *ovl, buffer_t *buf, int n)
{
	int i = 0, m;
	uint32_t raw_n;
	uint64_t bits_cat = 0;

	if (!buf->read_n) {
		buf->read_n = fread(buf->buffer, sizeof(uint8_t), DCBUFSIZE, fp);
		if (!buf->read_n) {
			buf->buffer_i = 0;
			return -1;
		}
	}

	while(1) {
		bits_cat = (bits_cat << 8) | buf->buffer[buf->buffer_i];

		if (buf->buffer[buf->buffer_i] < 128) {
			if (bits_cat < NUM_LIMIT)
				ovl[i++] = arr[bits_cat];
			else {
				m = 0, raw_n = 0;
				while (bits_cat) {
					raw_n |= ((bits_cat & 127) << m);
					m += 7;
					bits_cat >>= 8;
				}
				ovl[i++] = raw_n;
			}
			bits_cat = 0;
		}

		if (++(buf->buffer_i) == buf->read_n) {
			buf->read_n = fread(buf->buffer, sizeof(uint8_t), DCBUFSIZE, fp);
			buf->buffer_i = 0;
		}
		// printf("%u %u\n", ovl[0], ovl[1]);
		if (i == n) {
			if (ovl[1] & 0b10) ovl[0] = id->prev_qname - ovl[0];
			else ovl[0] = id->prev_qname + ovl[0];
			id->prev_qname = ovl[0];
			if (ovl[1] & 0b100) ovl[4] = id->prev_tname - ovl[4];
			else ovl[4] = id->prev_tname + ovl[4];
			id->prev_tname = ovl[4];
			if (ovl[1] & 0b1000) ovl[6] = ovl[5] + ovl[3] + ovl[6];
			else ovl[6] = ovl[5] + ovl[3] - ovl[6];
			ovl[1] &= 1;
			ovl[3] += ovl[2];
			return 1;
		}
	}
}

void encode_ovl_i(FILE *fp, bytes_t *arr, prev_t *id, overlap_i *ovl, buffer_t *buf)
{
	int i, j, m, tmp;
	uint32_t tlen, result[10];
	result[3] = ovl->qe - ovl->qs, tlen = ovl->te - ovl->ts;

	if (ovl->qname >= id->prev_qname) result[0] = ovl->qname - id->prev_qname;
	else ovl->rev |= 0b10, result[0] = id->prev_qname - ovl->qname;

	if (ovl->tname >= id->prev_tname)  result[4] = ovl->tname - id->prev_tname;
	else ovl->rev |= 0b100, result[4] = id->prev_tname - ovl->tname;

	if (ovl->qname == id->prev_qname) result[7] = 0;
	else result[7] = ovl->qlen;

	if (ovl->tname == id->prev_tname) result[8] = 0;
	else result[8] = ovl->tlen;
	
	if (result[3] >= tlen) result[6] = result[3] - tlen;
	else ovl->rev |= 0b1000, result[6] = tlen - result[3];
	
	id->prev_qname = ovl->qname;
	id->prev_tname = ovl->tname;

	result[1] = ovl->rev; result[2] = ovl->qs; 
	result[5] = ovl->ts; result[9] = ovl->identity;

	for(i = 0; i < 10; i++) {
		if (result[i] < NUM_LIMIT) {
			memcpy(buf->buffer + buf->buffer_i, arr[result[i]].bytes, arr[result[i]].n);
			buf->buffer_i += arr[result[i]].n;
		} else {
			j = 28, m = 0;
			uint8_t buffer1[5];
			while (j >= 0) {
				tmp = ((result[i] >> j) & 127);
				if (tmp > 0 || m > 0) buffer1[m++] = (tmp | 128);
				j -= 7;
			}
			buffer1[m - 1] &= 127;
			memcpy(buf->buffer + buf->buffer_i, buffer1, m);
			buf->buffer_i += m;
		}
	}
	if (buf->buffer_i >= 100000000) {
		fwrite(buf->buffer, sizeof(uint8_t), buf->buffer_i, fp);
		buf->buffer_i = 0;
	}
}

static void init_aln(ovlinfo_aln *s){
	s->alnm = INIT_ALNM;
	s->alns = (aln *) calloc(s->alnm, sizeof (aln));
}

static void reallocate_aln(ovlinfo_aln *s){
  s->alnm += INIT_ALNM;
  s->alns = (aln *) realloc(s->alns, s->alnm * sizeof(aln));
  memset(s->alns + s->alnm - INIT_ALNM, 0, INIT_ALNM * sizeof(aln));
}

static void merge_aln(ovlinfo_aln *s){
	uint16_t i, j;
	aln t;
	for (i = 1; i < s->alnm; i++){//insert sort
		t = s->alns[i];
		for (j = i; j > 0 && s->alns[j-1].s > t.s; j--) s->alns[j] = s->alns[j-1];
		s->alns[j] = t;
	}

	i = 0;
	while (i < s->alnm - 1){
		if (!s->alns[i].e) {
			i ++; 
			continue;
		}

		for (j = i + 1; j < s->alnm; j++){
		  	if (s->alns[j].e <= s->alns[i].e){
				s->alns[j].e = 0;
		 	 }else if (s->alns[j].s <= s->alns[i].e && s->alns[j].e >= s->alns[i].e){
				s->alns[i].e = s->alns[j].e;
				s->alns[j].e = 0;
		 	 }else{
				break;
		  	}
		}
		i = j;
	}
}

static uint16_t find_alni(ovlinfo_aln *s){//return s->alni
	uint16_t i;
	if (s->alni != s->alnm - 1){
		for (i = 0; i < s->alnm; i++){
	  		if (!s->alns[i].e) return i;
		}
	}
	merge_aln(s);
	for (i = 0; i < s->alnm; i++){
		if (!s->alns[i].e) return i;
	}
	reallocate_aln(s);
	return i;
}

uint32_t find_alnse(ovlinfo_aln *s, aln *aln){
	uint32_t i, j = 0;
	merge_aln(s);
	for (i = 0; i < s->alnm; i++){//TODO: update reversed loop
		if (s->alns[i].e){
			aln->s = s->alns[i].s - EDGEBACKLEN;
			aln->e = s->alns[i].e + EDGEBACKLEN;
			j ++;
		}
	}
	return j;
}

static void fill_aln(ovlinfo_aln *s, uint32_t as, uint32_t ae){
	if (s->con < MAX_CON){
		s->alni = find_alni(s);
		s->alns[s->alni].s = as + EDGEBACKLEN;//redefine aln edge
		s->alns[s->alni].e = ae - EDGEBACKLEN;
  	}
}

static void fill_alnl(ovlinfo_aln *s, uint32_t as, uint32_t ae){
	if (s->con < MAX_CON && ae - as > s->alnl.e - s->alnl.s){
		s->alnl.s = as;
		s->alnl.e = ae;
  	}
}

void out_bl(khash_t(ovlh_) *os, FILE *f){
	int i;
	khiter_t k;
	ovlinfo_aln *ovl;
	for (k = kh_begin(os); k != kh_end(os); ++k){
		if (!kh_exist(os, k)) continue;
		ovl = &kh_val(os, k);
		if (ovl->con < MAX_CON){
			fprintf(f, "%u\t%u\t%hu\t%hu\t%u\t%u\t%u\t%u\t%u\t%u\t%u", kh_key(os, k), ovl->con,\
		 		ovl->lc, ovl->rc, ovl->lim, ovl->rim, ovl->llm, ovl->rlm, ovl->len, ovl->alnl.s, ovl->alnl.e);
			merge_aln(ovl);
			for (i = 0; i < ovl->alnm; i++){
				if (ovl->alns[i].e){
					fprintf(f, "\t%u\t%u", ovl->alns[i].s - EDGEBACKLEN, ovl->alns[i].e + EDGEBACKLEN);
		  		}
		  	}
		  	fprintf(f, "\n");
			free(ovl->alns);
		}else{
			fprintf(f, "%u\t%u\n", kh_key(os, k), ovl->con);
		}
	}
}

void read_bl(khash_t(ovlh_) *os, char *file, list *l){
	FILE *fp = fopen(file, "r");
	if (!fp){
		fprintf(stderr, "Failed open file %s!\n", file);
		exit(1);
	}

	khiter_t k;
	int absent;
	ovlinfo_aln *oli = NULL;
	char *line = NULL;
	size_t len = 0;
	ssize_t nread;
	size_t ols = sizeof(ovlinfo_aln);

	char *token;
	int i;
	uint64_t s;
	uint32_t n, e, ovl[11] = {0};
	while ((nread = getline(&line, &len, fp)) != -1) {
		line[nread - 1 ] = '\0';
		token = strtok(line, "\t");
		i = 0;
		ovl[i++] = atoll(token);
		while (i < 11){
			token = strtok(NULL, "\t");
			if (!token) break;
			ovl[i++] = atoll(token);
		}

		n = ovl[0];
		if (n + 1 > l->im) reallocate_list(l, (uint64_t) n + 1024 - l->im);
		if (!l->data[n]) n = l->data[n] = ++l->i;//skip the first ele (0)
		else n = l->data[n];

		k = kh_get(ovlh_, os, n);
		if (k != kh_end(os)){
			oli = &kh_val(os, k);
			if (oli->con >= MAX_CON) {
				continue;
			}else if (ovl[1]){
				oli->con += ovl[1];
				if (oli->con >= MAX_CON) {
					free (oli->alns);
					continue;
				}
			}
		}else{
			k = kh_put(ovlh_, os, n, &absent);
			oli = &kh_val(os, k);
			kh_key(os, k) = n;
			memset(oli, 0, ols);
			if (ovl[1]) {
				oli->con = ovl[1];
				if (oli->con >= MAX_CON) continue;
			}
			init_aln(oli);
			oli->len = ovl[8];
		}

		if (oli->lc > UINT16_MAX - ovl[2]) oli->lc = UINT16_MAX;
		else oli->lc += ovl[2];
		if (oli->rc > UINT16_MAX - ovl[3]) oli->rc = UINT16_MAX;
		else oli->rc += ovl[3];
		
		if (oli->lim < ovl[4]) oli->lim = ovl[4];
		if (oli->rim < ovl[5]) oli->rim = ovl[5];
		if (oli->llm < ovl[6]) oli->llm = ovl[6];
		if (oli->rlm < ovl[7]) oli->rlm = ovl[7];
		fill_alnl(oli, ovl[9], ovl[10]);
		
		while (1){
			token = strtok(NULL, "\t");
			if (!token) break;
			s = atoll(token);
			token = strtok(NULL, "\t");
			e = atoll(token);
			// if (ovl[1]) assert(s != e);
			if (ovl[1] && s == e) { s = 10; e = oli->len - 10; }
			fill_aln(oli, s, e);
		}
	}
	if (line) free (line);
	fclose (fp);
}

int filter_ovl(overlap_i *ovl, khash_t(ovlh_) *os, int32_t maxhan1, int32_t maxhan2){
  uint32_t alnlen;
  khint_t k;
  int absent;
  ovlinfo_aln *loli = NULL, *roli = NULL;
  size_t ols = sizeof(ovlinfo_aln);

  k = kh_get(ovlh_, os, ovl->qname);
  if (k != kh_end(os)){
	loli = &kh_val(os, k);
	// if (!loli->con){
	if (loli->con < MAX_CON){
	  if (ovl->qs <= maxhan2 && loli->lc < UINT16_MAX) loli->lc++;
	  if (ovl->qlen - ovl->qe <= maxhan2 && loli->rc < UINT16_MAX) loli->rc++;
	}
  }else{
	k = kh_put(ovlh_, os, ovl->qname, &absent);
	loli = &kh_val(os, k);
	memset(loli, 0, ols);
	loli->len = ovl->qlen;
	kh_key(os, k) = ovl->qname;
	init_aln(loli);
	if (ovl->qs <= maxhan2) loli->lc++;
	if (ovl->qlen - ovl->qe <= maxhan2) loli->rc++;
  }

  k = kh_get(ovlh_, os, ovl->tname);
  if (k != kh_end(os)){
	roli = &kh_val(os, k);
	// if (!roli->con){
	if (roli->con < MAX_CON){
	  if (ovl->ts <= maxhan2 && roli->rc < UINT16_MAX) roli->lc++;
	  if (ovl->tlen - ovl->te <= maxhan2 && roli->rc < UINT16_MAX) roli->rc++;
	}
  }else{
	k = kh_put(ovlh_, os, ovl->tname, &absent);
	roli = &kh_val(os, k);
	memset(roli, 0, ols);
	roli->len = ovl->tlen;
	kh_key(os, k) = ovl->tname;
	init_aln(roli);
	if (ovl->ts <= maxhan2) roli->lc++;
	if (ovl->tlen - ovl->te <= maxhan2) roli->rc ++;
  }

  loli = &kh_val(os, kh_get(ovlh_, os, ovl->qname));//the lolic address will be changed after realloc when insertting ovl->tname
  
  fill_aln(loli, ovl->qs, ovl->qe);
  fill_aln(roli, ovl->ts, ovl->te);

  if (loli->con < MAX_CON && ovl->qs <= maxhan2 && ovl->qe + maxhan2 >= ovl->qlen) {
	if (++loli->con >= MAX_CON){
		free (loli->alns);
		loli->alnm = 0;
	}
	return 0;
  }
  if (roli->con < MAX_CON && ovl->ts <= maxhan2 && ovl->te + maxhan2 >= ovl->tlen) {
	if (++roli->con >= MAX_CON){
		free (roli->alns);
		roli->alnm = 0;
	}
	return 0;
  }
  if (loli->con >= MAX_CON || roli->con >= MAX_CON) return 0;

  if (ovl->rev){
	if (ovl->qs <= maxhan1 && ovl->ts <= maxhan1){
	  if (ovl->qs <= maxhan2 && ovl->ts <= maxhan2){
		alnlen = max(ovl->qe - ovl->qs, ovl->te - ovl->ts);
		if (alnlen > loli->llm) loli->llm = alnlen;
		if (alnlen > roli->llm) roli->llm = alnlen;
		if (ovl->identity > loli->lim) loli->lim = ovl->identity;
		if (ovl->identity > roli->lim) roli->lim = ovl->identity;
	  }
	  return 1;
	}else if (ovl->qlen - ovl->qe <= maxhan1 && ovl->tlen - ovl->te <= maxhan1){
	  if (ovl->qlen - ovl->qe <= maxhan2 && ovl->tlen - ovl->te <= maxhan2){
		alnlen = max(ovl->qe - ovl->qs, ovl->te - ovl->ts);
		if (alnlen > loli->rlm) loli->rlm = alnlen;
		if (alnlen > roli->rlm) roli->rlm = alnlen;
		if (ovl->identity > loli->rim) loli->rim = ovl->identity;
		if (ovl->identity > roli->rim) roli->rim = ovl->identity;
	  }
	  return 1;
	}
  }else{
	if (ovl->qlen - ovl->qe <= maxhan1 && ovl->ts <= maxhan1){
	  if (ovl->qlen - ovl->qe <= maxhan2 && ovl->ts <= maxhan2){
		alnlen = max(ovl->qe - ovl->qs, ovl->te - ovl->ts);
		if (alnlen > loli->rlm) loli->rlm = alnlen;
		if (alnlen > roli->llm) roli->llm = alnlen;
		if (ovl->identity > loli->rim) loli->rim = ovl->identity;
		if (ovl->identity > roli->lim) roli->lim = ovl->identity;
	  }
	  return 1;
	}else if (ovl->qs <= maxhan1 && ovl->tlen - ovl->te <= maxhan1){
	  if (ovl->qs <= maxhan2 && ovl->tlen - ovl->te <= maxhan2){
		alnlen = max(ovl->qe - ovl->qs, ovl->te - ovl->ts);
		if (alnlen > loli->llm) loli->llm = alnlen;
		if (alnlen > roli->rlm) roli->rlm = alnlen;
		if (ovl->identity > loli->lim) loli->lim = ovl->identity;
		if (ovl->identity > roli->rim) roli->rim = ovl->identity;
	  }
	  return 1;
	}
  }
  // should avoid del reads that be contained after edges clipped
  if (ovl->qs <= maxhan1 && ovl->qe + maxhan1 >= ovl->qlen) return 1;
  if (ovl->ts <= maxhan1 && ovl->te + maxhan1 >= ovl->tlen) return 1;

  fill_alnl(loli, ovl->qs, ovl->qe);// should avoid del reads that contain each other
  fill_alnl(roli, ovl->ts, ovl->te);// donot save contain & dovetail ovls
  return 0;
}

int init_list(list *l){
  l->i = 0;
  l->im = 65536;
  l->data = (uint32_t *) calloc(l->im, sizeof(uint32_t));
  return l->data ? 0 : -1;
}

int reallocate_list(list *l, uint64_t s){
  l->im += s;
  uint32_t *data = (uint32_t *) realloc(l->data, l->im * sizeof(uint32_t));
  if (!data) {
    l->im -= s;
    return -1;
  }
  l->data = data;
  memset(l->data + (l->im - s), 0, s * sizeof(uint32_t));
  return 0;
}

void destroy_list(list *l){
  free (l->data);
}
