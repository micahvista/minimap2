#ifndef CMAPPY_H
#define CMAPPY_H

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
#include "ksw2.h"
KSEQ_DECLARE(gzFile)

typedef struct {
	const char *ctg;
	int32_t ctg_start, ctg_end;
	int32_t qry_start, qry_end;
	int32_t blen, mlen, NM, ctg_len;
	uint8_t mapq, is_primary;
	int8_t strand, trans_strand;
	int32_t seg_id;
	int32_t n_cigar32;
	uint32_t *cigar32;
} mm_hitpy_t;

static inline void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
{
	h->ctg = mi->seq[r->rid].name;
	h->ctg_len = mi->seq[r->rid].len;
	h->ctg_start = r->rs, h->ctg_end = r->re;
	h->qry_start = r->qs, h->qry_end = r->qe;
	h->strand = r->rev? -1 : 1;
	h->mapq = r->mapq;
	h->mlen = r->mlen;
	h->blen = r->blen;
	h->NM = r->blen - r->mlen + r->p->n_ambi;
	h->trans_strand = r->p->trans_strand == 1? 1 : r->p->trans_strand == 2? -1 : 0;
	h->is_primary = (r->id == r->parent);
	h->seg_id = r->seg_id;
	h->n_cigar32 = r->p->n_cigar;
	h->cigar32 = r->p->cigar;
}

static inline void mm_free_reg1(mm_reg1_t *r)
{
	free(r->p);
}

static inline kseq_t *mm_fastx_open(const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-") != 0? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	return kseq_init(fp);
}

static inline void mm_fastx_close(kseq_t *ks)
{
	gzFile fp;
	fp = ks->f->f;
	kseq_destroy(ks);
	gzclose(fp);
}

static inline int mm_verbose_level(int v)
{
	if (v >= 0) mm_verbose = v;
	return mm_verbose;
}

static inline void mm_reset_timer(void)
{
	extern double realtime(void);
	mm_realtime0 = realtime();
}

extern unsigned char seq_comp_table[256];
static inline mm_reg1_t *mm_map_aux(const mm_idx_t *mi, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt)
{
	mm_reg1_t *r;

	Py_BEGIN_ALLOW_THREADS
	if (seq2 == 0) {
		r = mm_map(mi, strlen(seq1), seq1, n_regs, b, opt, NULL);
	} else {
		int _n_regs[2];
		mm_reg1_t *regs[2];
		char *seq[2];
		int i, len[2];

		len[0] = strlen(seq1);
		len[1] = strlen(seq2);
		seq[0] = (char*)seq1;
		seq[1] = strdup(seq2);
		for (i = 0; i < len[1]>>1; ++i) {
			int t = seq[1][len[1] - i - 1];
			seq[1][len[1] - i - 1] = seq_comp_table[(uint8_t)seq[1][i]];
			seq[1][i] = seq_comp_table[t];
		}
		if (len[1]&1) seq[1][len[1]>>1] = seq_comp_table[(uint8_t)seq[1][len[1]>>1]];
		mm_map_frag(mi, 2, len, (const char**)seq, _n_regs, regs, b, opt, NULL);
		for (i = 0; i < _n_regs[1]; ++i)
			regs[1][i].rev = !regs[1][i].rev;
		*n_regs = _n_regs[0] + _n_regs[1];
		regs[0] = (mm_reg1_t*)realloc(regs[0], sizeof(mm_reg1_t) * (*n_regs));
		memcpy(&regs[0][_n_regs[0]], regs[1], _n_regs[1] * sizeof(mm_reg1_t));
		free(regs[1]);
		r = regs[0];
	}
	Py_END_ALLOW_THREADS

	return r;
}

static inline char *mappy_revcomp(int len, const uint8_t *seq)
{
	int i;
	char *rev;
	rev = (char*)malloc(len + 1);
	for (i = 0; i < len; ++i)
		rev[len - i - 1] = seq_comp_table[seq[i]];
	rev[len] = 0;
	return rev;
}

static char *mappy_fetch_seq(const mm_idx_t *mi, const char *name, int st, int en, int *len)
{
	int i, rid;
	char *s;
	*len = 0;
	rid = mm_idx_name2id(mi, name);
	if (rid < 0) return 0;
	if ((uint32_t)st >= mi->seq[rid].len || st >= en) return 0;
	if (en < 0 || (uint32_t)en > mi->seq[rid].len)
		en = mi->seq[rid].len;
	s = (char*)malloc(en - st + 1);
	*len = mm_idx_getseq(mi, rid, st, en, (uint8_t*)s);
	for (i = 0; i < *len; ++i)
		s[i] = "ACGTN"[(uint8_t)s[i]];
	s[*len] = 0;
	return s;
}

static mm_idx_t *mappy_idx_seq(int w, int k, int is_hpc, int bucket_bits, const char *seq, int len)
{
	const char *fake_name = "N/A";
	char *s;
	mm_idx_t *mi;
	s = (char*)calloc(len + 1, 1);
	memcpy(s, seq, len);
	mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
	free(s);
	return mi;
}
typedef struct { int32_t h, e; } eh_t;

void test_ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int zdrop, int flag, ksw_extz_t *ez)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int32_t i, j, k, max_j = 0, gapoe = gapo + gape, n_col, *off = 0, with_cigar = !(flag&KSW_EZ_SCORE_ONLY);
	uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

	ksw_reset_extz(ez);

	// allocate memory
	if (w < 0) w = tlen > qlen? tlen : qlen;
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	qp = (int8_t*)kmalloc(km, qlen * m);
	eh = (eh_t*)kcalloc(km, qlen + 1, 8);
	if (with_cigar) {
		z = (uint8_t*)kmalloc(km, (size_t)n_col * tlen);
		off = (int32_t*)kcalloc(km, tlen, 4);
	}

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = 0, eh[0].e = -gapoe - gapoe;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(gapoe + gape * (j - 1)), eh[j].e = -(gapoe + gapoe + gape * j);
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = KSW_NEG_INF; // everything is -inf outside the band

	// DP loop
	for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
		int32_t f, h1, st, en, max = KSW_NEG_INF;
		int8_t *q = &qp[target[i] * qlen];
		st = i > w? i - w : 0;
		en = i + w < qlen - 1? i + w : qlen - 1;
		h1 = st > 0? KSW_NEG_INF : -(gapoe + gape * i);
		f  = st > 0? KSW_NEG_INF : -(gapoe + gapoe + gape * i);
		if (!with_cigar) {
			for (j = st; j <= en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				p->h = h1;
				h += q[j];
				h = h >= e? h : e;
				h = h >= f? h : f;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				e  = e > h? e : h;
				p->e = e;
				f -= gape;
				f  = f > h? f : h;
			}
		} else if (!(flag&KSW_EZ_RIGHT)) {
			uint8_t *zi = &z[(long)i * n_col];
			off[i] = st;
			for (j = st; j <= en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				h += q[j];
				d = h >= e? 0 : 1;
				h = h >= e? h : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				max_j = max > h? max_j : j;
				max   = max > h? max   : h;
				h -= gapoe;
				e -= gape;
				d |= e > h? 0x08 : 0;
				e  = e > h? e    : h;
				p->e = e;
				f -= gape;
				d |= f > h? 0x10 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > h? f    : h;
				zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			uint8_t *zi = &z[(long)i * n_col];
			off[i] = st;
			for (j = st; j <= en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				h += q[j];
				d = h > e? 0 : 1;
				h = h > e? h : e;
				d = h > f? d : 2;
				h = h > f? h : f;
				h1 = h;
				max_j = max >= h? max_j : j;
				max   = max >= h? max   : h;
				h -= gapoe;
				e -= gape;
				d |= e >= h? 0x08 : 0;
				e  = e >= h? e    : h;
				p->e = e;
				f -= gape;
				d |= f >= h? 0x10 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f >= h? f    : h;
				zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		}
		eh[j].h = h1, eh[j].e = KSW_NEG_INF;
		// update ez
		if (en == qlen - 1 && eh[qlen].h > ez->mqe)
			ez->mqe = eh[qlen].h, ez->mqe_t = i;
		if (i == tlen - 1)
			ez->mte = max, ez->mte_q = max_j;
		if (ksw_apply_zdrop(ez, 0, max, i, max_j, zdrop, gape)) break;
		if (i == tlen - 1 && en == qlen - 1)
			ez->score = eh[qlen].h;
	}
	kfree(km, qp); kfree(km, eh);
	if (with_cigar) {
		int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
		if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY))
			ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		else if (ez->max_t >= 0 && ez->max_q >= 0)
			ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
		kfree(km, z); kfree(km, off);
	}
}

ksw_extz_t test_align(ksw_extz_t ez, const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape)
{
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	//ksw_extz_t ez;

	
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	test_ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
	free(ts); free(qs);
	return ez;
}

#endif
