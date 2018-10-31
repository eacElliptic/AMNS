#include "add_mult_poly.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ Computes pa(X)*pb(X) mod(X^n - c)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7960169084095952184UL) + ((((uint64_t)op[1] * 100325685030622257UL) + ((uint64_t)op[2] * 1114802945301216390UL) + ((uint64_t)op[3] * 10130638288373760383UL) + ((uint64_t)op[4] * 1622033291604652726UL) + ((uint64_t)op[5] * 125919318967444085UL) + ((uint64_t)op[6] * 14836007502709785738UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 14836007502709785738UL) + ((uint64_t)op[1] * 7960169084095952184UL) + ((((uint64_t)op[2] * 100325685030622257UL) + ((uint64_t)op[3] * 1114802945301216390UL) + ((uint64_t)op[4] * 10130638288373760383UL) + ((uint64_t)op[5] * 1622033291604652726UL) + ((uint64_t)op[6] * 125919318967444085UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 125919318967444085UL) + ((uint64_t)op[1] * 14836007502709785738UL) + ((uint64_t)op[2] * 7960169084095952184UL) + ((((uint64_t)op[3] * 100325685030622257UL) + ((uint64_t)op[4] * 1114802945301216390UL) + ((uint64_t)op[5] * 10130638288373760383UL) + ((uint64_t)op[6] * 1622033291604652726UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 1622033291604652726UL) + ((uint64_t)op[1] * 125919318967444085UL) + ((uint64_t)op[2] * 14836007502709785738UL) + ((uint64_t)op[3] * 7960169084095952184UL) + ((((uint64_t)op[4] * 100325685030622257UL) + ((uint64_t)op[5] * 1114802945301216390UL) + ((uint64_t)op[6] * 10130638288373760383UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 10130638288373760383UL) + ((uint64_t)op[1] * 1622033291604652726UL) + ((uint64_t)op[2] * 125919318967444085UL) + ((uint64_t)op[3] * 14836007502709785738UL) + ((uint64_t)op[4] * 7960169084095952184UL) + ((((uint64_t)op[5] * 100325685030622257UL) + ((uint64_t)op[6] * 1114802945301216390UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 1114802945301216390UL) + ((uint64_t)op[1] * 10130638288373760383UL) + ((uint64_t)op[2] * 1622033291604652726UL) + ((uint64_t)op[3] * 125919318967444085UL) + ((uint64_t)op[4] * 14836007502709785738UL) + ((uint64_t)op[5] * 7960169084095952184UL) + ((uint64_t)op[6] * 501628425153111285UL);
	tmp_q[6] = ((uint64_t)op[0] * 100325685030622257UL) + ((uint64_t)op[1] * 1114802945301216390UL) + ((uint64_t)op[2] * 10130638288373760383UL) + ((uint64_t)op[3] * 1622033291604652726UL) + ((uint64_t)op[4] * 125919318967444085UL) + ((uint64_t)op[5] * 14836007502709785738UL) + ((uint64_t)op[6] * 7960169084095952184UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 8521771674L) + ((-((int128)tmp_q[1] * 25924599484L) + ((int128)tmp_q[2] * 45080879679L) + ((int128)tmp_q[3] * 8791195451L) + ((int128)tmp_q[4] * 22539957737L) + ((int128)tmp_q[5] * 4786737387L) - ((int128)tmp_q[6] * 19876612267L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 19876612267L) + ((int128)tmp_q[1] * 8521771674L) + ((-((int128)tmp_q[2] * 25924599484L) + ((int128)tmp_q[3] * 45080879679L) + ((int128)tmp_q[4] * 8791195451L) + ((int128)tmp_q[5] * 22539957737L) + ((int128)tmp_q[6] * 4786737387L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 4786737387L) - ((int128)tmp_q[1] * 19876612267L) + ((int128)tmp_q[2] * 8521771674L) + ((-((int128)tmp_q[3] * 25924599484L) + ((int128)tmp_q[4] * 45080879679L) + ((int128)tmp_q[5] * 8791195451L) + ((int128)tmp_q[6] * 22539957737L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 22539957737L) + ((int128)tmp_q[1] * 4786737387L) - ((int128)tmp_q[2] * 19876612267L) + ((int128)tmp_q[3] * 8521771674L) + ((-((int128)tmp_q[4] * 25924599484L) + ((int128)tmp_q[5] * 45080879679L) + ((int128)tmp_q[6] * 8791195451L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 8791195451L) + ((int128)tmp_q[1] * 22539957737L) + ((int128)tmp_q[2] * 4786737387L) - ((int128)tmp_q[3] * 19876612267L) + ((int128)tmp_q[4] * 8521771674L) + ((-((int128)tmp_q[5] * 25924599484L) + ((int128)tmp_q[6] * 45080879679L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 45080879679L) + ((int128)tmp_q[1] * 8791195451L) + ((int128)tmp_q[2] * 22539957737L) + ((int128)tmp_q[3] * 4786737387L) - ((int128)tmp_q[4] * 19876612267L) + ((int128)tmp_q[5] * 8521771674L) - ((int128)tmp_q[6] * 129622997420L);
	tmp_zero[6] = -((int128)tmp_q[0] * 25924599484L) + ((int128)tmp_q[1] * 45080879679L) + ((int128)tmp_q[2] * 8791195451L) + ((int128)tmp_q[3] * 22539957737L) + ((int128)tmp_q[4] * 4786737387L) - ((int128)tmp_q[5] * 19876612267L) + ((int128)tmp_q[6] * 8521771674L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

