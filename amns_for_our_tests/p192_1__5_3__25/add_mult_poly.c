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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8320643660829855810UL) + ((((uint64_t)op[1] * 18186633952674430655UL) + ((uint64_t)op[2] * 12812547944107205958UL) + ((uint64_t)op[3] * 16178720128374071534UL) + ((uint64_t)op[4] * 8150242636114708270UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 8150242636114708270UL) + ((uint64_t)op[1] * 8320643660829855810UL) + ((((uint64_t)op[2] * 18186633952674430655UL) + ((uint64_t)op[3] * 12812547944107205958UL) + ((uint64_t)op[4] * 16178720128374071534UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 16178720128374071534UL) + ((uint64_t)op[1] * 8150242636114708270UL) + ((uint64_t)op[2] * 8320643660829855810UL) + ((((uint64_t)op[3] * 18186633952674430655UL) + ((uint64_t)op[4] * 12812547944107205958UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 12812547944107205958UL) + ((uint64_t)op[1] * 16178720128374071534UL) + ((uint64_t)op[2] * 8150242636114708270UL) + ((uint64_t)op[3] * 8320643660829855810UL) + ((uint64_t)op[4] * 17666413710604188733UL);
	tmp_q[4] = ((uint64_t)op[0] * 18186633952674430655UL) + ((uint64_t)op[1] * 12812547944107205958UL) + ((uint64_t)op[2] * 16178720128374071534UL) + ((uint64_t)op[3] * 8150242636114708270UL) + ((uint64_t)op[4] * 8320643660829855810UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 409872543046L) + ((((int128)tmp_q[1] * 78489367650L) + ((int128)tmp_q[2] * 35242734858L) + ((int128)tmp_q[3] * 74729831622L) + ((int128)tmp_q[4] * 111945439227L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 111945439227L) + ((int128)tmp_q[1] * 409872543046L) + ((((int128)tmp_q[2] * 78489367650L) + ((int128)tmp_q[3] * 35242734858L) + ((int128)tmp_q[4] * 74729831622L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 74729831622L) + ((int128)tmp_q[1] * 111945439227L) + ((int128)tmp_q[2] * 409872543046L) + ((((int128)tmp_q[3] * 78489367650L) + ((int128)tmp_q[4] * 35242734858L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 35242734858L) + ((int128)tmp_q[1] * 74729831622L) + ((int128)tmp_q[2] * 111945439227L) + ((int128)tmp_q[3] * 409872543046L) + ((int128)tmp_q[4] * 235468102950L);
	tmp_zero[4] = ((int128)tmp_q[0] * 78489367650L) + ((int128)tmp_q[1] * 35242734858L) + ((int128)tmp_q[2] * 74729831622L) + ((int128)tmp_q[3] * 111945439227L) + ((int128)tmp_q[4] * 409872543046L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

