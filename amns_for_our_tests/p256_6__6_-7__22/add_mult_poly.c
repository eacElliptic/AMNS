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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11571623691295764595UL) + ((((uint64_t)op[1] * 4753541634270011146UL) + ((uint64_t)op[2] * 15234026176489576114UL) + ((uint64_t)op[3] * 5045221442573186928UL) + ((uint64_t)op[4] * 14546478585944903972UL) + ((uint64_t)op[5] * 7388307173379955564UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 7388307173379955564UL) + ((uint64_t)op[1] * 11571623691295764595UL) + ((((uint64_t)op[2] * 4753541634270011146UL) + ((uint64_t)op[3] * 15234026176489576114UL) + ((uint64_t)op[4] * 5045221442573186928UL) + ((uint64_t)op[5] * 14546478585944903972UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 14546478585944903972UL) + ((uint64_t)op[1] * 7388307173379955564UL) + ((uint64_t)op[2] * 11571623691295764595UL) + ((((uint64_t)op[3] * 4753541634270011146UL) + ((uint64_t)op[4] * 15234026176489576114UL) + ((uint64_t)op[5] * 5045221442573186928UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 5045221442573186928UL) + ((uint64_t)op[1] * 14546478585944903972UL) + ((uint64_t)op[2] * 7388307173379955564UL) + ((uint64_t)op[3] * 11571623691295764595UL) + ((((uint64_t)op[4] * 4753541634270011146UL) + ((uint64_t)op[5] * 15234026176489576114UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 15234026176489576114UL) + ((uint64_t)op[1] * 5045221442573186928UL) + ((uint64_t)op[2] * 14546478585944903972UL) + ((uint64_t)op[3] * 7388307173379955564UL) + ((uint64_t)op[4] * 11571623691295764595UL) + ((uint64_t)op[5] * 3618696707529025210UL);
	tmp_q[5] = ((uint64_t)op[0] * 4753541634270011146UL) + ((uint64_t)op[1] * 15234026176489576114UL) + ((uint64_t)op[2] * 5045221442573186928UL) + ((uint64_t)op[3] * 14546478585944903972UL) + ((uint64_t)op[4] * 7388307173379955564UL) + ((uint64_t)op[5] * 11571623691295764595UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2016319288035L) - ((((int128)tmp_q[1] * 4712632383386L) - ((int128)tmp_q[2] * 467619866746L) + ((int128)tmp_q[3] * 5952321416016L) - ((int128)tmp_q[4] * 2899337885936L) - ((int128)tmp_q[5] * 782124043628L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 782124043628L) - ((int128)tmp_q[1] * 2016319288035L) - ((((int128)tmp_q[2] * 4712632383386L) - ((int128)tmp_q[3] * 467619866746L) + ((int128)tmp_q[4] * 5952321416016L) - ((int128)tmp_q[5] * 2899337885936L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 2899337885936L) - ((int128)tmp_q[1] * 782124043628L) - ((int128)tmp_q[2] * 2016319288035L) - ((((int128)tmp_q[3] * 4712632383386L) - ((int128)tmp_q[4] * 467619866746L) + ((int128)tmp_q[5] * 5952321416016L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 5952321416016L) - ((int128)tmp_q[1] * 2899337885936L) - ((int128)tmp_q[2] * 782124043628L) - ((int128)tmp_q[3] * 2016319288035L) - ((((int128)tmp_q[4] * 4712632383386L) - ((int128)tmp_q[5] * 467619866746L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 467619866746L) + ((int128)tmp_q[1] * 5952321416016L) - ((int128)tmp_q[2] * 2899337885936L) - ((int128)tmp_q[3] * 782124043628L) - ((int128)tmp_q[4] * 2016319288035L) - ((int128)tmp_q[5] * 32988426683702L);
	tmp_zero[5] = ((int128)tmp_q[0] * 4712632383386L) - ((int128)tmp_q[1] * 467619866746L) + ((int128)tmp_q[2] * 5952321416016L) - ((int128)tmp_q[3] * 2899337885936L) - ((int128)tmp_q[4] * 782124043628L) - ((int128)tmp_q[5] * 2016319288035L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

