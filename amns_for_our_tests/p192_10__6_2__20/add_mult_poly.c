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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4283557503389965009UL) + ((((uint64_t)op[1] * 17546129512724576963UL) + ((uint64_t)op[2] * 11437858242989251663UL) + ((uint64_t)op[3] * 4930640501488296063UL) + ((uint64_t)op[4] * 14144014756451612365UL) + ((uint64_t)op[5] * 13102825759472624626UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 13102825759472624626UL) + ((uint64_t)op[1] * 4283557503389965009UL) + ((((uint64_t)op[2] * 17546129512724576963UL) + ((uint64_t)op[3] * 11437858242989251663UL) + ((uint64_t)op[4] * 4930640501488296063UL) + ((uint64_t)op[5] * 14144014756451612365UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 14144014756451612365UL) + ((uint64_t)op[1] * 13102825759472624626UL) + ((uint64_t)op[2] * 4283557503389965009UL) + ((((uint64_t)op[3] * 17546129512724576963UL) + ((uint64_t)op[4] * 11437858242989251663UL) + ((uint64_t)op[5] * 4930640501488296063UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 4930640501488296063UL) + ((uint64_t)op[1] * 14144014756451612365UL) + ((uint64_t)op[2] * 13102825759472624626UL) + ((uint64_t)op[3] * 4283557503389965009UL) + ((((uint64_t)op[4] * 17546129512724576963UL) + ((uint64_t)op[5] * 11437858242989251663UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 11437858242989251663UL) + ((uint64_t)op[1] * 4930640501488296063UL) + ((uint64_t)op[2] * 14144014756451612365UL) + ((uint64_t)op[3] * 13102825759472624626UL) + ((uint64_t)op[4] * 4283557503389965009UL) + ((uint64_t)op[5] * 16645514951739602310UL);
	tmp_q[5] = ((uint64_t)op[0] * 17546129512724576963UL) + ((uint64_t)op[1] * 11437858242989251663UL) + ((uint64_t)op[2] * 4930640501488296063UL) + ((uint64_t)op[3] * 14144014756451612365UL) + ((uint64_t)op[4] * 13102825759472624626UL) + ((uint64_t)op[5] * 4283557503389965009UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 822868225L) + ((((int128)tmp_q[1] * 2540373889L) + ((int128)tmp_q[2] * 763402178L) + ((int128)tmp_q[3] * 96173871L) - ((int128)tmp_q[4] * 2371674687L) + ((int128)tmp_q[5] * 402424724L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 402424724L) - ((int128)tmp_q[1] * 822868225L) + ((((int128)tmp_q[2] * 2540373889L) + ((int128)tmp_q[3] * 763402178L) + ((int128)tmp_q[4] * 96173871L) - ((int128)tmp_q[5] * 2371674687L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 2371674687L) + ((int128)tmp_q[1] * 402424724L) - ((int128)tmp_q[2] * 822868225L) + ((((int128)tmp_q[3] * 2540373889L) + ((int128)tmp_q[4] * 763402178L) + ((int128)tmp_q[5] * 96173871L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 96173871L) - ((int128)tmp_q[1] * 2371674687L) + ((int128)tmp_q[2] * 402424724L) - ((int128)tmp_q[3] * 822868225L) + ((((int128)tmp_q[4] * 2540373889L) + ((int128)tmp_q[5] * 763402178L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 763402178L) + ((int128)tmp_q[1] * 96173871L) - ((int128)tmp_q[2] * 2371674687L) + ((int128)tmp_q[3] * 402424724L) - ((int128)tmp_q[4] * 822868225L) + ((int128)tmp_q[5] * 5080747778L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2540373889L) + ((int128)tmp_q[1] * 763402178L) + ((int128)tmp_q[2] * 96173871L) - ((int128)tmp_q[3] * 2371674687L) + ((int128)tmp_q[4] * 402424724L) - ((int128)tmp_q[5] * 822868225L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

