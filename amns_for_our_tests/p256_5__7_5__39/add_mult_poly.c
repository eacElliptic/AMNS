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
	tmp_q[0] = ((uint64_t)op[0] * 4884016891187191485UL) + ((((uint64_t)op[1] * 5021725240399061340UL) + ((uint64_t)op[2] * 3434194271859835896UL) + ((uint64_t)op[3] * 12822131132697692655UL) + ((uint64_t)op[4] * 3196902914108992143UL) + ((uint64_t)op[5] * 13405109528684442954UL) + ((uint64_t)op[6] * 15717221962585689794UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 15717221962585689794UL) + ((uint64_t)op[1] * 4884016891187191485UL) + ((((uint64_t)op[2] * 5021725240399061340UL) + ((uint64_t)op[3] * 3434194271859835896UL) + ((uint64_t)op[4] * 12822131132697692655UL) + ((uint64_t)op[5] * 3196902914108992143UL) + ((uint64_t)op[6] * 13405109528684442954UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13405109528684442954UL) + ((uint64_t)op[1] * 15717221962585689794UL) + ((uint64_t)op[2] * 4884016891187191485UL) + ((((uint64_t)op[3] * 5021725240399061340UL) + ((uint64_t)op[4] * 3434194271859835896UL) + ((uint64_t)op[5] * 12822131132697692655UL) + ((uint64_t)op[6] * 3196902914108992143UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 3196902914108992143UL) + ((uint64_t)op[1] * 13405109528684442954UL) + ((uint64_t)op[2] * 15717221962585689794UL) + ((uint64_t)op[3] * 4884016891187191485UL) + ((((uint64_t)op[4] * 5021725240399061340UL) + ((uint64_t)op[5] * 3434194271859835896UL) + ((uint64_t)op[6] * 12822131132697692655UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 12822131132697692655UL) + ((uint64_t)op[1] * 3196902914108992143UL) + ((uint64_t)op[2] * 13405109528684442954UL) + ((uint64_t)op[3] * 15717221962585689794UL) + ((uint64_t)op[4] * 4884016891187191485UL) + ((((uint64_t)op[5] * 5021725240399061340UL) + ((uint64_t)op[6] * 3434194271859835896UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 3434194271859835896UL) + ((uint64_t)op[1] * 12822131132697692655UL) + ((uint64_t)op[2] * 3196902914108992143UL) + ((uint64_t)op[3] * 13405109528684442954UL) + ((uint64_t)op[4] * 15717221962585689794UL) + ((uint64_t)op[5] * 4884016891187191485UL) + ((uint64_t)op[6] * 6661882128285755084UL);
	tmp_q[6] = ((uint64_t)op[0] * 5021725240399061340UL) + ((uint64_t)op[1] * 3434194271859835896UL) + ((uint64_t)op[2] * 12822131132697692655UL) + ((uint64_t)op[3] * 3196902914108992143UL) + ((uint64_t)op[4] * 13405109528684442954UL) + ((uint64_t)op[5] * 15717221962585689794UL) + ((uint64_t)op[6] * 4884016891187191485UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 41953868079L) + ((((int128)tmp_q[1] * 66248773414L) + ((int128)tmp_q[2] * 22800675547L) + ((int128)tmp_q[3] * 50204106399L) + ((int128)tmp_q[4] * 4753983671L) - ((int128)tmp_q[5] * 29092464825L) - ((int128)tmp_q[6] * 82266177844L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 82266177844L) + ((int128)tmp_q[1] * 41953868079L) + ((((int128)tmp_q[2] * 66248773414L) + ((int128)tmp_q[3] * 22800675547L) + ((int128)tmp_q[4] * 50204106399L) + ((int128)tmp_q[5] * 4753983671L) - ((int128)tmp_q[6] * 29092464825L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 29092464825L) - ((int128)tmp_q[1] * 82266177844L) + ((int128)tmp_q[2] * 41953868079L) + ((((int128)tmp_q[3] * 66248773414L) + ((int128)tmp_q[4] * 22800675547L) + ((int128)tmp_q[5] * 50204106399L) + ((int128)tmp_q[6] * 4753983671L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 4753983671L) - ((int128)tmp_q[1] * 29092464825L) - ((int128)tmp_q[2] * 82266177844L) + ((int128)tmp_q[3] * 41953868079L) + ((((int128)tmp_q[4] * 66248773414L) + ((int128)tmp_q[5] * 22800675547L) + ((int128)tmp_q[6] * 50204106399L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 50204106399L) + ((int128)tmp_q[1] * 4753983671L) - ((int128)tmp_q[2] * 29092464825L) - ((int128)tmp_q[3] * 82266177844L) + ((int128)tmp_q[4] * 41953868079L) + ((((int128)tmp_q[5] * 66248773414L) + ((int128)tmp_q[6] * 22800675547L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 22800675547L) + ((int128)tmp_q[1] * 50204106399L) + ((int128)tmp_q[2] * 4753983671L) - ((int128)tmp_q[3] * 29092464825L) - ((int128)tmp_q[4] * 82266177844L) + ((int128)tmp_q[5] * 41953868079L) + ((int128)tmp_q[6] * 331243867070L);
	tmp_zero[6] = ((int128)tmp_q[0] * 66248773414L) + ((int128)tmp_q[1] * 22800675547L) + ((int128)tmp_q[2] * 50204106399L) + ((int128)tmp_q[3] * 4753983671L) - ((int128)tmp_q[4] * 29092464825L) - ((int128)tmp_q[5] * 82266177844L) + ((int128)tmp_q[6] * 41953868079L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

