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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7359953313762005121UL) + ((((uint64_t)op[1] * 4245400261967413267UL) + ((uint64_t)op[2] * 17188879505621137441UL) + ((uint64_t)op[3] * 56420038528176179UL) + ((uint64_t)op[4] * 374257521551947979UL) + ((uint64_t)op[5] * 14131647670644102884UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 14131647670644102884UL) + ((uint64_t)op[1] * 7359953313762005121UL) + ((((uint64_t)op[2] * 4245400261967413267UL) + ((uint64_t)op[3] * 17188879505621137441UL) + ((uint64_t)op[4] * 56420038528176179UL) + ((uint64_t)op[5] * 374257521551947979UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 374257521551947979UL) + ((uint64_t)op[1] * 14131647670644102884UL) + ((uint64_t)op[2] * 7359953313762005121UL) + ((((uint64_t)op[3] * 4245400261967413267UL) + ((uint64_t)op[4] * 17188879505621137441UL) + ((uint64_t)op[5] * 56420038528176179UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 56420038528176179UL) + ((uint64_t)op[1] * 374257521551947979UL) + ((uint64_t)op[2] * 14131647670644102884UL) + ((uint64_t)op[3] * 7359953313762005121UL) + ((((uint64_t)op[4] * 4245400261967413267UL) + ((uint64_t)op[5] * 17188879505621137441UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 17188879505621137441UL) + ((uint64_t)op[1] * 56420038528176179UL) + ((uint64_t)op[2] * 374257521551947979UL) + ((uint64_t)op[3] * 14131647670644102884UL) + ((uint64_t)op[4] * 7359953313762005121UL) + ((uint64_t)op[5] * 9955943549774725082UL);
	tmp_q[5] = ((uint64_t)op[0] * 4245400261967413267UL) + ((uint64_t)op[1] * 17188879505621137441UL) + ((uint64_t)op[2] * 56420038528176179UL) + ((uint64_t)op[3] * 374257521551947979UL) + ((uint64_t)op[4] * 14131647670644102884UL) + ((uint64_t)op[5] * 7359953313762005121UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1063686395L) - ((((int128)tmp_q[1] * 811415611L) - ((int128)tmp_q[2] * 1106054064L) - ((int128)tmp_q[3] * 2118695753L) - ((int128)tmp_q[4] * 636573269L) - ((int128)tmp_q[5] * 1635991690L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1635991690L) + ((int128)tmp_q[1] * 1063686395L) - ((((int128)tmp_q[2] * 811415611L) - ((int128)tmp_q[3] * 1106054064L) - ((int128)tmp_q[4] * 2118695753L) - ((int128)tmp_q[5] * 636573269L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 636573269L) - ((int128)tmp_q[1] * 1635991690L) + ((int128)tmp_q[2] * 1063686395L) - ((((int128)tmp_q[3] * 811415611L) - ((int128)tmp_q[4] * 1106054064L) - ((int128)tmp_q[5] * 2118695753L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2118695753L) - ((int128)tmp_q[1] * 636573269L) - ((int128)tmp_q[2] * 1635991690L) + ((int128)tmp_q[3] * 1063686395L) - ((((int128)tmp_q[4] * 811415611L) - ((int128)tmp_q[5] * 1106054064L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1106054064L) - ((int128)tmp_q[1] * 2118695753L) - ((int128)tmp_q[2] * 636573269L) - ((int128)tmp_q[3] * 1635991690L) + ((int128)tmp_q[4] * 1063686395L) - ((int128)tmp_q[5] * 1622831222L);
	tmp_zero[5] = ((int128)tmp_q[0] * 811415611L) - ((int128)tmp_q[1] * 1106054064L) - ((int128)tmp_q[2] * 2118695753L) - ((int128)tmp_q[3] * 636573269L) - ((int128)tmp_q[4] * 1635991690L) + ((int128)tmp_q[5] * 1063686395L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

