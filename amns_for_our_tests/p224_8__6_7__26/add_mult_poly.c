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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3371410063712349152UL) + ((((uint64_t)op[1] * 661691391107292818UL) + ((uint64_t)op[2] * 2179876808387580023UL) + ((uint64_t)op[3] * 6751829466565173992UL) + ((uint64_t)op[4] * 6924131884169442910UL) + ((uint64_t)op[5] * 1279135205599464244UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 1279135205599464244UL) + ((uint64_t)op[1] * 3371410063712349152UL) + ((((uint64_t)op[2] * 661691391107292818UL) + ((uint64_t)op[3] * 2179876808387580023UL) + ((uint64_t)op[4] * 6751829466565173992UL) + ((uint64_t)op[5] * 6924131884169442910UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 6924131884169442910UL) + ((uint64_t)op[1] * 1279135205599464244UL) + ((uint64_t)op[2] * 3371410063712349152UL) + ((((uint64_t)op[3] * 661691391107292818UL) + ((uint64_t)op[4] * 2179876808387580023UL) + ((uint64_t)op[5] * 6751829466565173992UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 6751829466565173992UL) + ((uint64_t)op[1] * 6924131884169442910UL) + ((uint64_t)op[2] * 1279135205599464244UL) + ((uint64_t)op[3] * 3371410063712349152UL) + ((((uint64_t)op[4] * 661691391107292818UL) + ((uint64_t)op[5] * 2179876808387580023UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 2179876808387580023UL) + ((uint64_t)op[1] * 6751829466565173992UL) + ((uint64_t)op[2] * 6924131884169442910UL) + ((uint64_t)op[3] * 1279135205599464244UL) + ((uint64_t)op[4] * 3371410063712349152UL) + ((uint64_t)op[5] * 4631839737751049726UL);
	tmp_q[5] = ((uint64_t)op[0] * 661691391107292818UL) + ((uint64_t)op[1] * 2179876808387580023UL) + ((uint64_t)op[2] * 6751829466565173992UL) + ((uint64_t)op[3] * 6924131884169442910UL) + ((uint64_t)op[4] * 1279135205599464244UL) + ((uint64_t)op[5] * 3371410063712349152UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 7452886750L) + ((((int128)tmp_q[1] * 60285274532L) - ((int128)tmp_q[2] * 12800931680L) - ((int128)tmp_q[3] * 6207615986L) + ((int128)tmp_q[4] * 1597346687L) - ((int128)tmp_q[5] * 2616552800L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 2616552800L) - ((int128)tmp_q[1] * 7452886750L) + ((((int128)tmp_q[2] * 60285274532L) - ((int128)tmp_q[3] * 12800931680L) - ((int128)tmp_q[4] * 6207615986L) + ((int128)tmp_q[5] * 1597346687L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 1597346687L) - ((int128)tmp_q[1] * 2616552800L) - ((int128)tmp_q[2] * 7452886750L) + ((((int128)tmp_q[3] * 60285274532L) - ((int128)tmp_q[4] * 12800931680L) - ((int128)tmp_q[5] * 6207615986L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 6207615986L) + ((int128)tmp_q[1] * 1597346687L) - ((int128)tmp_q[2] * 2616552800L) - ((int128)tmp_q[3] * 7452886750L) + ((((int128)tmp_q[4] * 60285274532L) - ((int128)tmp_q[5] * 12800931680L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 12800931680L) - ((int128)tmp_q[1] * 6207615986L) + ((int128)tmp_q[2] * 1597346687L) - ((int128)tmp_q[3] * 2616552800L) - ((int128)tmp_q[4] * 7452886750L) + ((int128)tmp_q[5] * 421996921724L);
	tmp_zero[5] = ((int128)tmp_q[0] * 60285274532L) - ((int128)tmp_q[1] * 12800931680L) - ((int128)tmp_q[2] * 6207615986L) + ((int128)tmp_q[3] * 1597346687L) - ((int128)tmp_q[4] * 2616552800L) - ((int128)tmp_q[5] * 7452886750L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

