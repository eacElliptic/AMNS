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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9487596030508211967UL) + ((((uint64_t)op[1] * 17249325730905263310UL) + ((uint64_t)op[2] * 11298882668436615892UL) + ((uint64_t)op[3] * 1911044402704598200UL) + ((uint64_t)op[4] * 12421411268429470353UL) + ((uint64_t)op[5] * 15816883361539579985UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 15816883361539579985UL) + ((uint64_t)op[1] * 9487596030508211967UL) + ((((uint64_t)op[2] * 17249325730905263310UL) + ((uint64_t)op[3] * 11298882668436615892UL) + ((uint64_t)op[4] * 1911044402704598200UL) + ((uint64_t)op[5] * 12421411268429470353UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 12421411268429470353UL) + ((uint64_t)op[1] * 15816883361539579985UL) + ((uint64_t)op[2] * 9487596030508211967UL) + ((((uint64_t)op[3] * 17249325730905263310UL) + ((uint64_t)op[4] * 11298882668436615892UL) + ((uint64_t)op[5] * 1911044402704598200UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 1911044402704598200UL) + ((uint64_t)op[1] * 12421411268429470353UL) + ((uint64_t)op[2] * 15816883361539579985UL) + ((uint64_t)op[3] * 9487596030508211967UL) + ((((uint64_t)op[4] * 17249325730905263310UL) + ((uint64_t)op[5] * 11298882668436615892UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 11298882668436615892UL) + ((uint64_t)op[1] * 1911044402704598200UL) + ((uint64_t)op[2] * 12421411268429470353UL) + ((uint64_t)op[3] * 15816883361539579985UL) + ((uint64_t)op[4] * 9487596030508211967UL) + ((uint64_t)op[5] * 4789673371217153224UL);
	tmp_q[5] = ((uint64_t)op[0] * 17249325730905263310UL) + ((uint64_t)op[1] * 11298882668436615892UL) + ((uint64_t)op[2] * 1911044402704598200UL) + ((uint64_t)op[3] * 12421411268429470353UL) + ((uint64_t)op[4] * 15816883361539579985UL) + ((uint64_t)op[5] * 9487596030508211967UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 158251249053L) - ((((int128)tmp_q[1] * 2616215368414L) + ((int128)tmp_q[2] * 2904913808725L) - ((int128)tmp_q[3] * 4738625528929L) + ((int128)tmp_q[4] * 4340369110234L) + ((int128)tmp_q[5] * 2440647045621L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 2440647045621L) + ((int128)tmp_q[1] * 158251249053L) - ((((int128)tmp_q[2] * 2616215368414L) + ((int128)tmp_q[3] * 2904913808725L) - ((int128)tmp_q[4] * 4738625528929L) + ((int128)tmp_q[5] * 4340369110234L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 4340369110234L) + ((int128)tmp_q[1] * 2440647045621L) + ((int128)tmp_q[2] * 158251249053L) - ((((int128)tmp_q[3] * 2616215368414L) + ((int128)tmp_q[4] * 2904913808725L) - ((int128)tmp_q[5] * 4738625528929L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 4738625528929L) + ((int128)tmp_q[1] * 4340369110234L) + ((int128)tmp_q[2] * 2440647045621L) + ((int128)tmp_q[3] * 158251249053L) - ((((int128)tmp_q[4] * 2616215368414L) + ((int128)tmp_q[5] * 2904913808725L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 2904913808725L) - ((int128)tmp_q[1] * 4738625528929L) + ((int128)tmp_q[2] * 4340369110234L) + ((int128)tmp_q[3] * 2440647045621L) + ((int128)tmp_q[4] * 158251249053L) - ((int128)tmp_q[5] * 10464861473656L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2616215368414L) + ((int128)tmp_q[1] * 2904913808725L) - ((int128)tmp_q[2] * 4738625528929L) + ((int128)tmp_q[3] * 4340369110234L) + ((int128)tmp_q[4] * 2440647045621L) + ((int128)tmp_q[5] * 158251249053L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

