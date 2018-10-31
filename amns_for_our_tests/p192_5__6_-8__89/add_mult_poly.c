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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5987658044980762543UL) + ((((uint64_t)op[1] * 9959544152377391898UL) + ((uint64_t)op[2] * 5793129926552854188UL) + ((uint64_t)op[3] * 12369749884086775555UL) + ((uint64_t)op[4] * 13930509852550628339UL) + ((uint64_t)op[5] * 13079259703407799812UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 13079259703407799812UL) + ((uint64_t)op[1] * 5987658044980762543UL) + ((((uint64_t)op[2] * 9959544152377391898UL) + ((uint64_t)op[3] * 5793129926552854188UL) + ((uint64_t)op[4] * 12369749884086775555UL) + ((uint64_t)op[5] * 13930509852550628339UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 13930509852550628339UL) + ((uint64_t)op[1] * 13079259703407799812UL) + ((uint64_t)op[2] * 5987658044980762543UL) + ((((uint64_t)op[3] * 9959544152377391898UL) + ((uint64_t)op[4] * 5793129926552854188UL) + ((uint64_t)op[5] * 12369749884086775555UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 12369749884086775555UL) + ((uint64_t)op[1] * 13930509852550628339UL) + ((uint64_t)op[2] * 13079259703407799812UL) + ((uint64_t)op[3] * 5987658044980762543UL) + ((((uint64_t)op[4] * 9959544152377391898UL) + ((uint64_t)op[5] * 5793129926552854188UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 5793129926552854188UL) + ((uint64_t)op[1] * 12369749884086775555UL) + ((uint64_t)op[2] * 13930509852550628339UL) + ((uint64_t)op[3] * 13079259703407799812UL) + ((uint64_t)op[4] * 5987658044980762543UL) + ((uint64_t)op[5] * 12557367149528622896UL);
	tmp_q[5] = ((uint64_t)op[0] * 9959544152377391898UL) + ((uint64_t)op[1] * 5793129926552854188UL) + ((uint64_t)op[2] * 12369749884086775555UL) + ((uint64_t)op[3] * 13930509852550628339UL) + ((uint64_t)op[4] * 13079259703407799812UL) + ((uint64_t)op[5] * 5987658044980762543UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 40019269370897L) - ((((int128)tmp_q[1] * 63448610686592L) - ((int128)tmp_q[2] * 76941030861467L) + ((int128)tmp_q[3] * 10313531636787L) + ((int128)tmp_q[4] * 41568023311059L) - ((int128)tmp_q[5] * 10705386853252L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 10705386853252L) + ((int128)tmp_q[1] * 40019269370897L) - ((((int128)tmp_q[2] * 63448610686592L) - ((int128)tmp_q[3] * 76941030861467L) + ((int128)tmp_q[4] * 10313531636787L) + ((int128)tmp_q[5] * 41568023311059L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 41568023311059L) - ((int128)tmp_q[1] * 10705386853252L) + ((int128)tmp_q[2] * 40019269370897L) - ((((int128)tmp_q[3] * 63448610686592L) - ((int128)tmp_q[4] * 76941030861467L) + ((int128)tmp_q[5] * 10313531636787L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 10313531636787L) + ((int128)tmp_q[1] * 41568023311059L) - ((int128)tmp_q[2] * 10705386853252L) + ((int128)tmp_q[3] * 40019269370897L) - ((((int128)tmp_q[4] * 63448610686592L) - ((int128)tmp_q[5] * 76941030861467L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 76941030861467L) + ((int128)tmp_q[1] * 10313531636787L) + ((int128)tmp_q[2] * 41568023311059L) - ((int128)tmp_q[3] * 10705386853252L) + ((int128)tmp_q[4] * 40019269370897L) - ((int128)tmp_q[5] * 507588885492736L);
	tmp_zero[5] = ((int128)tmp_q[0] * 63448610686592L) - ((int128)tmp_q[1] * 76941030861467L) + ((int128)tmp_q[2] * 10313531636787L) + ((int128)tmp_q[3] * 41568023311059L) - ((int128)tmp_q[4] * 10705386853252L) + ((int128)tmp_q[5] * 40019269370897L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

