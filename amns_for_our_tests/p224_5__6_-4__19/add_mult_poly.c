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
	tmp_q[0] = ((uint64_t)op[0] * 4968151464882148265UL) + ((((uint64_t)op[1] * 7682024467895800247UL) + ((uint64_t)op[2] * 16264602047014690269UL) + ((uint64_t)op[3] * 2885305659408419870UL) + ((uint64_t)op[4] * 8603599507878770886UL) + ((uint64_t)op[5] * 13043734881951223879UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 13043734881951223879UL) + ((uint64_t)op[1] * 4968151464882148265UL) + ((((uint64_t)op[2] * 7682024467895800247UL) + ((uint64_t)op[3] * 16264602047014690269UL) + ((uint64_t)op[4] * 2885305659408419870UL) + ((uint64_t)op[5] * 8603599507878770886UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 8603599507878770886UL) + ((uint64_t)op[1] * 13043734881951223879UL) + ((uint64_t)op[2] * 4968151464882148265UL) + ((((uint64_t)op[3] * 7682024467895800247UL) + ((uint64_t)op[4] * 16264602047014690269UL) + ((uint64_t)op[5] * 2885305659408419870UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 2885305659408419870UL) + ((uint64_t)op[1] * 8603599507878770886UL) + ((uint64_t)op[2] * 13043734881951223879UL) + ((uint64_t)op[3] * 4968151464882148265UL) + ((((uint64_t)op[4] * 7682024467895800247UL) + ((uint64_t)op[5] * 16264602047014690269UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 16264602047014690269UL) + ((uint64_t)op[1] * 2885305659408419870UL) + ((uint64_t)op[2] * 8603599507878770886UL) + ((uint64_t)op[3] * 13043734881951223879UL) + ((uint64_t)op[4] * 4968151464882148265UL) + ((uint64_t)op[5] * 6165390275835902244UL);
	tmp_q[5] = ((uint64_t)op[0] * 7682024467895800247UL) + ((uint64_t)op[1] * 16264602047014690269UL) + ((uint64_t)op[2] * 2885305659408419870UL) + ((uint64_t)op[3] * 8603599507878770886UL) + ((uint64_t)op[4] * 13043734881951223879UL) + ((uint64_t)op[5] * 4968151464882148265UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 26033196671L) - ((((int128)tmp_q[1] * 85969012922L) - ((int128)tmp_q[2] * 67310267086L) + ((int128)tmp_q[3] * 4346088653L) + ((int128)tmp_q[4] * 46863112537L) - ((int128)tmp_q[5] * 35961203665L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 35961203665L) + ((int128)tmp_q[1] * 26033196671L) - ((((int128)tmp_q[2] * 85969012922L) - ((int128)tmp_q[3] * 67310267086L) + ((int128)tmp_q[4] * 4346088653L) + ((int128)tmp_q[5] * 46863112537L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 46863112537L) - ((int128)tmp_q[1] * 35961203665L) + ((int128)tmp_q[2] * 26033196671L) - ((((int128)tmp_q[3] * 85969012922L) - ((int128)tmp_q[4] * 67310267086L) + ((int128)tmp_q[5] * 4346088653L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 4346088653L) + ((int128)tmp_q[1] * 46863112537L) - ((int128)tmp_q[2] * 35961203665L) + ((int128)tmp_q[3] * 26033196671L) - ((((int128)tmp_q[4] * 85969012922L) - ((int128)tmp_q[5] * 67310267086L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 67310267086L) + ((int128)tmp_q[1] * 4346088653L) + ((int128)tmp_q[2] * 46863112537L) - ((int128)tmp_q[3] * 35961203665L) + ((int128)tmp_q[4] * 26033196671L) - ((int128)tmp_q[5] * 343876051688L);
	tmp_zero[5] = ((int128)tmp_q[0] * 85969012922L) - ((int128)tmp_q[1] * 67310267086L) + ((int128)tmp_q[2] * 4346088653L) + ((int128)tmp_q[3] * 46863112537L) - ((int128)tmp_q[4] * 35961203665L) + ((int128)tmp_q[5] * 26033196671L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

