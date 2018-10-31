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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14828926053474120907UL) + ((((uint64_t)op[1] * 3447899831015103639UL) + ((uint64_t)op[2] * 9449824516962330875UL) + ((uint64_t)op[3] * 3287485000762709375UL) + ((uint64_t)op[4] * 17422756373506763957UL) + ((uint64_t)op[5] * 5464618139875026028UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 5464618139875026028UL) + ((uint64_t)op[1] * 14828926053474120907UL) + ((((uint64_t)op[2] * 3447899831015103639UL) + ((uint64_t)op[3] * 9449824516962330875UL) + ((uint64_t)op[4] * 3287485000762709375UL) + ((uint64_t)op[5] * 17422756373506763957UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 17422756373506763957UL) + ((uint64_t)op[1] * 5464618139875026028UL) + ((uint64_t)op[2] * 14828926053474120907UL) + ((((uint64_t)op[3] * 3447899831015103639UL) + ((uint64_t)op[4] * 9449824516962330875UL) + ((uint64_t)op[5] * 3287485000762709375UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 3287485000762709375UL) + ((uint64_t)op[1] * 17422756373506763957UL) + ((uint64_t)op[2] * 5464618139875026028UL) + ((uint64_t)op[3] * 14828926053474120907UL) + ((((uint64_t)op[4] * 3447899831015103639UL) + ((uint64_t)op[5] * 9449824516962330875UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 9449824516962330875UL) + ((uint64_t)op[1] * 3287485000762709375UL) + ((uint64_t)op[2] * 17422756373506763957UL) + ((uint64_t)op[3] * 5464618139875026028UL) + ((uint64_t)op[4] * 14828926053474120907UL) + ((uint64_t)op[5] * 17239499155075518195UL);
	tmp_q[5] = ((uint64_t)op[0] * 3447899831015103639UL) + ((uint64_t)op[1] * 9449824516962330875UL) + ((uint64_t)op[2] * 3287485000762709375UL) + ((uint64_t)op[3] * 17422756373506763957UL) + ((uint64_t)op[4] * 5464618139875026028UL) + ((uint64_t)op[5] * 14828926053474120907UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 931626306757L) + ((((int128)tmp_q[1] * 397795114156L) + ((int128)tmp_q[2] * 2556658929397L) + ((int128)tmp_q[3] * 1257242578425L) + ((int128)tmp_q[4] * 1491687606705L) - ((int128)tmp_q[5] * 3007090031339L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 3007090031339L) - ((int128)tmp_q[1] * 931626306757L) + ((((int128)tmp_q[2] * 397795114156L) + ((int128)tmp_q[3] * 2556658929397L) + ((int128)tmp_q[4] * 1257242578425L) + ((int128)tmp_q[5] * 1491687606705L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 1491687606705L) - ((int128)tmp_q[1] * 3007090031339L) - ((int128)tmp_q[2] * 931626306757L) + ((((int128)tmp_q[3] * 397795114156L) + ((int128)tmp_q[4] * 2556658929397L) + ((int128)tmp_q[5] * 1257242578425L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1257242578425L) + ((int128)tmp_q[1] * 1491687606705L) - ((int128)tmp_q[2] * 3007090031339L) - ((int128)tmp_q[3] * 931626306757L) + ((((int128)tmp_q[4] * 397795114156L) + ((int128)tmp_q[5] * 2556658929397L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 2556658929397L) + ((int128)tmp_q[1] * 1257242578425L) + ((int128)tmp_q[2] * 1491687606705L) - ((int128)tmp_q[3] * 3007090031339L) - ((int128)tmp_q[4] * 931626306757L) + ((int128)tmp_q[5] * 1988975570780L);
	tmp_zero[5] = ((int128)tmp_q[0] * 397795114156L) + ((int128)tmp_q[1] * 2556658929397L) + ((int128)tmp_q[2] * 1257242578425L) + ((int128)tmp_q[3] * 1491687606705L) - ((int128)tmp_q[4] * 3007090031339L) - ((int128)tmp_q[5] * 931626306757L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

