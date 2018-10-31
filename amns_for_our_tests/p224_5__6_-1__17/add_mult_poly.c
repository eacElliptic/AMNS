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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9324273273924526264UL) + ((((uint64_t)op[1] * 1731251275376707637UL) + ((uint64_t)op[2] * 9164510055787929996UL) + ((uint64_t)op[3] * 1498538187436803241UL) + ((uint64_t)op[4] * 3731388070744814967UL) + ((uint64_t)op[5] * 10835333356285826574UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 10835333356285826574UL) + ((uint64_t)op[1] * 9324273273924526264UL) + ((((uint64_t)op[2] * 1731251275376707637UL) + ((uint64_t)op[3] * 9164510055787929996UL) + ((uint64_t)op[4] * 1498538187436803241UL) + ((uint64_t)op[5] * 3731388070744814967UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 3731388070744814967UL) + ((uint64_t)op[1] * 10835333356285826574UL) + ((uint64_t)op[2] * 9324273273924526264UL) + ((((uint64_t)op[3] * 1731251275376707637UL) + ((uint64_t)op[4] * 9164510055787929996UL) + ((uint64_t)op[5] * 1498538187436803241UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 1498538187436803241UL) + ((uint64_t)op[1] * 3731388070744814967UL) + ((uint64_t)op[2] * 10835333356285826574UL) + ((uint64_t)op[3] * 9324273273924526264UL) + ((((uint64_t)op[4] * 1731251275376707637UL) + ((uint64_t)op[5] * 9164510055787929996UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 9164510055787929996UL) + ((uint64_t)op[1] * 1498538187436803241UL) + ((uint64_t)op[2] * 3731388070744814967UL) + ((uint64_t)op[3] * 10835333356285826574UL) + ((uint64_t)op[4] * 9324273273924526264UL) + ((uint64_t)op[5] * 16715492798332843979UL);
	tmp_q[5] = ((uint64_t)op[0] * 1731251275376707637UL) + ((uint64_t)op[1] * 9164510055787929996UL) + ((uint64_t)op[2] * 1498538187436803241UL) + ((uint64_t)op[3] * 3731388070744814967UL) + ((uint64_t)op[4] * 10835333356285826574UL) + ((uint64_t)op[5] * 9324273273924526264UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3076395023784178L) - (((int128)tmp_q[1] * 22103927292677259L) + ((int128)tmp_q[2] * 24426670602789218L) - ((int128)tmp_q[3] * 6938303930640773L) + ((int128)tmp_q[4] * 27503065626573397L) - ((int128)tmp_q[5] * 29042231223318030L));
	tmp_zero[1] = -((int128)tmp_q[0] * 29042231223318030L) + ((int128)tmp_q[1] * 3076395023784178L) - (((int128)tmp_q[2] * 22103927292677259L) + ((int128)tmp_q[3] * 24426670602789218L) - ((int128)tmp_q[4] * 6938303930640773L) + ((int128)tmp_q[5] * 27503065626573397L));
	tmp_zero[2] = ((int128)tmp_q[0] * 27503065626573397L) - ((int128)tmp_q[1] * 29042231223318030L) + ((int128)tmp_q[2] * 3076395023784178L) - (((int128)tmp_q[3] * 22103927292677259L) + ((int128)tmp_q[4] * 24426670602789218L) - ((int128)tmp_q[5] * 6938303930640773L));
	tmp_zero[3] = -((int128)tmp_q[0] * 6938303930640773L) + ((int128)tmp_q[1] * 27503065626573397L) - ((int128)tmp_q[2] * 29042231223318030L) + ((int128)tmp_q[3] * 3076395023784178L) - (((int128)tmp_q[4] * 22103927292677259L) + ((int128)tmp_q[5] * 24426670602789218L));
	tmp_zero[4] = ((int128)tmp_q[0] * 24426670602789218L) - ((int128)tmp_q[1] * 6938303930640773L) + ((int128)tmp_q[2] * 27503065626573397L) - ((int128)tmp_q[3] * 29042231223318030L) + ((int128)tmp_q[4] * 3076395023784178L) - ((int128)tmp_q[5] * 22103927292677259L);
	tmp_zero[5] = ((int128)tmp_q[0] * 22103927292677259L) + ((int128)tmp_q[1] * 24426670602789218L) - ((int128)tmp_q[2] * 6938303930640773L) + ((int128)tmp_q[3] * 27503065626573397L) - ((int128)tmp_q[4] * 29042231223318030L) + ((int128)tmp_q[5] * 3076395023784178L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

