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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9622615381895470154UL) + ((((uint64_t)op[1] * 14453630968845673888UL) + ((uint64_t)op[2] * 3983384586424273138UL) + ((uint64_t)op[3] * 12260376471481048288UL) + ((uint64_t)op[4] * 7677732396438064095UL) + ((uint64_t)op[5] * 468581808193135702UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 468581808193135702UL) + ((uint64_t)op[1] * 9622615381895470154UL) + ((((uint64_t)op[2] * 14453630968845673888UL) + ((uint64_t)op[3] * 3983384586424273138UL) + ((uint64_t)op[4] * 12260376471481048288UL) + ((uint64_t)op[5] * 7677732396438064095UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 7677732396438064095UL) + ((uint64_t)op[1] * 468581808193135702UL) + ((uint64_t)op[2] * 9622615381895470154UL) + ((((uint64_t)op[3] * 14453630968845673888UL) + ((uint64_t)op[4] * 3983384586424273138UL) + ((uint64_t)op[5] * 12260376471481048288UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 12260376471481048288UL) + ((uint64_t)op[1] * 7677732396438064095UL) + ((uint64_t)op[2] * 468581808193135702UL) + ((uint64_t)op[3] * 9622615381895470154UL) + ((((uint64_t)op[4] * 14453630968845673888UL) + ((uint64_t)op[5] * 3983384586424273138UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 3983384586424273138UL) + ((uint64_t)op[1] * 12260376471481048288UL) + ((uint64_t)op[2] * 7677732396438064095UL) + ((uint64_t)op[3] * 468581808193135702UL) + ((uint64_t)op[4] * 9622615381895470154UL) + ((uint64_t)op[5] * 1518821450609837024UL);
	tmp_q[5] = ((uint64_t)op[0] * 14453630968845673888UL) + ((uint64_t)op[1] * 3983384586424273138UL) + ((uint64_t)op[2] * 12260376471481048288UL) + ((uint64_t)op[3] * 7677732396438064095UL) + ((uint64_t)op[4] * 468581808193135702UL) + ((uint64_t)op[5] * 9622615381895470154UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6164358242L) - ((-((int128)tmp_q[1] * 135866863504L) - ((int128)tmp_q[2] * 29809126749L) - ((int128)tmp_q[3] * 2234065198L) - ((int128)tmp_q[4] * 156072978618L) + ((int128)tmp_q[5] * 18292648888L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 18292648888L) - ((int128)tmp_q[1] * 6164358242L) - ((-((int128)tmp_q[2] * 135866863504L) - ((int128)tmp_q[3] * 29809126749L) - ((int128)tmp_q[4] * 2234065198L) - ((int128)tmp_q[5] * 156072978618L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 156072978618L) + ((int128)tmp_q[1] * 18292648888L) - ((int128)tmp_q[2] * 6164358242L) - ((-((int128)tmp_q[3] * 135866863504L) - ((int128)tmp_q[4] * 29809126749L) - ((int128)tmp_q[5] * 2234065198L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2234065198L) - ((int128)tmp_q[1] * 156072978618L) + ((int128)tmp_q[2] * 18292648888L) - ((int128)tmp_q[3] * 6164358242L) - ((-((int128)tmp_q[4] * 135866863504L) - ((int128)tmp_q[5] * 29809126749L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 29809126749L) - ((int128)tmp_q[1] * 2234065198L) - ((int128)tmp_q[2] * 156072978618L) + ((int128)tmp_q[3] * 18292648888L) - ((int128)tmp_q[4] * 6164358242L) + ((int128)tmp_q[5] * 679334317520L);
	tmp_zero[5] = -((int128)tmp_q[0] * 135866863504L) - ((int128)tmp_q[1] * 29809126749L) - ((int128)tmp_q[2] * 2234065198L) - ((int128)tmp_q[3] * 156072978618L) + ((int128)tmp_q[4] * 18292648888L) - ((int128)tmp_q[5] * 6164358242L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

