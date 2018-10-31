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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15782746288515005607UL) + ((((uint64_t)op[1] * 11197633266556689677UL) + ((uint64_t)op[2] * 12751960239000173549UL) + ((uint64_t)op[3] * 3178761799426139819UL) + ((uint64_t)op[4] * 16419063794149399756UL) + ((uint64_t)op[5] * 11565093070865070856UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 11565093070865070856UL) + ((uint64_t)op[1] * 15782746288515005607UL) + ((((uint64_t)op[2] * 11197633266556689677UL) + ((uint64_t)op[3] * 12751960239000173549UL) + ((uint64_t)op[4] * 3178761799426139819UL) + ((uint64_t)op[5] * 16419063794149399756UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 16419063794149399756UL) + ((uint64_t)op[1] * 11565093070865070856UL) + ((uint64_t)op[2] * 15782746288515005607UL) + ((((uint64_t)op[3] * 11197633266556689677UL) + ((uint64_t)op[4] * 12751960239000173549UL) + ((uint64_t)op[5] * 3178761799426139819UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 3178761799426139819UL) + ((uint64_t)op[1] * 16419063794149399756UL) + ((uint64_t)op[2] * 11565093070865070856UL) + ((uint64_t)op[3] * 15782746288515005607UL) + ((((uint64_t)op[4] * 11197633266556689677UL) + ((uint64_t)op[5] * 12751960239000173549UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 12751960239000173549UL) + ((uint64_t)op[1] * 3178761799426139819UL) + ((uint64_t)op[2] * 16419063794149399756UL) + ((uint64_t)op[3] * 11565093070865070856UL) + ((uint64_t)op[4] * 15782746288515005607UL) + ((uint64_t)op[5] * 6601176695498068402UL);
	tmp_q[5] = ((uint64_t)op[0] * 11197633266556689677UL) + ((uint64_t)op[1] * 12751960239000173549UL) + ((uint64_t)op[2] * 3178761799426139819UL) + ((uint64_t)op[3] * 16419063794149399756UL) + ((uint64_t)op[4] * 11565093070865070856UL) + ((uint64_t)op[5] * 15782746288515005607UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6396094037L) - ((((int128)tmp_q[1] * 4475236841353L) + ((int128)tmp_q[2] * 516346770477L) + ((int128)tmp_q[3] * 2735063663053L) - ((int128)tmp_q[4] * 666335413098L) - ((int128)tmp_q[5] * 839817252500L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 839817252500L) - ((int128)tmp_q[1] * 6396094037L) - ((((int128)tmp_q[2] * 4475236841353L) + ((int128)tmp_q[3] * 516346770477L) + ((int128)tmp_q[4] * 2735063663053L) - ((int128)tmp_q[5] * 666335413098L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 666335413098L) - ((int128)tmp_q[1] * 839817252500L) - ((int128)tmp_q[2] * 6396094037L) - ((((int128)tmp_q[3] * 4475236841353L) + ((int128)tmp_q[4] * 516346770477L) + ((int128)tmp_q[5] * 2735063663053L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 2735063663053L) - ((int128)tmp_q[1] * 666335413098L) - ((int128)tmp_q[2] * 839817252500L) - ((int128)tmp_q[3] * 6396094037L) - ((((int128)tmp_q[4] * 4475236841353L) + ((int128)tmp_q[5] * 516346770477L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 516346770477L) + ((int128)tmp_q[1] * 2735063663053L) - ((int128)tmp_q[2] * 666335413098L) - ((int128)tmp_q[3] * 839817252500L) - ((int128)tmp_q[4] * 6396094037L) - ((int128)tmp_q[5] * 26851421048118L);
	tmp_zero[5] = ((int128)tmp_q[0] * 4475236841353L) + ((int128)tmp_q[1] * 516346770477L) + ((int128)tmp_q[2] * 2735063663053L) - ((int128)tmp_q[3] * 666335413098L) - ((int128)tmp_q[4] * 839817252500L) - ((int128)tmp_q[5] * 6396094037L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

