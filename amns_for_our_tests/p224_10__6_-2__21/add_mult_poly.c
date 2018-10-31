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
	tmp_q[0] = ((uint64_t)op[0] * 7204952325942082783UL) + ((((uint64_t)op[1] * 15555711153733723382UL) + ((uint64_t)op[2] * 9147462079789080739UL) + ((uint64_t)op[3] * 6968577415544813315UL) + ((uint64_t)op[4] * 4854303368299124121UL) + ((uint64_t)op[5] * 4748022948862617604UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 4748022948862617604UL) + ((uint64_t)op[1] * 7204952325942082783UL) + ((((uint64_t)op[2] * 15555711153733723382UL) + ((uint64_t)op[3] * 9147462079789080739UL) + ((uint64_t)op[4] * 6968577415544813315UL) + ((uint64_t)op[5] * 4854303368299124121UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4854303368299124121UL) + ((uint64_t)op[1] * 4748022948862617604UL) + ((uint64_t)op[2] * 7204952325942082783UL) + ((((uint64_t)op[3] * 15555711153733723382UL) + ((uint64_t)op[4] * 9147462079789080739UL) + ((uint64_t)op[5] * 6968577415544813315UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 6968577415544813315UL) + ((uint64_t)op[1] * 4854303368299124121UL) + ((uint64_t)op[2] * 4748022948862617604UL) + ((uint64_t)op[3] * 7204952325942082783UL) + ((((uint64_t)op[4] * 15555711153733723382UL) + ((uint64_t)op[5] * 9147462079789080739UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 9147462079789080739UL) + ((uint64_t)op[1] * 6968577415544813315UL) + ((uint64_t)op[2] * 4854303368299124121UL) + ((uint64_t)op[3] * 4748022948862617604UL) + ((uint64_t)op[4] * 7204952325942082783UL) + ((uint64_t)op[5] * 5782065839951656468UL);
	tmp_q[5] = ((uint64_t)op[0] * 15555711153733723382UL) + ((uint64_t)op[1] * 9147462079789080739UL) + ((uint64_t)op[2] * 6968577415544813315UL) + ((uint64_t)op[3] * 4854303368299124121UL) + ((uint64_t)op[4] * 4748022948862617604UL) + ((uint64_t)op[5] * 7204952325942082783UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 62535176911L) - ((-((int128)tmp_q[1] * 59151735352L) + ((int128)tmp_q[2] * 18395888242L) + ((int128)tmp_q[3] * 36716126989L) - ((int128)tmp_q[4] * 88517202935L) - ((int128)tmp_q[5] * 16748853906L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 16748853906L) - ((int128)tmp_q[1] * 62535176911L) - ((-((int128)tmp_q[2] * 59151735352L) + ((int128)tmp_q[3] * 18395888242L) + ((int128)tmp_q[4] * 36716126989L) - ((int128)tmp_q[5] * 88517202935L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 88517202935L) - ((int128)tmp_q[1] * 16748853906L) - ((int128)tmp_q[2] * 62535176911L) - ((-((int128)tmp_q[3] * 59151735352L) + ((int128)tmp_q[4] * 18395888242L) + ((int128)tmp_q[5] * 36716126989L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 36716126989L) - ((int128)tmp_q[1] * 88517202935L) - ((int128)tmp_q[2] * 16748853906L) - ((int128)tmp_q[3] * 62535176911L) - ((-((int128)tmp_q[4] * 59151735352L) + ((int128)tmp_q[5] * 18395888242L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 18395888242L) + ((int128)tmp_q[1] * 36716126989L) - ((int128)tmp_q[2] * 88517202935L) - ((int128)tmp_q[3] * 16748853906L) - ((int128)tmp_q[4] * 62535176911L) + ((int128)tmp_q[5] * 118303470704L);
	tmp_zero[5] = -((int128)tmp_q[0] * 59151735352L) + ((int128)tmp_q[1] * 18395888242L) + ((int128)tmp_q[2] * 36716126989L) - ((int128)tmp_q[3] * 88517202935L) - ((int128)tmp_q[4] * 16748853906L) - ((int128)tmp_q[5] * 62535176911L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

