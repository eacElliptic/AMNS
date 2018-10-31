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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16756762169862383997UL) + ((((uint64_t)op[1] * 8850356329034668241UL) + ((uint64_t)op[2] * 17700995956139354475UL) + ((uint64_t)op[3] * 1458588786198674528UL) + ((uint64_t)op[4] * 7243308964266069457UL) + ((uint64_t)op[5] * 14991621472898168578UL) + ((uint64_t)op[6] * 7449119548593488340UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 7449119548593488340UL) + ((uint64_t)op[1] * 16756762169862383997UL) + ((((uint64_t)op[2] * 8850356329034668241UL) + ((uint64_t)op[3] * 17700995956139354475UL) + ((uint64_t)op[4] * 1458588786198674528UL) + ((uint64_t)op[5] * 7243308964266069457UL) + ((uint64_t)op[6] * 14991621472898168578UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 14991621472898168578UL) + ((uint64_t)op[1] * 7449119548593488340UL) + ((uint64_t)op[2] * 16756762169862383997UL) + ((((uint64_t)op[3] * 8850356329034668241UL) + ((uint64_t)op[4] * 17700995956139354475UL) + ((uint64_t)op[5] * 1458588786198674528UL) + ((uint64_t)op[6] * 7243308964266069457UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 7243308964266069457UL) + ((uint64_t)op[1] * 14991621472898168578UL) + ((uint64_t)op[2] * 7449119548593488340UL) + ((uint64_t)op[3] * 16756762169862383997UL) + ((((uint64_t)op[4] * 8850356329034668241UL) + ((uint64_t)op[5] * 17700995956139354475UL) + ((uint64_t)op[6] * 1458588786198674528UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 1458588786198674528UL) + ((uint64_t)op[1] * 7243308964266069457UL) + ((uint64_t)op[2] * 14991621472898168578UL) + ((uint64_t)op[3] * 7449119548593488340UL) + ((uint64_t)op[4] * 16756762169862383997UL) + ((((uint64_t)op[5] * 8850356329034668241UL) + ((uint64_t)op[6] * 17700995956139354475UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 17700995956139354475UL) + ((uint64_t)op[1] * 1458588786198674528UL) + ((uint64_t)op[2] * 7243308964266069457UL) + ((uint64_t)op[3] * 14991621472898168578UL) + ((uint64_t)op[4] * 7449119548593488340UL) + ((uint64_t)op[5] * 16756762169862383997UL) + ((uint64_t)op[6] * 746031415640215134UL);
	tmp_q[6] = ((uint64_t)op[0] * 8850356329034668241UL) + ((uint64_t)op[1] * 17700995956139354475UL) + ((uint64_t)op[2] * 1458588786198674528UL) + ((uint64_t)op[3] * 7243308964266069457UL) + ((uint64_t)op[4] * 14991621472898168578UL) + ((uint64_t)op[5] * 7449119548593488340UL) + ((uint64_t)op[6] * 16756762169862383997UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 61744767021L) - ((-((int128)tmp_q[1] * 56590647542L) - ((int128)tmp_q[2] * 21192021079L) + ((int128)tmp_q[3] * 40569272314L) - ((int128)tmp_q[4] * 30984609277L) + ((int128)tmp_q[5] * 16908894352L) + ((int128)tmp_q[6] * 28045092736L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 28045092736L) - ((int128)tmp_q[1] * 61744767021L) - ((-((int128)tmp_q[2] * 56590647542L) - ((int128)tmp_q[3] * 21192021079L) + ((int128)tmp_q[4] * 40569272314L) - ((int128)tmp_q[5] * 30984609277L) + ((int128)tmp_q[6] * 16908894352L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 16908894352L) + ((int128)tmp_q[1] * 28045092736L) - ((int128)tmp_q[2] * 61744767021L) - ((-((int128)tmp_q[3] * 56590647542L) - ((int128)tmp_q[4] * 21192021079L) + ((int128)tmp_q[5] * 40569272314L) - ((int128)tmp_q[6] * 30984609277L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 30984609277L) + ((int128)tmp_q[1] * 16908894352L) + ((int128)tmp_q[2] * 28045092736L) - ((int128)tmp_q[3] * 61744767021L) - ((-((int128)tmp_q[4] * 56590647542L) - ((int128)tmp_q[5] * 21192021079L) + ((int128)tmp_q[6] * 40569272314L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 40569272314L) - ((int128)tmp_q[1] * 30984609277L) + ((int128)tmp_q[2] * 16908894352L) + ((int128)tmp_q[3] * 28045092736L) - ((int128)tmp_q[4] * 61744767021L) - ((-((int128)tmp_q[5] * 56590647542L) - ((int128)tmp_q[6] * 21192021079L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 21192021079L) + ((int128)tmp_q[1] * 40569272314L) - ((int128)tmp_q[2] * 30984609277L) + ((int128)tmp_q[3] * 16908894352L) + ((int128)tmp_q[4] * 28045092736L) - ((int128)tmp_q[5] * 61744767021L) + ((int128)tmp_q[6] * 113181295084L);
	tmp_zero[6] = -((int128)tmp_q[0] * 56590647542L) - ((int128)tmp_q[1] * 21192021079L) + ((int128)tmp_q[2] * 40569272314L) - ((int128)tmp_q[3] * 30984609277L) + ((int128)tmp_q[4] * 16908894352L) + ((int128)tmp_q[5] * 28045092736L) - ((int128)tmp_q[6] * 61744767021L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

