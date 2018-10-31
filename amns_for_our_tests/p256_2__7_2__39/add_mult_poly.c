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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14442135405329045595UL) + ((((uint64_t)op[1] * 976706897008409451UL) + ((uint64_t)op[2] * 9255199444158421220UL) + ((uint64_t)op[3] * 10683906425504536543UL) + ((uint64_t)op[4] * 8330064425871501239UL) + ((uint64_t)op[5] * 9684666751195260401UL) + ((uint64_t)op[6] * 16764280232001760870UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 16764280232001760870UL) + ((uint64_t)op[1] * 14442135405329045595UL) + ((((uint64_t)op[2] * 976706897008409451UL) + ((uint64_t)op[3] * 9255199444158421220UL) + ((uint64_t)op[4] * 10683906425504536543UL) + ((uint64_t)op[5] * 8330064425871501239UL) + ((uint64_t)op[6] * 9684666751195260401UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 9684666751195260401UL) + ((uint64_t)op[1] * 16764280232001760870UL) + ((uint64_t)op[2] * 14442135405329045595UL) + ((((uint64_t)op[3] * 976706897008409451UL) + ((uint64_t)op[4] * 9255199444158421220UL) + ((uint64_t)op[5] * 10683906425504536543UL) + ((uint64_t)op[6] * 8330064425871501239UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 8330064425871501239UL) + ((uint64_t)op[1] * 9684666751195260401UL) + ((uint64_t)op[2] * 16764280232001760870UL) + ((uint64_t)op[3] * 14442135405329045595UL) + ((((uint64_t)op[4] * 976706897008409451UL) + ((uint64_t)op[5] * 9255199444158421220UL) + ((uint64_t)op[6] * 10683906425504536543UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 10683906425504536543UL) + ((uint64_t)op[1] * 8330064425871501239UL) + ((uint64_t)op[2] * 9684666751195260401UL) + ((uint64_t)op[3] * 16764280232001760870UL) + ((uint64_t)op[4] * 14442135405329045595UL) + ((((uint64_t)op[5] * 976706897008409451UL) + ((uint64_t)op[6] * 9255199444158421220UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 9255199444158421220UL) + ((uint64_t)op[1] * 10683906425504536543UL) + ((uint64_t)op[2] * 8330064425871501239UL) + ((uint64_t)op[3] * 9684666751195260401UL) + ((uint64_t)op[4] * 16764280232001760870UL) + ((uint64_t)op[5] * 14442135405329045595UL) + ((uint64_t)op[6] * 1953413794016818902UL);
	tmp_q[6] = ((uint64_t)op[0] * 976706897008409451UL) + ((uint64_t)op[1] * 9255199444158421220UL) + ((uint64_t)op[2] * 10683906425504536543UL) + ((uint64_t)op[3] * 8330064425871501239UL) + ((uint64_t)op[4] * 9684666751195260401UL) + ((uint64_t)op[5] * 16764280232001760870UL) + ((uint64_t)op[6] * 14442135405329045595UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 34712657623L) + ((-((int128)tmp_q[1] * 1422386513L) + ((int128)tmp_q[2] * 27676619428L) + ((int128)tmp_q[3] * 5421569620L) + ((int128)tmp_q[4] * 22436969755L) - ((int128)tmp_q[5] * 43559442093L) + ((int128)tmp_q[6] * 66119655690L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 66119655690L) + ((int128)tmp_q[1] * 34712657623L) + ((-((int128)tmp_q[2] * 1422386513L) + ((int128)tmp_q[3] * 27676619428L) + ((int128)tmp_q[4] * 5421569620L) + ((int128)tmp_q[5] * 22436969755L) - ((int128)tmp_q[6] * 43559442093L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 43559442093L) + ((int128)tmp_q[1] * 66119655690L) + ((int128)tmp_q[2] * 34712657623L) + ((-((int128)tmp_q[3] * 1422386513L) + ((int128)tmp_q[4] * 27676619428L) + ((int128)tmp_q[5] * 5421569620L) + ((int128)tmp_q[6] * 22436969755L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 22436969755L) - ((int128)tmp_q[1] * 43559442093L) + ((int128)tmp_q[2] * 66119655690L) + ((int128)tmp_q[3] * 34712657623L) + ((-((int128)tmp_q[4] * 1422386513L) + ((int128)tmp_q[5] * 27676619428L) + ((int128)tmp_q[6] * 5421569620L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 5421569620L) + ((int128)tmp_q[1] * 22436969755L) - ((int128)tmp_q[2] * 43559442093L) + ((int128)tmp_q[3] * 66119655690L) + ((int128)tmp_q[4] * 34712657623L) + ((-((int128)tmp_q[5] * 1422386513L) + ((int128)tmp_q[6] * 27676619428L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 27676619428L) + ((int128)tmp_q[1] * 5421569620L) + ((int128)tmp_q[2] * 22436969755L) - ((int128)tmp_q[3] * 43559442093L) + ((int128)tmp_q[4] * 66119655690L) + ((int128)tmp_q[5] * 34712657623L) - ((int128)tmp_q[6] * 2844773026L);
	tmp_zero[6] = -((int128)tmp_q[0] * 1422386513L) + ((int128)tmp_q[1] * 27676619428L) + ((int128)tmp_q[2] * 5421569620L) + ((int128)tmp_q[3] * 22436969755L) - ((int128)tmp_q[4] * 43559442093L) + ((int128)tmp_q[5] * 66119655690L) + ((int128)tmp_q[6] * 34712657623L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

