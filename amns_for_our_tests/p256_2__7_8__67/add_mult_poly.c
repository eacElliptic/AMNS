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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14379721310095346885UL) + ((((uint64_t)op[1] * 10116834308986863750UL) + ((uint64_t)op[2] * 1064583963861362970UL) + ((uint64_t)op[3] * 17314980960237593514UL) + ((uint64_t)op[4] * 12646402793355509000UL) + ((uint64_t)op[5] * 8927120968073414276UL) + ((uint64_t)op[6] * 12464694409136629642UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 12464694409136629642UL) + ((uint64_t)op[1] * 14379721310095346885UL) + ((((uint64_t)op[2] * 10116834308986863750UL) + ((uint64_t)op[3] * 1064583963861362970UL) + ((uint64_t)op[4] * 17314980960237593514UL) + ((uint64_t)op[5] * 12646402793355509000UL) + ((uint64_t)op[6] * 8927120968073414276UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 8927120968073414276UL) + ((uint64_t)op[1] * 12464694409136629642UL) + ((uint64_t)op[2] * 14379721310095346885UL) + ((((uint64_t)op[3] * 10116834308986863750UL) + ((uint64_t)op[4] * 1064583963861362970UL) + ((uint64_t)op[5] * 17314980960237593514UL) + ((uint64_t)op[6] * 12646402793355509000UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 12646402793355509000UL) + ((uint64_t)op[1] * 8927120968073414276UL) + ((uint64_t)op[2] * 12464694409136629642UL) + ((uint64_t)op[3] * 14379721310095346885UL) + ((((uint64_t)op[4] * 10116834308986863750UL) + ((uint64_t)op[5] * 1064583963861362970UL) + ((uint64_t)op[6] * 17314980960237593514UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 17314980960237593514UL) + ((uint64_t)op[1] * 12646402793355509000UL) + ((uint64_t)op[2] * 8927120968073414276UL) + ((uint64_t)op[3] * 12464694409136629642UL) + ((uint64_t)op[4] * 14379721310095346885UL) + ((((uint64_t)op[5] * 10116834308986863750UL) + ((uint64_t)op[6] * 1064583963861362970UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 1064583963861362970UL) + ((uint64_t)op[1] * 17314980960237593514UL) + ((uint64_t)op[2] * 12646402793355509000UL) + ((uint64_t)op[3] * 8927120968073414276UL) + ((uint64_t)op[4] * 12464694409136629642UL) + ((uint64_t)op[5] * 14379721310095346885UL) + ((uint64_t)op[6] * 7147698177056703536UL);
	tmp_q[6] = ((uint64_t)op[0] * 10116834308986863750UL) + ((uint64_t)op[1] * 1064583963861362970UL) + ((uint64_t)op[2] * 17314980960237593514UL) + ((uint64_t)op[3] * 12646402793355509000UL) + ((uint64_t)op[4] * 8927120968073414276UL) + ((uint64_t)op[5] * 12464694409136629642UL) + ((uint64_t)op[6] * 14379721310095346885UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 18591691379L) + ((((int128)tmp_q[1] * 57151693942L) - ((int128)tmp_q[2] * 7663205886L) + ((int128)tmp_q[3] * 207390506L) - ((int128)tmp_q[4] * 6805028416L) + ((int128)tmp_q[5] * 59394440048L) - ((int128)tmp_q[6] * 39648503494L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 39648503494L) + ((int128)tmp_q[1] * 18591691379L) + ((((int128)tmp_q[2] * 57151693942L) - ((int128)tmp_q[3] * 7663205886L) + ((int128)tmp_q[4] * 207390506L) - ((int128)tmp_q[5] * 6805028416L) + ((int128)tmp_q[6] * 59394440048L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 59394440048L) - ((int128)tmp_q[1] * 39648503494L) + ((int128)tmp_q[2] * 18591691379L) + ((((int128)tmp_q[3] * 57151693942L) - ((int128)tmp_q[4] * 7663205886L) + ((int128)tmp_q[5] * 207390506L) - ((int128)tmp_q[6] * 6805028416L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 6805028416L) + ((int128)tmp_q[1] * 59394440048L) - ((int128)tmp_q[2] * 39648503494L) + ((int128)tmp_q[3] * 18591691379L) + ((((int128)tmp_q[4] * 57151693942L) - ((int128)tmp_q[5] * 7663205886L) + ((int128)tmp_q[6] * 207390506L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 207390506L) - ((int128)tmp_q[1] * 6805028416L) + ((int128)tmp_q[2] * 59394440048L) - ((int128)tmp_q[3] * 39648503494L) + ((int128)tmp_q[4] * 18591691379L) + ((((int128)tmp_q[5] * 57151693942L) - ((int128)tmp_q[6] * 7663205886L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 7663205886L) + ((int128)tmp_q[1] * 207390506L) - ((int128)tmp_q[2] * 6805028416L) + ((int128)tmp_q[3] * 59394440048L) - ((int128)tmp_q[4] * 39648503494L) + ((int128)tmp_q[5] * 18591691379L) + ((int128)tmp_q[6] * 457213551536L);
	tmp_zero[6] = ((int128)tmp_q[0] * 57151693942L) - ((int128)tmp_q[1] * 7663205886L) + ((int128)tmp_q[2] * 207390506L) - ((int128)tmp_q[3] * 6805028416L) + ((int128)tmp_q[4] * 59394440048L) - ((int128)tmp_q[5] * 39648503494L) + ((int128)tmp_q[6] * 18591691379L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

