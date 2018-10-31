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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17719395207415537163UL) + ((((uint64_t)op[1] * 16941647354830115063UL) + ((uint64_t)op[2] * 8189001411382086272UL) + ((uint64_t)op[3] * 5681639607865276847UL) + ((uint64_t)op[4] * 6217902790011530934UL) + ((uint64_t)op[5] * 7533449277630506959UL) + ((uint64_t)op[6] * 11315839380184739796UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 11315839380184739796UL) + ((uint64_t)op[1] * 17719395207415537163UL) + ((((uint64_t)op[2] * 16941647354830115063UL) + ((uint64_t)op[3] * 8189001411382086272UL) + ((uint64_t)op[4] * 5681639607865276847UL) + ((uint64_t)op[5] * 6217902790011530934UL) + ((uint64_t)op[6] * 7533449277630506959UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 7533449277630506959UL) + ((uint64_t)op[1] * 11315839380184739796UL) + ((uint64_t)op[2] * 17719395207415537163UL) + ((((uint64_t)op[3] * 16941647354830115063UL) + ((uint64_t)op[4] * 8189001411382086272UL) + ((uint64_t)op[5] * 5681639607865276847UL) + ((uint64_t)op[6] * 6217902790011530934UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 6217902790011530934UL) + ((uint64_t)op[1] * 7533449277630506959UL) + ((uint64_t)op[2] * 11315839380184739796UL) + ((uint64_t)op[3] * 17719395207415537163UL) + ((((uint64_t)op[4] * 16941647354830115063UL) + ((uint64_t)op[5] * 8189001411382086272UL) + ((uint64_t)op[6] * 5681639607865276847UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 5681639607865276847UL) + ((uint64_t)op[1] * 6217902790011530934UL) + ((uint64_t)op[2] * 7533449277630506959UL) + ((uint64_t)op[3] * 11315839380184739796UL) + ((uint64_t)op[4] * 17719395207415537163UL) + ((((uint64_t)op[5] * 16941647354830115063UL) + ((uint64_t)op[6] * 8189001411382086272UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 8189001411382086272UL) + ((uint64_t)op[1] * 5681639607865276847UL) + ((uint64_t)op[2] * 6217902790011530934UL) + ((uint64_t)op[3] * 7533449277630506959UL) + ((uint64_t)op[4] * 11315839380184739796UL) + ((uint64_t)op[5] * 17719395207415537163UL) + ((uint64_t)op[6] * 12040773751035492424UL);
	tmp_q[6] = ((uint64_t)op[0] * 16941647354830115063UL) + ((uint64_t)op[1] * 8189001411382086272UL) + ((uint64_t)op[2] * 5681639607865276847UL) + ((uint64_t)op[3] * 6217902790011530934UL) + ((uint64_t)op[4] * 7533449277630506959UL) + ((uint64_t)op[5] * 11315839380184739796UL) + ((uint64_t)op[6] * 17719395207415537163UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 66576895405L) - ((((int128)tmp_q[1] * 1549745900L) + ((int128)tmp_q[2] * 50606876184L) + ((int128)tmp_q[3] * 22549711148L) - ((int128)tmp_q[4] * 55725930938L) - ((int128)tmp_q[5] * 61568993513L) - ((int128)tmp_q[6] * 56346071156L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 56346071156L) + ((int128)tmp_q[1] * 66576895405L) - ((((int128)tmp_q[2] * 1549745900L) + ((int128)tmp_q[3] * 50606876184L) + ((int128)tmp_q[4] * 22549711148L) - ((int128)tmp_q[5] * 55725930938L) - ((int128)tmp_q[6] * 61568993513L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 61568993513L) - ((int128)tmp_q[1] * 56346071156L) + ((int128)tmp_q[2] * 66576895405L) - ((((int128)tmp_q[3] * 1549745900L) + ((int128)tmp_q[4] * 50606876184L) + ((int128)tmp_q[5] * 22549711148L) - ((int128)tmp_q[6] * 55725930938L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 55725930938L) - ((int128)tmp_q[1] * 61568993513L) - ((int128)tmp_q[2] * 56346071156L) + ((int128)tmp_q[3] * 66576895405L) - ((((int128)tmp_q[4] * 1549745900L) + ((int128)tmp_q[5] * 50606876184L) + ((int128)tmp_q[6] * 22549711148L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 22549711148L) - ((int128)tmp_q[1] * 55725930938L) - ((int128)tmp_q[2] * 61568993513L) - ((int128)tmp_q[3] * 56346071156L) + ((int128)tmp_q[4] * 66576895405L) - ((((int128)tmp_q[5] * 1549745900L) + ((int128)tmp_q[6] * 50606876184L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 50606876184L) + ((int128)tmp_q[1] * 22549711148L) - ((int128)tmp_q[2] * 55725930938L) - ((int128)tmp_q[3] * 61568993513L) - ((int128)tmp_q[4] * 56346071156L) + ((int128)tmp_q[5] * 66576895405L) - ((int128)tmp_q[6] * 12397967200L);
	tmp_zero[6] = ((int128)tmp_q[0] * 1549745900L) + ((int128)tmp_q[1] * 50606876184L) + ((int128)tmp_q[2] * 22549711148L) - ((int128)tmp_q[3] * 55725930938L) - ((int128)tmp_q[4] * 61568993513L) - ((int128)tmp_q[5] * 56346071156L) + ((int128)tmp_q[6] * 66576895405L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

