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
	tmp_q[0] = ((uint64_t)op[0] * 11077485357039528955UL) + ((((uint64_t)op[1] * 12128753823721858337UL) + ((uint64_t)op[2] * 16648503626517690548UL) + ((uint64_t)op[3] * 4052095835143552829UL) + ((uint64_t)op[4] * 7807170007713711862UL) + ((uint64_t)op[5] * 18307062554990472418UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 18307062554990472418UL) + ((uint64_t)op[1] * 11077485357039528955UL) + ((((uint64_t)op[2] * 12128753823721858337UL) + ((uint64_t)op[3] * 16648503626517690548UL) + ((uint64_t)op[4] * 4052095835143552829UL) + ((uint64_t)op[5] * 7807170007713711862UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 7807170007713711862UL) + ((uint64_t)op[1] * 18307062554990472418UL) + ((uint64_t)op[2] * 11077485357039528955UL) + ((((uint64_t)op[3] * 12128753823721858337UL) + ((uint64_t)op[4] * 16648503626517690548UL) + ((uint64_t)op[5] * 4052095835143552829UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 4052095835143552829UL) + ((uint64_t)op[1] * 7807170007713711862UL) + ((uint64_t)op[2] * 18307062554990472418UL) + ((uint64_t)op[3] * 11077485357039528955UL) + ((((uint64_t)op[4] * 12128753823721858337UL) + ((uint64_t)op[5] * 16648503626517690548UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 16648503626517690548UL) + ((uint64_t)op[1] * 4052095835143552829UL) + ((uint64_t)op[2] * 7807170007713711862UL) + ((uint64_t)op[3] * 18307062554990472418UL) + ((uint64_t)op[4] * 11077485357039528955UL) + ((uint64_t)op[5] * 5303536897480636837UL);
	tmp_q[5] = ((uint64_t)op[0] * 12128753823721858337UL) + ((uint64_t)op[1] * 16648503626517690548UL) + ((uint64_t)op[2] * 4052095835143552829UL) + ((uint64_t)op[3] * 7807170007713711862UL) + ((uint64_t)op[4] * 18307062554990472418UL) + ((uint64_t)op[5] * 11077485357039528955UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 550809442L) + ((((int128)tmp_q[1] * 575012525L) - ((int128)tmp_q[2] * 3486350862L) - ((int128)tmp_q[3] * 2300456764L) - ((int128)tmp_q[4] * 1362126553L) - ((int128)tmp_q[5] * 2563860981L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 2563860981L) - ((int128)tmp_q[1] * 550809442L) + ((((int128)tmp_q[2] * 575012525L) - ((int128)tmp_q[3] * 3486350862L) - ((int128)tmp_q[4] * 2300456764L) - ((int128)tmp_q[5] * 1362126553L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 1362126553L) - ((int128)tmp_q[1] * 2563860981L) - ((int128)tmp_q[2] * 550809442L) + ((((int128)tmp_q[3] * 575012525L) - ((int128)tmp_q[4] * 3486350862L) - ((int128)tmp_q[5] * 2300456764L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2300456764L) - ((int128)tmp_q[1] * 1362126553L) - ((int128)tmp_q[2] * 2563860981L) - ((int128)tmp_q[3] * 550809442L) + ((((int128)tmp_q[4] * 575012525L) - ((int128)tmp_q[5] * 3486350862L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 3486350862L) - ((int128)tmp_q[1] * 2300456764L) - ((int128)tmp_q[2] * 1362126553L) - ((int128)tmp_q[3] * 2563860981L) - ((int128)tmp_q[4] * 550809442L) + ((int128)tmp_q[5] * 2875062625L);
	tmp_zero[5] = ((int128)tmp_q[0] * 575012525L) - ((int128)tmp_q[1] * 3486350862L) - ((int128)tmp_q[2] * 2300456764L) - ((int128)tmp_q[3] * 1362126553L) - ((int128)tmp_q[4] * 2563860981L) - ((int128)tmp_q[5] * 550809442L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

