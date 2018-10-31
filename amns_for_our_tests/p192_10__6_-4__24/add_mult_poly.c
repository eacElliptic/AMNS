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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10602366616913909561UL) + ((((uint64_t)op[1] * 10066194094458858820UL) + ((uint64_t)op[2] * 17337206635493878141UL) + ((uint64_t)op[3] * 6102118086133997917UL) + ((uint64_t)op[4] * 16184293018460075553UL) + ((uint64_t)op[5] * 9340004877159931554UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 9340004877159931554UL) + ((uint64_t)op[1] * 10602366616913909561UL) + ((((uint64_t)op[2] * 10066194094458858820UL) + ((uint64_t)op[3] * 17337206635493878141UL) + ((uint64_t)op[4] * 6102118086133997917UL) + ((uint64_t)op[5] * 16184293018460075553UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 16184293018460075553UL) + ((uint64_t)op[1] * 9340004877159931554UL) + ((uint64_t)op[2] * 10602366616913909561UL) + ((((uint64_t)op[3] * 10066194094458858820UL) + ((uint64_t)op[4] * 17337206635493878141UL) + ((uint64_t)op[5] * 6102118086133997917UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 6102118086133997917UL) + ((uint64_t)op[1] * 16184293018460075553UL) + ((uint64_t)op[2] * 9340004877159931554UL) + ((uint64_t)op[3] * 10602366616913909561UL) + ((((uint64_t)op[4] * 10066194094458858820UL) + ((uint64_t)op[5] * 17337206635493878141UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 17337206635493878141UL) + ((uint64_t)op[1] * 6102118086133997917UL) + ((uint64_t)op[2] * 16184293018460075553UL) + ((uint64_t)op[3] * 9340004877159931554UL) + ((uint64_t)op[4] * 10602366616913909561UL) + ((uint64_t)op[5] * 15075455843293219568UL);
	tmp_q[5] = ((uint64_t)op[0] * 10066194094458858820UL) + ((uint64_t)op[1] * 17337206635493878141UL) + ((uint64_t)op[2] * 6102118086133997917UL) + ((uint64_t)op[3] * 16184293018460075553UL) + ((uint64_t)op[4] * 9340004877159931554UL) + ((uint64_t)op[5] * 10602366616913909561UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1556973279L) - ((((int128)tmp_q[1] * 2512466984L) - ((int128)tmp_q[2] * 964776488L) + ((int128)tmp_q[3] * 1273278245L) + ((int128)tmp_q[4] * 677977005L) - ((int128)tmp_q[5] * 972926682L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 972926682L) + ((int128)tmp_q[1] * 1556973279L) - ((((int128)tmp_q[2] * 2512466984L) - ((int128)tmp_q[3] * 964776488L) + ((int128)tmp_q[4] * 1273278245L) + ((int128)tmp_q[5] * 677977005L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 677977005L) - ((int128)tmp_q[1] * 972926682L) + ((int128)tmp_q[2] * 1556973279L) - ((((int128)tmp_q[3] * 2512466984L) - ((int128)tmp_q[4] * 964776488L) + ((int128)tmp_q[5] * 1273278245L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 1273278245L) + ((int128)tmp_q[1] * 677977005L) - ((int128)tmp_q[2] * 972926682L) + ((int128)tmp_q[3] * 1556973279L) - ((((int128)tmp_q[4] * 2512466984L) - ((int128)tmp_q[5] * 964776488L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 964776488L) + ((int128)tmp_q[1] * 1273278245L) + ((int128)tmp_q[2] * 677977005L) - ((int128)tmp_q[3] * 972926682L) + ((int128)tmp_q[4] * 1556973279L) - ((int128)tmp_q[5] * 10049867936L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2512466984L) - ((int128)tmp_q[1] * 964776488L) + ((int128)tmp_q[2] * 1273278245L) + ((int128)tmp_q[3] * 677977005L) - ((int128)tmp_q[4] * 972926682L) + ((int128)tmp_q[5] * 1556973279L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

