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
	tmp_q[0] = ((uint64_t)op[0] * 16950787534563821033UL) + ((((uint64_t)op[1] * 4611714075407134313UL) + ((uint64_t)op[2] * 5307564123049521121UL) + ((uint64_t)op[3] * 9341691724871993531UL) + ((uint64_t)op[4] * 4367203912113899393UL) + ((uint64_t)op[5] * 10502031650334300791UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 10502031650334300791UL) + ((uint64_t)op[1] * 16950787534563821033UL) + ((((uint64_t)op[2] * 4611714075407134313UL) + ((uint64_t)op[3] * 5307564123049521121UL) + ((uint64_t)op[4] * 9341691724871993531UL) + ((uint64_t)op[5] * 4367203912113899393UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 4367203912113899393UL) + ((uint64_t)op[1] * 10502031650334300791UL) + ((uint64_t)op[2] * 16950787534563821033UL) + ((((uint64_t)op[3] * 4611714075407134313UL) + ((uint64_t)op[4] * 5307564123049521121UL) + ((uint64_t)op[5] * 9341691724871993531UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 9341691724871993531UL) + ((uint64_t)op[1] * 4367203912113899393UL) + ((uint64_t)op[2] * 10502031650334300791UL) + ((uint64_t)op[3] * 16950787534563821033UL) + ((((uint64_t)op[4] * 4611714075407134313UL) + ((uint64_t)op[5] * 5307564123049521121UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 5307564123049521121UL) + ((uint64_t)op[1] * 9341691724871993531UL) + ((uint64_t)op[2] * 4367203912113899393UL) + ((uint64_t)op[3] * 10502031650334300791UL) + ((uint64_t)op[4] * 16950787534563821033UL) + ((uint64_t)op[5] * 9223203694976297354UL);
	tmp_q[5] = ((uint64_t)op[0] * 4611714075407134313UL) + ((uint64_t)op[1] * 5307564123049521121UL) + ((uint64_t)op[2] * 9341691724871993531UL) + ((uint64_t)op[3] * 4367203912113899393UL) + ((uint64_t)op[4] * 10502031650334300791UL) + ((uint64_t)op[5] * 16950787534563821033UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 31787318537L) - ((-((int128)tmp_q[1] * 27442123586L) - ((int128)tmp_q[2] * 102167785120L) - ((int128)tmp_q[3] * 127579608108L) - ((int128)tmp_q[4] * 43197920616L) + ((int128)tmp_q[5] * 6065772677L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 6065772677L) + ((int128)tmp_q[1] * 31787318537L) - ((-((int128)tmp_q[2] * 27442123586L) - ((int128)tmp_q[3] * 102167785120L) - ((int128)tmp_q[4] * 127579608108L) - ((int128)tmp_q[5] * 43197920616L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 43197920616L) + ((int128)tmp_q[1] * 6065772677L) + ((int128)tmp_q[2] * 31787318537L) - ((-((int128)tmp_q[3] * 27442123586L) - ((int128)tmp_q[4] * 102167785120L) - ((int128)tmp_q[5] * 127579608108L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 127579608108L) - ((int128)tmp_q[1] * 43197920616L) + ((int128)tmp_q[2] * 6065772677L) + ((int128)tmp_q[3] * 31787318537L) - ((-((int128)tmp_q[4] * 27442123586L) - ((int128)tmp_q[5] * 102167785120L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 102167785120L) - ((int128)tmp_q[1] * 127579608108L) - ((int128)tmp_q[2] * 43197920616L) + ((int128)tmp_q[3] * 6065772677L) + ((int128)tmp_q[4] * 31787318537L) + ((int128)tmp_q[5] * 164652741516L);
	tmp_zero[5] = -((int128)tmp_q[0] * 27442123586L) - ((int128)tmp_q[1] * 102167785120L) - ((int128)tmp_q[2] * 127579608108L) - ((int128)tmp_q[3] * 43197920616L) + ((int128)tmp_q[4] * 6065772677L) + ((int128)tmp_q[5] * 31787318537L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

