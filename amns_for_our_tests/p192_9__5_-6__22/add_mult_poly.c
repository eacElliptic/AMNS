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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6194278786511156767UL) + ((((uint64_t)op[1] * 10856539538578120459UL) + ((uint64_t)op[2] * 731903941403698019UL) + ((uint64_t)op[3] * 1955046748724901458UL) + ((uint64_t)op[4] * 5022757894361737293UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 5022757894361737293UL) + ((uint64_t)op[1] * 6194278786511156767UL) + ((((uint64_t)op[2] * 10856539538578120459UL) + ((uint64_t)op[3] * 731903941403698019UL) + ((uint64_t)op[4] * 1955046748724901458UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 1955046748724901458UL) + ((uint64_t)op[1] * 5022757894361737293UL) + ((uint64_t)op[2] * 6194278786511156767UL) + ((((uint64_t)op[3] * 10856539538578120459UL) + ((uint64_t)op[4] * 731903941403698019UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 731903941403698019UL) + ((uint64_t)op[1] * 1955046748724901458UL) + ((uint64_t)op[2] * 5022757894361737293UL) + ((uint64_t)op[3] * 6194278786511156767UL) + ((uint64_t)op[4] * 8647739063369483710UL);
	tmp_q[4] = ((uint64_t)op[0] * 10856539538578120459UL) + ((uint64_t)op[1] * 731903941403698019UL) + ((uint64_t)op[2] * 1955046748724901458UL) + ((uint64_t)op[3] * 5022757894361737293UL) + ((uint64_t)op[4] * 6194278786511156767UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 46386256119L) - ((-((int128)tmp_q[1] * 149403556728L) - ((int128)tmp_q[2] * 145093111498L) - ((int128)tmp_q[3] * 141089818047L) + ((int128)tmp_q[4] * 31033177383L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 31033177383L) - ((int128)tmp_q[1] * 46386256119L) - ((-((int128)tmp_q[2] * 149403556728L) - ((int128)tmp_q[3] * 145093111498L) - ((int128)tmp_q[4] * 141089818047L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 141089818047L) + ((int128)tmp_q[1] * 31033177383L) - ((int128)tmp_q[2] * 46386256119L) - ((-((int128)tmp_q[3] * 149403556728L) - ((int128)tmp_q[4] * 145093111498L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 145093111498L) - ((int128)tmp_q[1] * 141089818047L) + ((int128)tmp_q[2] * 31033177383L) - ((int128)tmp_q[3] * 46386256119L) + ((int128)tmp_q[4] * 896421340368L);
	tmp_zero[4] = -((int128)tmp_q[0] * 149403556728L) - ((int128)tmp_q[1] * 145093111498L) - ((int128)tmp_q[2] * 141089818047L) + ((int128)tmp_q[3] * 31033177383L) - ((int128)tmp_q[4] * 46386256119L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

