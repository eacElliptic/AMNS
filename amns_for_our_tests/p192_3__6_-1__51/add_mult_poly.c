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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3806539322401247713UL) + ((((uint64_t)op[1] * 661503553662950216UL) + ((uint64_t)op[2] * 13750813731419070508UL) + ((uint64_t)op[3] * 14876336996198233992UL) + ((uint64_t)op[4] * 13868004239078407898UL) + ((uint64_t)op[5] * 14214833442535283776UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 14214833442535283776UL) + ((uint64_t)op[1] * 3806539322401247713UL) + ((((uint64_t)op[2] * 661503553662950216UL) + ((uint64_t)op[3] * 13750813731419070508UL) + ((uint64_t)op[4] * 14876336996198233992UL) + ((uint64_t)op[5] * 13868004239078407898UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 13868004239078407898UL) + ((uint64_t)op[1] * 14214833442535283776UL) + ((uint64_t)op[2] * 3806539322401247713UL) + ((((uint64_t)op[3] * 661503553662950216UL) + ((uint64_t)op[4] * 13750813731419070508UL) + ((uint64_t)op[5] * 14876336996198233992UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 14876336996198233992UL) + ((uint64_t)op[1] * 13868004239078407898UL) + ((uint64_t)op[2] * 14214833442535283776UL) + ((uint64_t)op[3] * 3806539322401247713UL) + ((((uint64_t)op[4] * 661503553662950216UL) + ((uint64_t)op[5] * 13750813731419070508UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 13750813731419070508UL) + ((uint64_t)op[1] * 14876336996198233992UL) + ((uint64_t)op[2] * 13868004239078407898UL) + ((uint64_t)op[3] * 14214833442535283776UL) + ((uint64_t)op[4] * 3806539322401247713UL) + ((uint64_t)op[5] * 17785240520046601400UL);
	tmp_q[5] = ((uint64_t)op[0] * 661503553662950216UL) + ((uint64_t)op[1] * 13750813731419070508UL) + ((uint64_t)op[2] * 14876336996198233992UL) + ((uint64_t)op[3] * 13868004239078407898UL) + ((uint64_t)op[4] * 14214833442535283776UL) + ((uint64_t)op[5] * 3806539322401247713UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 137235661266855L) - (((int128)tmp_q[1] * 34407772882728L) + ((int128)tmp_q[2] * 25198163024456L) + ((int128)tmp_q[3] * 154909577207144L) + ((int128)tmp_q[4] * 162433824291306L) + ((int128)tmp_q[5] * 120501804324416L));
	tmp_zero[1] = ((int128)tmp_q[0] * 120501804324416L) + ((int128)tmp_q[1] * 137235661266855L) - (((int128)tmp_q[2] * 34407772882728L) + ((int128)tmp_q[3] * 25198163024456L) + ((int128)tmp_q[4] * 154909577207144L) + ((int128)tmp_q[5] * 162433824291306L));
	tmp_zero[2] = ((int128)tmp_q[0] * 162433824291306L) + ((int128)tmp_q[1] * 120501804324416L) + ((int128)tmp_q[2] * 137235661266855L) - (((int128)tmp_q[3] * 34407772882728L) + ((int128)tmp_q[4] * 25198163024456L) + ((int128)tmp_q[5] * 154909577207144L));
	tmp_zero[3] = ((int128)tmp_q[0] * 154909577207144L) + ((int128)tmp_q[1] * 162433824291306L) + ((int128)tmp_q[2] * 120501804324416L) + ((int128)tmp_q[3] * 137235661266855L) - (((int128)tmp_q[4] * 34407772882728L) + ((int128)tmp_q[5] * 25198163024456L));
	tmp_zero[4] = ((int128)tmp_q[0] * 25198163024456L) + ((int128)tmp_q[1] * 154909577207144L) + ((int128)tmp_q[2] * 162433824291306L) + ((int128)tmp_q[3] * 120501804324416L) + ((int128)tmp_q[4] * 137235661266855L) - ((int128)tmp_q[5] * 34407772882728L);
	tmp_zero[5] = ((int128)tmp_q[0] * 34407772882728L) + ((int128)tmp_q[1] * 25198163024456L) + ((int128)tmp_q[2] * 154909577207144L) + ((int128)tmp_q[3] * 162433824291306L) + ((int128)tmp_q[4] * 120501804324416L) + ((int128)tmp_q[5] * 137235661266855L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

