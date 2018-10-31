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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 44236387551597973UL) + ((((uint64_t)op[1] * 5856732181630200620UL) + ((uint64_t)op[2] * 7624203822849195209UL) + ((uint64_t)op[3] * 4211459956208755992UL) + ((uint64_t)op[4] * 799726287998118787UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 799726287998118787UL) + ((uint64_t)op[1] * 44236387551597973UL) + ((((uint64_t)op[2] * 5856732181630200620UL) + ((uint64_t)op[3] * 7624203822849195209UL) + ((uint64_t)op[4] * 4211459956208755992UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 4211459956208755992UL) + ((uint64_t)op[1] * 799726287998118787UL) + ((uint64_t)op[2] * 44236387551597973UL) + ((((uint64_t)op[3] * 5856732181630200620UL) + ((uint64_t)op[4] * 7624203822849195209UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7624203822849195209UL) + ((uint64_t)op[1] * 4211459956208755992UL) + ((uint64_t)op[2] * 799726287998118787UL) + ((uint64_t)op[3] * 44236387551597973UL) + ((uint64_t)op[4] * 876547528818949756UL);
	tmp_q[4] = ((uint64_t)op[0] * 5856732181630200620UL) + ((uint64_t)op[1] * 7624203822849195209UL) + ((uint64_t)op[2] * 4211459956208755992UL) + ((uint64_t)op[3] * 799726287998118787UL) + ((uint64_t)op[4] * 44236387551597973UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 995331113204668L) - ((((int128)tmp_q[1] * 1547808945463230L) + ((int128)tmp_q[2] * 1345095180770909L) - ((int128)tmp_q[3] * 399176888929967L) - ((int128)tmp_q[4] * 42667340589037L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 42667340589037L) + ((int128)tmp_q[1] * 995331113204668L) - ((((int128)tmp_q[2] * 1547808945463230L) + ((int128)tmp_q[3] * 1345095180770909L) - ((int128)tmp_q[4] * 399176888929967L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 399176888929967L) - ((int128)tmp_q[1] * 42667340589037L) + ((int128)tmp_q[2] * 995331113204668L) - ((((int128)tmp_q[3] * 1547808945463230L) + ((int128)tmp_q[4] * 1345095180770909L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1345095180770909L) - ((int128)tmp_q[1] * 399176888929967L) - ((int128)tmp_q[2] * 42667340589037L) + ((int128)tmp_q[3] * 995331113204668L) - ((int128)tmp_q[4] * 4643426836389690L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1547808945463230L) + ((int128)tmp_q[1] * 1345095180770909L) - ((int128)tmp_q[2] * 399176888929967L) - ((int128)tmp_q[3] * 42667340589037L) + ((int128)tmp_q[4] * 995331113204668L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

