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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15074114780743224916UL) + ((((uint64_t)op[1] * 17456093136004733486UL) + ((uint64_t)op[2] * 11666957055543619038UL) + ((uint64_t)op[3] * 2978177403483349772UL) + ((uint64_t)op[4] * 18132958339317484171UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 18132958339317484171UL) + ((uint64_t)op[1] * 15074114780743224916UL) + ((((uint64_t)op[2] * 17456093136004733486UL) + ((uint64_t)op[3] * 11666957055543619038UL) + ((uint64_t)op[4] * 2978177403483349772UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 2978177403483349772UL) + ((uint64_t)op[1] * 18132958339317484171UL) + ((uint64_t)op[2] * 15074114780743224916UL) + ((((uint64_t)op[3] * 17456093136004733486UL) + ((uint64_t)op[4] * 11666957055543619038UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 11666957055543619038UL) + ((uint64_t)op[1] * 2978177403483349772UL) + ((uint64_t)op[2] * 18132958339317484171UL) + ((uint64_t)op[3] * 15074114780743224916UL) + ((uint64_t)op[4] * 13493489385185460966UL);
	tmp_q[4] = ((uint64_t)op[0] * 17456093136004733486UL) + ((uint64_t)op[1] * 11666957055543619038UL) + ((uint64_t)op[2] * 2978177403483349772UL) + ((uint64_t)op[3] * 18132958339317484171UL) + ((uint64_t)op[4] * 15074114780743224916UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 7752903738744L) + ((((int128)tmp_q[1] * 6549187710625L) + ((int128)tmp_q[2] * 10636706045824L) - ((int128)tmp_q[3] * 1088339713914L) + ((int128)tmp_q[4] * 990307140406L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 990307140406L) + ((int128)tmp_q[1] * 7752903738744L) + ((((int128)tmp_q[2] * 6549187710625L) + ((int128)tmp_q[3] * 10636706045824L) - ((int128)tmp_q[4] * 1088339713914L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 1088339713914L) + ((int128)tmp_q[1] * 990307140406L) + ((int128)tmp_q[2] * 7752903738744L) + ((((int128)tmp_q[3] * 6549187710625L) + ((int128)tmp_q[4] * 10636706045824L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 10636706045824L) - ((int128)tmp_q[1] * 1088339713914L) + ((int128)tmp_q[2] * 990307140406L) + ((int128)tmp_q[3] * 7752903738744L) + ((int128)tmp_q[4] * 32745938553125L);
	tmp_zero[4] = ((int128)tmp_q[0] * 6549187710625L) + ((int128)tmp_q[1] * 10636706045824L) - ((int128)tmp_q[2] * 1088339713914L) + ((int128)tmp_q[3] * 990307140406L) + ((int128)tmp_q[4] * 7752903738744L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

