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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1684600514116274903UL) + ((((uint64_t)op[1] * 13144262246046298999UL) + ((uint64_t)op[2] * 9824216774462399520UL) + ((uint64_t)op[3] * 12527445394453431098UL) + ((uint64_t)op[4] * 5247405719592778160UL) + ((uint64_t)op[5] * 13487228194178435258UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 13487228194178435258UL) + ((uint64_t)op[1] * 1684600514116274903UL) + ((((uint64_t)op[2] * 13144262246046298999UL) + ((uint64_t)op[3] * 9824216774462399520UL) + ((uint64_t)op[4] * 12527445394453431098UL) + ((uint64_t)op[5] * 5247405719592778160UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 5247405719592778160UL) + ((uint64_t)op[1] * 13487228194178435258UL) + ((uint64_t)op[2] * 1684600514116274903UL) + ((((uint64_t)op[3] * 13144262246046298999UL) + ((uint64_t)op[4] * 9824216774462399520UL) + ((uint64_t)op[5] * 12527445394453431098UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 12527445394453431098UL) + ((uint64_t)op[1] * 5247405719592778160UL) + ((uint64_t)op[2] * 13487228194178435258UL) + ((uint64_t)op[3] * 1684600514116274903UL) + ((((uint64_t)op[4] * 13144262246046298999UL) + ((uint64_t)op[5] * 9824216774462399520UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 9824216774462399520UL) + ((uint64_t)op[1] * 12527445394453431098UL) + ((uint64_t)op[2] * 5247405719592778160UL) + ((uint64_t)op[3] * 13487228194178435258UL) + ((uint64_t)op[4] * 1684600514116274903UL) + ((uint64_t)op[5] * 7841780418383046382UL);
	tmp_q[5] = ((uint64_t)op[0] * 13144262246046298999UL) + ((uint64_t)op[1] * 9824216774462399520UL) + ((uint64_t)op[2] * 12527445394453431098UL) + ((uint64_t)op[3] * 5247405719592778160UL) + ((uint64_t)op[4] * 13487228194178435258UL) + ((uint64_t)op[5] * 1684600514116274903UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4339348743L) + ((((int128)tmp_q[1] * 387713812979L) - ((int128)tmp_q[2] * 3856976200054L) - ((int128)tmp_q[3] * 3426224242130L) - ((int128)tmp_q[4] * 3769084595884L) - ((int128)tmp_q[5] * 562208982758L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 562208982758L) - ((int128)tmp_q[1] * 4339348743L) + ((((int128)tmp_q[2] * 387713812979L) - ((int128)tmp_q[3] * 3856976200054L) - ((int128)tmp_q[4] * 3426224242130L) - ((int128)tmp_q[5] * 3769084595884L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 3769084595884L) - ((int128)tmp_q[1] * 562208982758L) - ((int128)tmp_q[2] * 4339348743L) + ((((int128)tmp_q[3] * 387713812979L) - ((int128)tmp_q[4] * 3856976200054L) - ((int128)tmp_q[5] * 3426224242130L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 3426224242130L) - ((int128)tmp_q[1] * 3769084595884L) - ((int128)tmp_q[2] * 562208982758L) - ((int128)tmp_q[3] * 4339348743L) + ((((int128)tmp_q[4] * 387713812979L) - ((int128)tmp_q[5] * 3856976200054L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 3856976200054L) - ((int128)tmp_q[1] * 3426224242130L) - ((int128)tmp_q[2] * 3769084595884L) - ((int128)tmp_q[3] * 562208982758L) - ((int128)tmp_q[4] * 4339348743L) + ((int128)tmp_q[5] * 775427625958L);
	tmp_zero[5] = ((int128)tmp_q[0] * 387713812979L) - ((int128)tmp_q[1] * 3856976200054L) - ((int128)tmp_q[2] * 3426224242130L) - ((int128)tmp_q[3] * 3769084595884L) - ((int128)tmp_q[4] * 562208982758L) - ((int128)tmp_q[5] * 4339348743L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

