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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8798384236751667831UL) + ((((uint64_t)op[1] * 2375489651270638486UL) + ((uint64_t)op[2] * 8559686193993668353UL) + ((uint64_t)op[3] * 6794698905767241366UL) + ((uint64_t)op[4] * 17979959636563198395UL) + ((uint64_t)op[5] * 7560771998832064037UL) + ((uint64_t)op[6] * 10195881469798428581UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 10195881469798428581UL) + ((uint64_t)op[1] * 8798384236751667831UL) + ((((uint64_t)op[2] * 2375489651270638486UL) + ((uint64_t)op[3] * 8559686193993668353UL) + ((uint64_t)op[4] * 6794698905767241366UL) + ((uint64_t)op[5] * 17979959636563198395UL) + ((uint64_t)op[6] * 7560771998832064037UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 7560771998832064037UL) + ((uint64_t)op[1] * 10195881469798428581UL) + ((uint64_t)op[2] * 8798384236751667831UL) + ((((uint64_t)op[3] * 2375489651270638486UL) + ((uint64_t)op[4] * 8559686193993668353UL) + ((uint64_t)op[5] * 6794698905767241366UL) + ((uint64_t)op[6] * 17979959636563198395UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 17979959636563198395UL) + ((uint64_t)op[1] * 7560771998832064037UL) + ((uint64_t)op[2] * 10195881469798428581UL) + ((uint64_t)op[3] * 8798384236751667831UL) + ((((uint64_t)op[4] * 2375489651270638486UL) + ((uint64_t)op[5] * 8559686193993668353UL) + ((uint64_t)op[6] * 6794698905767241366UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 6794698905767241366UL) + ((uint64_t)op[1] * 17979959636563198395UL) + ((uint64_t)op[2] * 7560771998832064037UL) + ((uint64_t)op[3] * 10195881469798428581UL) + ((uint64_t)op[4] * 8798384236751667831UL) + ((((uint64_t)op[5] * 2375489651270638486UL) + ((uint64_t)op[6] * 8559686193993668353UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 8559686193993668353UL) + ((uint64_t)op[1] * 6794698905767241366UL) + ((uint64_t)op[2] * 17979959636563198395UL) + ((uint64_t)op[3] * 7560771998832064037UL) + ((uint64_t)op[4] * 10195881469798428581UL) + ((uint64_t)op[5] * 8798384236751667831UL) + ((uint64_t)op[6] * 11877448256353192430UL);
	tmp_q[6] = ((uint64_t)op[0] * 2375489651270638486UL) + ((uint64_t)op[1] * 8559686193993668353UL) + ((uint64_t)op[2] * 6794698905767241366UL) + ((uint64_t)op[3] * 17979959636563198395UL) + ((uint64_t)op[4] * 7560771998832064037UL) + ((uint64_t)op[5] * 10195881469798428581UL) + ((uint64_t)op[6] * 8798384236751667831UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 20483968482L) + ((((int128)tmp_q[1] * 8564363847L) - ((int128)tmp_q[2] * 54483180365L) - ((int128)tmp_q[3] * 62581446090L) - ((int128)tmp_q[4] * 29773988374L) + ((int128)tmp_q[5] * 51528725729L) + ((int128)tmp_q[6] * 9167873186L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 9167873186L) + ((int128)tmp_q[1] * 20483968482L) + ((((int128)tmp_q[2] * 8564363847L) - ((int128)tmp_q[3] * 54483180365L) - ((int128)tmp_q[4] * 62581446090L) - ((int128)tmp_q[5] * 29773988374L) + ((int128)tmp_q[6] * 51528725729L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 51528725729L) + ((int128)tmp_q[1] * 9167873186L) + ((int128)tmp_q[2] * 20483968482L) + ((((int128)tmp_q[3] * 8564363847L) - ((int128)tmp_q[4] * 54483180365L) - ((int128)tmp_q[5] * 62581446090L) - ((int128)tmp_q[6] * 29773988374L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 29773988374L) + ((int128)tmp_q[1] * 51528725729L) + ((int128)tmp_q[2] * 9167873186L) + ((int128)tmp_q[3] * 20483968482L) + ((((int128)tmp_q[4] * 8564363847L) - ((int128)tmp_q[5] * 54483180365L) - ((int128)tmp_q[6] * 62581446090L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 62581446090L) - ((int128)tmp_q[1] * 29773988374L) + ((int128)tmp_q[2] * 51528725729L) + ((int128)tmp_q[3] * 9167873186L) + ((int128)tmp_q[4] * 20483968482L) + ((((int128)tmp_q[5] * 8564363847L) - ((int128)tmp_q[6] * 54483180365L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 54483180365L) - ((int128)tmp_q[1] * 62581446090L) - ((int128)tmp_q[2] * 29773988374L) + ((int128)tmp_q[3] * 51528725729L) + ((int128)tmp_q[4] * 9167873186L) + ((int128)tmp_q[5] * 20483968482L) + ((int128)tmp_q[6] * 42821819235L);
	tmp_zero[6] = ((int128)tmp_q[0] * 8564363847L) - ((int128)tmp_q[1] * 54483180365L) - ((int128)tmp_q[2] * 62581446090L) - ((int128)tmp_q[3] * 29773988374L) + ((int128)tmp_q[4] * 51528725729L) + ((int128)tmp_q[5] * 9167873186L) + ((int128)tmp_q[6] * 20483968482L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

