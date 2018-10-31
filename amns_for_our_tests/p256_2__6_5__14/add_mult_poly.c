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
	tmp_q[0] = ((uint64_t)op[0] * 13926312840191840064UL) + ((((uint64_t)op[1] * 9070355932997661586UL) + ((uint64_t)op[2] * 11817312336139722252UL) + ((uint64_t)op[3] * 9060908607150898806UL) + ((uint64_t)op[4] * 13952929374307239951UL) + ((uint64_t)op[5] * 15950768096861196990UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 15950768096861196990UL) + ((uint64_t)op[1] * 13926312840191840064UL) + ((((uint64_t)op[2] * 9070355932997661586UL) + ((uint64_t)op[3] * 11817312336139722252UL) + ((uint64_t)op[4] * 9060908607150898806UL) + ((uint64_t)op[5] * 13952929374307239951UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13952929374307239951UL) + ((uint64_t)op[1] * 15950768096861196990UL) + ((uint64_t)op[2] * 13926312840191840064UL) + ((((uint64_t)op[3] * 9070355932997661586UL) + ((uint64_t)op[4] * 11817312336139722252UL) + ((uint64_t)op[5] * 9060908607150898806UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 9060908607150898806UL) + ((uint64_t)op[1] * 13952929374307239951UL) + ((uint64_t)op[2] * 15950768096861196990UL) + ((uint64_t)op[3] * 13926312840191840064UL) + ((((uint64_t)op[4] * 9070355932997661586UL) + ((uint64_t)op[5] * 11817312336139722252UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 11817312336139722252UL) + ((uint64_t)op[1] * 9060908607150898806UL) + ((uint64_t)op[2] * 13952929374307239951UL) + ((uint64_t)op[3] * 15950768096861196990UL) + ((uint64_t)op[4] * 13926312840191840064UL) + ((uint64_t)op[5] * 8458291517569204698UL);
	tmp_q[5] = ((uint64_t)op[0] * 9070355932997661586UL) + ((uint64_t)op[1] * 11817312336139722252UL) + ((uint64_t)op[2] * 9060908607150898806UL) + ((uint64_t)op[3] * 13952929374307239951UL) + ((uint64_t)op[4] * 15950768096861196990UL) + ((uint64_t)op[5] * 13926312840191840064UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 229357370312L) + ((((int128)tmp_q[1] * 307234144598L) - ((int128)tmp_q[2] * 2783499071607L) - ((int128)tmp_q[3] * 1411755179202L) + ((int128)tmp_q[4] * 5824325447932L) + ((int128)tmp_q[5] * 2842168597802L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 2842168597802L) + ((int128)tmp_q[1] * 229357370312L) + ((((int128)tmp_q[2] * 307234144598L) - ((int128)tmp_q[3] * 2783499071607L) - ((int128)tmp_q[4] * 1411755179202L) + ((int128)tmp_q[5] * 5824325447932L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 5824325447932L) + ((int128)tmp_q[1] * 2842168597802L) + ((int128)tmp_q[2] * 229357370312L) + ((((int128)tmp_q[3] * 307234144598L) - ((int128)tmp_q[4] * 2783499071607L) - ((int128)tmp_q[5] * 1411755179202L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 1411755179202L) + ((int128)tmp_q[1] * 5824325447932L) + ((int128)tmp_q[2] * 2842168597802L) + ((int128)tmp_q[3] * 229357370312L) + ((((int128)tmp_q[4] * 307234144598L) - ((int128)tmp_q[5] * 2783499071607L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 2783499071607L) - ((int128)tmp_q[1] * 1411755179202L) + ((int128)tmp_q[2] * 5824325447932L) + ((int128)tmp_q[3] * 2842168597802L) + ((int128)tmp_q[4] * 229357370312L) + ((int128)tmp_q[5] * 1536170722990L);
	tmp_zero[5] = ((int128)tmp_q[0] * 307234144598L) - ((int128)tmp_q[1] * 2783499071607L) - ((int128)tmp_q[2] * 1411755179202L) + ((int128)tmp_q[3] * 5824325447932L) + ((int128)tmp_q[4] * 2842168597802L) + ((int128)tmp_q[5] * 229357370312L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

