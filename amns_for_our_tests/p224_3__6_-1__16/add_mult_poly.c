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
	tmp_q[0] = ((uint64_t)op[0] * 17429974569012509599UL) + ((((uint64_t)op[1] * 17202853976549227139UL) + ((uint64_t)op[2] * 13293865996900695744UL) + ((uint64_t)op[3] * 7762041576603892835UL) + ((uint64_t)op[4] * 12277096492203653727UL) + ((uint64_t)op[5] * 9005931673764217313UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 9005931673764217313UL) + ((uint64_t)op[1] * 17429974569012509599UL) + ((((uint64_t)op[2] * 17202853976549227139UL) + ((uint64_t)op[3] * 13293865996900695744UL) + ((uint64_t)op[4] * 7762041576603892835UL) + ((uint64_t)op[5] * 12277096492203653727UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 12277096492203653727UL) + ((uint64_t)op[1] * 9005931673764217313UL) + ((uint64_t)op[2] * 17429974569012509599UL) + ((((uint64_t)op[3] * 17202853976549227139UL) + ((uint64_t)op[4] * 13293865996900695744UL) + ((uint64_t)op[5] * 7762041576603892835UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 7762041576603892835UL) + ((uint64_t)op[1] * 12277096492203653727UL) + ((uint64_t)op[2] * 9005931673764217313UL) + ((uint64_t)op[3] * 17429974569012509599UL) + ((((uint64_t)op[4] * 17202853976549227139UL) + ((uint64_t)op[5] * 13293865996900695744UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 13293865996900695744UL) + ((uint64_t)op[1] * 7762041576603892835UL) + ((uint64_t)op[2] * 12277096492203653727UL) + ((uint64_t)op[3] * 9005931673764217313UL) + ((uint64_t)op[4] * 17429974569012509599UL) + ((uint64_t)op[5] * 1243890097160324477UL);
	tmp_q[5] = ((uint64_t)op[0] * 17202853976549227139UL) + ((uint64_t)op[1] * 13293865996900695744UL) + ((uint64_t)op[2] * 7762041576603892835UL) + ((uint64_t)op[3] * 12277096492203653727UL) + ((uint64_t)op[4] * 9005931673764217313UL) + ((uint64_t)op[5] * 17429974569012509599UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 41523551107368041L) - (-((int128)tmp_q[1] * 6548926278699239L) + ((int128)tmp_q[2] * 26984817605288197L) - ((int128)tmp_q[3] * 18976886543795675L) - ((int128)tmp_q[4] * 14538733502079844L) - ((int128)tmp_q[5] * 12427960265096435L));
	tmp_zero[1] = -((int128)tmp_q[0] * 12427960265096435L) - ((int128)tmp_q[1] * 41523551107368041L) - (-((int128)tmp_q[2] * 6548926278699239L) + ((int128)tmp_q[3] * 26984817605288197L) - ((int128)tmp_q[4] * 18976886543795675L) - ((int128)tmp_q[5] * 14538733502079844L));
	tmp_zero[2] = -((int128)tmp_q[0] * 14538733502079844L) - ((int128)tmp_q[1] * 12427960265096435L) - ((int128)tmp_q[2] * 41523551107368041L) - (-((int128)tmp_q[3] * 6548926278699239L) + ((int128)tmp_q[4] * 26984817605288197L) - ((int128)tmp_q[5] * 18976886543795675L));
	tmp_zero[3] = -((int128)tmp_q[0] * 18976886543795675L) - ((int128)tmp_q[1] * 14538733502079844L) - ((int128)tmp_q[2] * 12427960265096435L) - ((int128)tmp_q[3] * 41523551107368041L) - (-((int128)tmp_q[4] * 6548926278699239L) + ((int128)tmp_q[5] * 26984817605288197L));
	tmp_zero[4] = ((int128)tmp_q[0] * 26984817605288197L) - ((int128)tmp_q[1] * 18976886543795675L) - ((int128)tmp_q[2] * 14538733502079844L) - ((int128)tmp_q[3] * 12427960265096435L) - ((int128)tmp_q[4] * 41523551107368041L) + ((int128)tmp_q[5] * 6548926278699239L);
	tmp_zero[5] = -((int128)tmp_q[0] * 6548926278699239L) + ((int128)tmp_q[1] * 26984817605288197L) - ((int128)tmp_q[2] * 18976886543795675L) - ((int128)tmp_q[3] * 14538733502079844L) - ((int128)tmp_q[4] * 12427960265096435L) - ((int128)tmp_q[5] * 41523551107368041L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

