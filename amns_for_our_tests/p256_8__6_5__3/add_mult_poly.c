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
	tmp_q[0] = ((uint64_t)op[0] * 14879611323584275976UL) + ((((uint64_t)op[1] * 1220930882643089859UL) + ((uint64_t)op[2] * 5656187835627325077UL) + ((uint64_t)op[3] * 17146775741178687312UL) + ((uint64_t)op[4] * 9530722883650885827UL) + ((uint64_t)op[5] * 726752163517882800UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 726752163517882800UL) + ((uint64_t)op[1] * 14879611323584275976UL) + ((((uint64_t)op[2] * 1220930882643089859UL) + ((uint64_t)op[3] * 5656187835627325077UL) + ((uint64_t)op[4] * 17146775741178687312UL) + ((uint64_t)op[5] * 9530722883650885827UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 9530722883650885827UL) + ((uint64_t)op[1] * 726752163517882800UL) + ((uint64_t)op[2] * 14879611323584275976UL) + ((((uint64_t)op[3] * 1220930882643089859UL) + ((uint64_t)op[4] * 5656187835627325077UL) + ((uint64_t)op[5] * 17146775741178687312UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 17146775741178687312UL) + ((uint64_t)op[1] * 9530722883650885827UL) + ((uint64_t)op[2] * 726752163517882800UL) + ((uint64_t)op[3] * 14879611323584275976UL) + ((((uint64_t)op[4] * 1220930882643089859UL) + ((uint64_t)op[5] * 5656187835627325077UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 5656187835627325077UL) + ((uint64_t)op[1] * 17146775741178687312UL) + ((uint64_t)op[2] * 9530722883650885827UL) + ((uint64_t)op[3] * 726752163517882800UL) + ((uint64_t)op[4] * 14879611323584275976UL) + ((uint64_t)op[5] * 6104654413215449295UL);
	tmp_q[5] = ((uint64_t)op[0] * 1220930882643089859UL) + ((uint64_t)op[1] * 5656187835627325077UL) + ((uint64_t)op[2] * 17146775741178687312UL) + ((uint64_t)op[3] * 9530722883650885827UL) + ((uint64_t)op[4] * 726752163517882800UL) + ((uint64_t)op[5] * 14879611323584275976UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 876752838807L) + ((-((int128)tmp_q[1] * 99548944584L) + ((int128)tmp_q[2] * 2807683380794L) - ((int128)tmp_q[3] * 220781087273L) - ((int128)tmp_q[4] * 4710838390137L) - ((int128)tmp_q[5] * 3582723287218L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 3582723287218L) + ((int128)tmp_q[1] * 876752838807L) + ((-((int128)tmp_q[2] * 99548944584L) + ((int128)tmp_q[3] * 2807683380794L) - ((int128)tmp_q[4] * 220781087273L) - ((int128)tmp_q[5] * 4710838390137L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 4710838390137L) - ((int128)tmp_q[1] * 3582723287218L) + ((int128)tmp_q[2] * 876752838807L) + ((-((int128)tmp_q[3] * 99548944584L) + ((int128)tmp_q[4] * 2807683380794L) - ((int128)tmp_q[5] * 220781087273L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 220781087273L) - ((int128)tmp_q[1] * 4710838390137L) - ((int128)tmp_q[2] * 3582723287218L) + ((int128)tmp_q[3] * 876752838807L) + ((-((int128)tmp_q[4] * 99548944584L) + ((int128)tmp_q[5] * 2807683380794L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 2807683380794L) - ((int128)tmp_q[1] * 220781087273L) - ((int128)tmp_q[2] * 4710838390137L) - ((int128)tmp_q[3] * 3582723287218L) + ((int128)tmp_q[4] * 876752838807L) - ((int128)tmp_q[5] * 497744722920L);
	tmp_zero[5] = -((int128)tmp_q[0] * 99548944584L) + ((int128)tmp_q[1] * 2807683380794L) - ((int128)tmp_q[2] * 220781087273L) - ((int128)tmp_q[3] * 4710838390137L) - ((int128)tmp_q[4] * 3582723287218L) + ((int128)tmp_q[5] * 876752838807L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

