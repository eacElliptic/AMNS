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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12073075291050678859UL) + ((((uint64_t)op[1] * 1414608646426415585UL) + ((uint64_t)op[2] * 5623877815878313452UL) + ((uint64_t)op[3] * 2797753346198364111UL) + ((uint64_t)op[4] * 15303028963962180661UL) + ((uint64_t)op[5] * 10619240348812197128UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 10619240348812197128UL) + ((uint64_t)op[1] * 12073075291050678859UL) + ((((uint64_t)op[2] * 1414608646426415585UL) + ((uint64_t)op[3] * 5623877815878313452UL) + ((uint64_t)op[4] * 2797753346198364111UL) + ((uint64_t)op[5] * 15303028963962180661UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 15303028963962180661UL) + ((uint64_t)op[1] * 10619240348812197128UL) + ((uint64_t)op[2] * 12073075291050678859UL) + ((((uint64_t)op[3] * 1414608646426415585UL) + ((uint64_t)op[4] * 5623877815878313452UL) + ((uint64_t)op[5] * 2797753346198364111UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 2797753346198364111UL) + ((uint64_t)op[1] * 15303028963962180661UL) + ((uint64_t)op[2] * 10619240348812197128UL) + ((uint64_t)op[3] * 12073075291050678859UL) + ((((uint64_t)op[4] * 1414608646426415585UL) + ((uint64_t)op[5] * 5623877815878313452UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 5623877815878313452UL) + ((uint64_t)op[1] * 2797753346198364111UL) + ((uint64_t)op[2] * 15303028963962180661UL) + ((uint64_t)op[3] * 10619240348812197128UL) + ((uint64_t)op[4] * 12073075291050678859UL) + ((uint64_t)op[5] * 12788309488003889276UL);
	tmp_q[5] = ((uint64_t)op[0] * 1414608646426415585UL) + ((uint64_t)op[1] * 5623877815878313452UL) + ((uint64_t)op[2] * 2797753346198364111UL) + ((uint64_t)op[3] * 15303028963962180661UL) + ((uint64_t)op[4] * 10619240348812197128UL) + ((uint64_t)op[5] * 12073075291050678859UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1881760027L) - ((((int128)tmp_q[1] * 1545812359L) + ((int128)tmp_q[2] * 775468337L) - ((int128)tmp_q[3] * 2935111953L) - ((int128)tmp_q[4] * 955820811L) - ((int128)tmp_q[5] * 2153039924L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 2153039924L) - ((int128)tmp_q[1] * 1881760027L) - ((((int128)tmp_q[2] * 1545812359L) + ((int128)tmp_q[3] * 775468337L) - ((int128)tmp_q[4] * 2935111953L) - ((int128)tmp_q[5] * 955820811L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 955820811L) - ((int128)tmp_q[1] * 2153039924L) - ((int128)tmp_q[2] * 1881760027L) - ((((int128)tmp_q[3] * 1545812359L) + ((int128)tmp_q[4] * 775468337L) - ((int128)tmp_q[5] * 2935111953L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 2935111953L) - ((int128)tmp_q[1] * 955820811L) - ((int128)tmp_q[2] * 2153039924L) - ((int128)tmp_q[3] * 1881760027L) - ((((int128)tmp_q[4] * 1545812359L) + ((int128)tmp_q[5] * 775468337L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 775468337L) - ((int128)tmp_q[1] * 2935111953L) - ((int128)tmp_q[2] * 955820811L) - ((int128)tmp_q[3] * 2153039924L) - ((int128)tmp_q[4] * 1881760027L) - ((int128)tmp_q[5] * 6183249436L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1545812359L) + ((int128)tmp_q[1] * 775468337L) - ((int128)tmp_q[2] * 2935111953L) - ((int128)tmp_q[3] * 955820811L) - ((int128)tmp_q[4] * 2153039924L) - ((int128)tmp_q[5] * 1881760027L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

