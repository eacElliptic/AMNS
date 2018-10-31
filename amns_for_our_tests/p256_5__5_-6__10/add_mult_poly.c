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
	tmp_q[0] = ((uint64_t)op[0] * 12461641886127938505UL) + ((((uint64_t)op[1] * 18104258453495208669UL) + ((uint64_t)op[2] * 11403744247532501670UL) + ((uint64_t)op[3] * 2867085358695916190UL) + ((uint64_t)op[4] * 241909599167204440UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 241909599167204440UL) + ((uint64_t)op[1] * 12461641886127938505UL) + ((((uint64_t)op[2] * 18104258453495208669UL) + ((uint64_t)op[3] * 11403744247532501670UL) + ((uint64_t)op[4] * 2867085358695916190UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 2867085358695916190UL) + ((uint64_t)op[1] * 241909599167204440UL) + ((uint64_t)op[2] * 12461641886127938505UL) + ((((uint64_t)op[3] * 18104258453495208669UL) + ((uint64_t)op[4] * 11403744247532501670UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 11403744247532501670UL) + ((uint64_t)op[1] * 2867085358695916190UL) + ((uint64_t)op[2] * 241909599167204440UL) + ((uint64_t)op[3] * 12461641886127938505UL) + ((uint64_t)op[4] * 2054913721286057682UL);
	tmp_q[4] = ((uint64_t)op[0] * 18104258453495208669UL) + ((uint64_t)op[1] * 11403744247532501670UL) + ((uint64_t)op[2] * 2867085358695916190UL) + ((uint64_t)op[3] * 241909599167204440UL) + ((uint64_t)op[4] * 12461641886127938505UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 689296833488895L) - ((((int128)tmp_q[1] * 1472557218767897L) + ((int128)tmp_q[2] * 187727790344100L) - ((int128)tmp_q[3] * 533958199014198L) + ((int128)tmp_q[4] * 240808560660232L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 240808560660232L) + ((int128)tmp_q[1] * 689296833488895L) - ((((int128)tmp_q[2] * 1472557218767897L) + ((int128)tmp_q[3] * 187727790344100L) - ((int128)tmp_q[4] * 533958199014198L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 533958199014198L) + ((int128)tmp_q[1] * 240808560660232L) + ((int128)tmp_q[2] * 689296833488895L) - ((((int128)tmp_q[3] * 1472557218767897L) + ((int128)tmp_q[4] * 187727790344100L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 187727790344100L) - ((int128)tmp_q[1] * 533958199014198L) + ((int128)tmp_q[2] * 240808560660232L) + ((int128)tmp_q[3] * 689296833488895L) - ((int128)tmp_q[4] * 8835343312607382L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1472557218767897L) + ((int128)tmp_q[1] * 187727790344100L) - ((int128)tmp_q[2] * 533958199014198L) + ((int128)tmp_q[3] * 240808560660232L) + ((int128)tmp_q[4] * 689296833488895L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

