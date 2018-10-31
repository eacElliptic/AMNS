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
	tmp_q[0] = ((uint64_t)op[0] * 1974998174148013177UL) + ((((uint64_t)op[1] * 17238638189948140744UL) + ((uint64_t)op[2] * 12338102725520629155UL) + ((uint64_t)op[3] * 18049035391935895102UL) + ((uint64_t)op[4] * 16524898227364016033UL) + ((uint64_t)op[5] * 11866916032701353744UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 11866916032701353744UL) + ((uint64_t)op[1] * 1974998174148013177UL) + ((((uint64_t)op[2] * 17238638189948140744UL) + ((uint64_t)op[3] * 12338102725520629155UL) + ((uint64_t)op[4] * 18049035391935895102UL) + ((uint64_t)op[5] * 16524898227364016033UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 16524898227364016033UL) + ((uint64_t)op[1] * 11866916032701353744UL) + ((uint64_t)op[2] * 1974998174148013177UL) + ((((uint64_t)op[3] * 17238638189948140744UL) + ((uint64_t)op[4] * 12338102725520629155UL) + ((uint64_t)op[5] * 18049035391935895102UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 18049035391935895102UL) + ((uint64_t)op[1] * 16524898227364016033UL) + ((uint64_t)op[2] * 11866916032701353744UL) + ((uint64_t)op[3] * 1974998174148013177UL) + ((((uint64_t)op[4] * 17238638189948140744UL) + ((uint64_t)op[5] * 12338102725520629155UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 12338102725520629155UL) + ((uint64_t)op[1] * 18049035391935895102UL) + ((uint64_t)op[2] * 16524898227364016033UL) + ((uint64_t)op[3] * 11866916032701353744UL) + ((uint64_t)op[4] * 1974998174148013177UL) + ((uint64_t)op[5] * 4832423535045643488UL);
	tmp_q[5] = ((uint64_t)op[0] * 17238638189948140744UL) + ((uint64_t)op[1] * 12338102725520629155UL) + ((uint64_t)op[2] * 18049035391935895102UL) + ((uint64_t)op[3] * 16524898227364016033UL) + ((uint64_t)op[4] * 11866916032701353744UL) + ((uint64_t)op[5] * 1974998174148013177UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1673135733L) - ((-((int128)tmp_q[1] * 910613028L) - ((int128)tmp_q[2] * 648280534L) - ((int128)tmp_q[3] * 1400502130L) - ((int128)tmp_q[4] * 348876523L) - ((int128)tmp_q[5] * 97522824L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 97522824L) - ((int128)tmp_q[1] * 1673135733L) - ((-((int128)tmp_q[2] * 910613028L) - ((int128)tmp_q[3] * 648280534L) - ((int128)tmp_q[4] * 1400502130L) - ((int128)tmp_q[5] * 348876523L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 348876523L) - ((int128)tmp_q[1] * 97522824L) - ((int128)tmp_q[2] * 1673135733L) - ((-((int128)tmp_q[3] * 910613028L) - ((int128)tmp_q[4] * 648280534L) - ((int128)tmp_q[5] * 1400502130L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 1400502130L) - ((int128)tmp_q[1] * 348876523L) - ((int128)tmp_q[2] * 97522824L) - ((int128)tmp_q[3] * 1673135733L) - ((-((int128)tmp_q[4] * 910613028L) - ((int128)tmp_q[5] * 648280534L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 648280534L) - ((int128)tmp_q[1] * 1400502130L) - ((int128)tmp_q[2] * 348876523L) - ((int128)tmp_q[3] * 97522824L) - ((int128)tmp_q[4] * 1673135733L) + ((int128)tmp_q[5] * 3642452112L);
	tmp_zero[5] = -((int128)tmp_q[0] * 910613028L) - ((int128)tmp_q[1] * 648280534L) - ((int128)tmp_q[2] * 1400502130L) - ((int128)tmp_q[3] * 348876523L) - ((int128)tmp_q[4] * 97522824L) - ((int128)tmp_q[5] * 1673135733L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

