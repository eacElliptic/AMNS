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
	tmp_q[0] = ((uint64_t)op[0] * 13160319514481032351UL) + ((((uint64_t)op[1] * 392094393858631340UL) + ((uint64_t)op[2] * 11854541789370785135UL) + ((uint64_t)op[3] * 5224494717092749485UL) + ((uint64_t)op[4] * 16662393843846482556UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16662393843846482556UL) + ((uint64_t)op[1] * 13160319514481032351UL) + ((((uint64_t)op[2] * 392094393858631340UL) + ((uint64_t)op[3] * 11854541789370785135UL) + ((uint64_t)op[4] * 5224494717092749485UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 5224494717092749485UL) + ((uint64_t)op[1] * 16662393843846482556UL) + ((uint64_t)op[2] * 13160319514481032351UL) + ((((uint64_t)op[3] * 392094393858631340UL) + ((uint64_t)op[4] * 11854541789370785135UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 11854541789370785135UL) + ((uint64_t)op[1] * 5224494717092749485UL) + ((uint64_t)op[2] * 16662393843846482556UL) + ((uint64_t)op[3] * 13160319514481032351UL) + ((uint64_t)op[4] * 16094177710557763576UL);
	tmp_q[4] = ((uint64_t)op[0] * 392094393858631340UL) + ((uint64_t)op[1] * 11854541789370785135UL) + ((uint64_t)op[2] * 5224494717092749485UL) + ((uint64_t)op[3] * 16662393843846482556UL) + ((uint64_t)op[4] * 13160319514481032351UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1150647777520033L) - ((-((int128)tmp_q[1] * 48934812257465L) + ((int128)tmp_q[2] * 989686011641495L) + ((int128)tmp_q[3] * 389113969946255L) - ((int128)tmp_q[4] * 1122565759793588L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 1122565759793588L) + ((int128)tmp_q[1] * 1150647777520033L) - ((-((int128)tmp_q[2] * 48934812257465L) + ((int128)tmp_q[3] * 989686011641495L) + ((int128)tmp_q[4] * 389113969946255L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 389113969946255L) - ((int128)tmp_q[1] * 1122565759793588L) + ((int128)tmp_q[2] * 1150647777520033L) - ((-((int128)tmp_q[3] * 48934812257465L) + ((int128)tmp_q[4] * 989686011641495L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 989686011641495L) + ((int128)tmp_q[1] * 389113969946255L) - ((int128)tmp_q[2] * 1122565759793588L) + ((int128)tmp_q[3] * 1150647777520033L) + ((int128)tmp_q[4] * 293608873544790L);
	tmp_zero[4] = -((int128)tmp_q[0] * 48934812257465L) + ((int128)tmp_q[1] * 989686011641495L) + ((int128)tmp_q[2] * 389113969946255L) - ((int128)tmp_q[3] * 1122565759793588L) + ((int128)tmp_q[4] * 1150647777520033L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

