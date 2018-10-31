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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15542351457083841051UL) + ((((uint64_t)op[1] * 9045280657387434794UL) + ((uint64_t)op[2] * 9088369879937050265UL) + ((uint64_t)op[3] * 10757925121006547447UL) + ((uint64_t)op[4] * 11463775441668647456UL) + ((uint64_t)op[5] * 12106642174692118734UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 12106642174692118734UL) + ((uint64_t)op[1] * 15542351457083841051UL) + ((((uint64_t)op[2] * 9045280657387434794UL) + ((uint64_t)op[3] * 9088369879937050265UL) + ((uint64_t)op[4] * 10757925121006547447UL) + ((uint64_t)op[5] * 11463775441668647456UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 11463775441668647456UL) + ((uint64_t)op[1] * 12106642174692118734UL) + ((uint64_t)op[2] * 15542351457083841051UL) + ((((uint64_t)op[3] * 9045280657387434794UL) + ((uint64_t)op[4] * 9088369879937050265UL) + ((uint64_t)op[5] * 10757925121006547447UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 10757925121006547447UL) + ((uint64_t)op[1] * 11463775441668647456UL) + ((uint64_t)op[2] * 12106642174692118734UL) + ((uint64_t)op[3] * 15542351457083841051UL) + ((((uint64_t)op[4] * 9045280657387434794UL) + ((uint64_t)op[5] * 9088369879937050265UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 9088369879937050265UL) + ((uint64_t)op[1] * 10757925121006547447UL) + ((uint64_t)op[2] * 11463775441668647456UL) + ((uint64_t)op[3] * 12106642174692118734UL) + ((uint64_t)op[4] * 15542351457083841051UL) + ((uint64_t)op[5] * 1068548276804046084UL);
	tmp_q[5] = ((uint64_t)op[0] * 9045280657387434794UL) + ((uint64_t)op[1] * 9088369879937050265UL) + ((uint64_t)op[2] * 10757925121006547447UL) + ((uint64_t)op[3] * 11463775441668647456UL) + ((uint64_t)op[4] * 12106642174692118734UL) + ((uint64_t)op[5] * 15542351457083841051UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 264873927143L) - ((((int128)tmp_q[1] * 2401250165304L) + ((int128)tmp_q[2] * 4043703569603L) + ((int128)tmp_q[3] * 3718945352109L) - ((int128)tmp_q[4] * 268166205378L) + ((int128)tmp_q[5] * 401365234502L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 401365234502L) + ((int128)tmp_q[1] * 264873927143L) - ((((int128)tmp_q[2] * 2401250165304L) + ((int128)tmp_q[3] * 4043703569603L) + ((int128)tmp_q[4] * 3718945352109L) - ((int128)tmp_q[5] * 268166205378L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 268166205378L) + ((int128)tmp_q[1] * 401365234502L) + ((int128)tmp_q[2] * 264873927143L) - ((((int128)tmp_q[3] * 2401250165304L) + ((int128)tmp_q[4] * 4043703569603L) + ((int128)tmp_q[5] * 3718945352109L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 3718945352109L) - ((int128)tmp_q[1] * 268166205378L) + ((int128)tmp_q[2] * 401365234502L) + ((int128)tmp_q[3] * 264873927143L) - ((((int128)tmp_q[4] * 2401250165304L) + ((int128)tmp_q[5] * 4043703569603L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 4043703569603L) + ((int128)tmp_q[1] * 3718945352109L) - ((int128)tmp_q[2] * 268166205378L) + ((int128)tmp_q[3] * 401365234502L) + ((int128)tmp_q[4] * 264873927143L) - ((int128)tmp_q[5] * 14407500991824L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2401250165304L) + ((int128)tmp_q[1] * 4043703569603L) + ((int128)tmp_q[2] * 3718945352109L) - ((int128)tmp_q[3] * 268166205378L) + ((int128)tmp_q[4] * 401365234502L) + ((int128)tmp_q[5] * 264873927143L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

