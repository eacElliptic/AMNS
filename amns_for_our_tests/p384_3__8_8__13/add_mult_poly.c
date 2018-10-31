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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 4);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15040021007212759643UL) + ((((uint64_t)op[1] * 16627617136122254772UL) + ((uint64_t)op[2] * 9318833035200925880UL) + ((uint64_t)op[3] * 3796260209979470480UL) + ((uint64_t)op[4] * 789377470070392001UL) + ((uint64_t)op[5] * 12451013210021303074UL) + ((uint64_t)op[6] * 13586668815379307294UL) + ((uint64_t)op[7] * 981403640068792724UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 981403640068792724UL) + ((uint64_t)op[1] * 15040021007212759643UL) + ((((uint64_t)op[2] * 16627617136122254772UL) + ((uint64_t)op[3] * 9318833035200925880UL) + ((uint64_t)op[4] * 3796260209979470480UL) + ((uint64_t)op[5] * 789377470070392001UL) + ((uint64_t)op[6] * 12451013210021303074UL) + ((uint64_t)op[7] * 13586668815379307294UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 13586668815379307294UL) + ((uint64_t)op[1] * 981403640068792724UL) + ((uint64_t)op[2] * 15040021007212759643UL) + ((((uint64_t)op[3] * 16627617136122254772UL) + ((uint64_t)op[4] * 9318833035200925880UL) + ((uint64_t)op[5] * 3796260209979470480UL) + ((uint64_t)op[6] * 789377470070392001UL) + ((uint64_t)op[7] * 12451013210021303074UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 12451013210021303074UL) + ((uint64_t)op[1] * 13586668815379307294UL) + ((uint64_t)op[2] * 981403640068792724UL) + ((uint64_t)op[3] * 15040021007212759643UL) + ((((uint64_t)op[4] * 16627617136122254772UL) + ((uint64_t)op[5] * 9318833035200925880UL) + ((uint64_t)op[6] * 3796260209979470480UL) + ((uint64_t)op[7] * 789377470070392001UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 789377470070392001UL) + ((uint64_t)op[1] * 12451013210021303074UL) + ((uint64_t)op[2] * 13586668815379307294UL) + ((uint64_t)op[3] * 981403640068792724UL) + ((uint64_t)op[4] * 15040021007212759643UL) + ((((uint64_t)op[5] * 16627617136122254772UL) + ((uint64_t)op[6] * 9318833035200925880UL) + ((uint64_t)op[7] * 3796260209979470480UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 3796260209979470480UL) + ((uint64_t)op[1] * 789377470070392001UL) + ((uint64_t)op[2] * 12451013210021303074UL) + ((uint64_t)op[3] * 13586668815379307294UL) + ((uint64_t)op[4] * 981403640068792724UL) + ((uint64_t)op[5] * 15040021007212759643UL) + ((((uint64_t)op[6] * 16627617136122254772UL) + ((uint64_t)op[7] * 9318833035200925880UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 9318833035200925880UL) + ((uint64_t)op[1] * 3796260209979470480UL) + ((uint64_t)op[2] * 789377470070392001UL) + ((uint64_t)op[3] * 12451013210021303074UL) + ((uint64_t)op[4] * 13586668815379307294UL) + ((uint64_t)op[5] * 981403640068792724UL) + ((uint64_t)op[6] * 15040021007212759643UL) + ((uint64_t)op[7] * 3893728573011176864UL);
	tmp_q[7] = ((uint64_t)op[0] * 16627617136122254772UL) + ((uint64_t)op[1] * 9318833035200925880UL) + ((uint64_t)op[2] * 3796260209979470480UL) + ((uint64_t)op[3] * 789377470070392001UL) + ((uint64_t)op[4] * 12451013210021303074UL) + ((uint64_t)op[5] * 13586668815379307294UL) + ((uint64_t)op[6] * 981403640068792724UL) + ((uint64_t)op[7] * 15040021007212759643UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 53068424624309L) + ((-((int128)tmp_q[1] * 85844991654496L) - ((int128)tmp_q[2] * 130876286931136L) + ((int128)tmp_q[3] * 26316620257664L) - ((int128)tmp_q[4] * 42699716747851L) + ((int128)tmp_q[5] * 212621329064210L) + ((int128)tmp_q[6] * 124618851440078L) - ((int128)tmp_q[7] * 10947908331116L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 10947908331116L) + ((int128)tmp_q[1] * 53068424624309L) + ((-((int128)tmp_q[2] * 85844991654496L) - ((int128)tmp_q[3] * 130876286931136L) + ((int128)tmp_q[4] * 26316620257664L) - ((int128)tmp_q[5] * 42699716747851L) + ((int128)tmp_q[6] * 212621329064210L) + ((int128)tmp_q[7] * 124618851440078L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 124618851440078L) - ((int128)tmp_q[1] * 10947908331116L) + ((int128)tmp_q[2] * 53068424624309L) + ((-((int128)tmp_q[3] * 85844991654496L) - ((int128)tmp_q[4] * 130876286931136L) + ((int128)tmp_q[5] * 26316620257664L) - ((int128)tmp_q[6] * 42699716747851L) + ((int128)tmp_q[7] * 212621329064210L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 212621329064210L) + ((int128)tmp_q[1] * 124618851440078L) - ((int128)tmp_q[2] * 10947908331116L) + ((int128)tmp_q[3] * 53068424624309L) + ((-((int128)tmp_q[4] * 85844991654496L) - ((int128)tmp_q[5] * 130876286931136L) + ((int128)tmp_q[6] * 26316620257664L) - ((int128)tmp_q[7] * 42699716747851L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 42699716747851L) + ((int128)tmp_q[1] * 212621329064210L) + ((int128)tmp_q[2] * 124618851440078L) - ((int128)tmp_q[3] * 10947908331116L) + ((int128)tmp_q[4] * 53068424624309L) + ((-((int128)tmp_q[5] * 85844991654496L) - ((int128)tmp_q[6] * 130876286931136L) + ((int128)tmp_q[7] * 26316620257664L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 26316620257664L) - ((int128)tmp_q[1] * 42699716747851L) + ((int128)tmp_q[2] * 212621329064210L) + ((int128)tmp_q[3] * 124618851440078L) - ((int128)tmp_q[4] * 10947908331116L) + ((int128)tmp_q[5] * 53068424624309L) + ((-((int128)tmp_q[6] * 85844991654496L) - ((int128)tmp_q[7] * 130876286931136L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 130876286931136L) + ((int128)tmp_q[1] * 26316620257664L) - ((int128)tmp_q[2] * 42699716747851L) + ((int128)tmp_q[3] * 212621329064210L) + ((int128)tmp_q[4] * 124618851440078L) - ((int128)tmp_q[5] * 10947908331116L) + ((int128)tmp_q[6] * 53068424624309L) - ((int128)tmp_q[7] * 686759933235968L);
	tmp_zero[7] = -((int128)tmp_q[0] * 85844991654496L) - ((int128)tmp_q[1] * 130876286931136L) + ((int128)tmp_q[2] * 26316620257664L) - ((int128)tmp_q[3] * 42699716747851L) + ((int128)tmp_q[4] * 212621329064210L) + ((int128)tmp_q[5] * 124618851440078L) - ((int128)tmp_q[6] * 10947908331116L) + ((int128)tmp_q[7] * 53068424624309L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

