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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) << 1);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) << 1);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16142903425487003975UL) + ((((uint64_t)op[1] * 11593094283000342024UL) + ((uint64_t)op[2] * 14371920072954873536UL) + ((uint64_t)op[3] * 3237228902217768122UL) + ((uint64_t)op[4] * 18215257407017648876UL) + ((uint64_t)op[5] * 16718196318350402300UL) + ((uint64_t)op[6] * 6711872113147943998UL) + ((uint64_t)op[7] * 17800426605072376579UL) + ((uint64_t)op[8] * 7532343272095280190UL) + ((uint64_t)op[9] * 13843665365790255483UL) + ((uint64_t)op[10] * 10947683593237729900UL) + ((uint64_t)op[11] * 9704048954052723790UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 9704048954052723790UL) + ((uint64_t)op[1] * 16142903425487003975UL) + ((((uint64_t)op[2] * 11593094283000342024UL) + ((uint64_t)op[3] * 14371920072954873536UL) + ((uint64_t)op[4] * 3237228902217768122UL) + ((uint64_t)op[5] * 18215257407017648876UL) + ((uint64_t)op[6] * 16718196318350402300UL) + ((uint64_t)op[7] * 6711872113147943998UL) + ((uint64_t)op[8] * 17800426605072376579UL) + ((uint64_t)op[9] * 7532343272095280190UL) + ((uint64_t)op[10] * 13843665365790255483UL) + ((uint64_t)op[11] * 10947683593237729900UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 10947683593237729900UL) + ((uint64_t)op[1] * 9704048954052723790UL) + ((uint64_t)op[2] * 16142903425487003975UL) + ((((uint64_t)op[3] * 11593094283000342024UL) + ((uint64_t)op[4] * 14371920072954873536UL) + ((uint64_t)op[5] * 3237228902217768122UL) + ((uint64_t)op[6] * 18215257407017648876UL) + ((uint64_t)op[7] * 16718196318350402300UL) + ((uint64_t)op[8] * 6711872113147943998UL) + ((uint64_t)op[9] * 17800426605072376579UL) + ((uint64_t)op[10] * 7532343272095280190UL) + ((uint64_t)op[11] * 13843665365790255483UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 13843665365790255483UL) + ((uint64_t)op[1] * 10947683593237729900UL) + ((uint64_t)op[2] * 9704048954052723790UL) + ((uint64_t)op[3] * 16142903425487003975UL) + ((((uint64_t)op[4] * 11593094283000342024UL) + ((uint64_t)op[5] * 14371920072954873536UL) + ((uint64_t)op[6] * 3237228902217768122UL) + ((uint64_t)op[7] * 18215257407017648876UL) + ((uint64_t)op[8] * 16718196318350402300UL) + ((uint64_t)op[9] * 6711872113147943998UL) + ((uint64_t)op[10] * 17800426605072376579UL) + ((uint64_t)op[11] * 7532343272095280190UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 7532343272095280190UL) + ((uint64_t)op[1] * 13843665365790255483UL) + ((uint64_t)op[2] * 10947683593237729900UL) + ((uint64_t)op[3] * 9704048954052723790UL) + ((uint64_t)op[4] * 16142903425487003975UL) + ((((uint64_t)op[5] * 11593094283000342024UL) + ((uint64_t)op[6] * 14371920072954873536UL) + ((uint64_t)op[7] * 3237228902217768122UL) + ((uint64_t)op[8] * 18215257407017648876UL) + ((uint64_t)op[9] * 16718196318350402300UL) + ((uint64_t)op[10] * 6711872113147943998UL) + ((uint64_t)op[11] * 17800426605072376579UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 17800426605072376579UL) + ((uint64_t)op[1] * 7532343272095280190UL) + ((uint64_t)op[2] * 13843665365790255483UL) + ((uint64_t)op[3] * 10947683593237729900UL) + ((uint64_t)op[4] * 9704048954052723790UL) + ((uint64_t)op[5] * 16142903425487003975UL) + ((((uint64_t)op[6] * 11593094283000342024UL) + ((uint64_t)op[7] * 14371920072954873536UL) + ((uint64_t)op[8] * 3237228902217768122UL) + ((uint64_t)op[9] * 18215257407017648876UL) + ((uint64_t)op[10] * 16718196318350402300UL) + ((uint64_t)op[11] * 6711872113147943998UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 6711872113147943998UL) + ((uint64_t)op[1] * 17800426605072376579UL) + ((uint64_t)op[2] * 7532343272095280190UL) + ((uint64_t)op[3] * 13843665365790255483UL) + ((uint64_t)op[4] * 10947683593237729900UL) + ((uint64_t)op[5] * 9704048954052723790UL) + ((uint64_t)op[6] * 16142903425487003975UL) + ((((uint64_t)op[7] * 11593094283000342024UL) + ((uint64_t)op[8] * 14371920072954873536UL) + ((uint64_t)op[9] * 3237228902217768122UL) + ((uint64_t)op[10] * 18215257407017648876UL) + ((uint64_t)op[11] * 16718196318350402300UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 16718196318350402300UL) + ((uint64_t)op[1] * 6711872113147943998UL) + ((uint64_t)op[2] * 17800426605072376579UL) + ((uint64_t)op[3] * 7532343272095280190UL) + ((uint64_t)op[4] * 13843665365790255483UL) + ((uint64_t)op[5] * 10947683593237729900UL) + ((uint64_t)op[6] * 9704048954052723790UL) + ((uint64_t)op[7] * 16142903425487003975UL) + ((((uint64_t)op[8] * 11593094283000342024UL) + ((uint64_t)op[9] * 14371920072954873536UL) + ((uint64_t)op[10] * 3237228902217768122UL) + ((uint64_t)op[11] * 18215257407017648876UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 18215257407017648876UL) + ((uint64_t)op[1] * 16718196318350402300UL) + ((uint64_t)op[2] * 6711872113147943998UL) + ((uint64_t)op[3] * 17800426605072376579UL) + ((uint64_t)op[4] * 7532343272095280190UL) + ((uint64_t)op[5] * 13843665365790255483UL) + ((uint64_t)op[6] * 10947683593237729900UL) + ((uint64_t)op[7] * 9704048954052723790UL) + ((uint64_t)op[8] * 16142903425487003975UL) + ((((uint64_t)op[9] * 11593094283000342024UL) + ((uint64_t)op[10] * 14371920072954873536UL) + ((uint64_t)op[11] * 3237228902217768122UL)) * 18446744073709551614);
	tmp_q[9] = ((uint64_t)op[0] * 3237228902217768122UL) + ((uint64_t)op[1] * 18215257407017648876UL) + ((uint64_t)op[2] * 16718196318350402300UL) + ((uint64_t)op[3] * 6711872113147943998UL) + ((uint64_t)op[4] * 17800426605072376579UL) + ((uint64_t)op[5] * 7532343272095280190UL) + ((uint64_t)op[6] * 13843665365790255483UL) + ((uint64_t)op[7] * 10947683593237729900UL) + ((uint64_t)op[8] * 9704048954052723790UL) + ((uint64_t)op[9] * 16142903425487003975UL) + ((((uint64_t)op[10] * 11593094283000342024UL) + ((uint64_t)op[11] * 14371920072954873536UL)) * 18446744073709551614);
	tmp_q[10] = ((uint64_t)op[0] * 14371920072954873536UL) + ((uint64_t)op[1] * 3237228902217768122UL) + ((uint64_t)op[2] * 18215257407017648876UL) + ((uint64_t)op[3] * 16718196318350402300UL) + ((uint64_t)op[4] * 6711872113147943998UL) + ((uint64_t)op[5] * 17800426605072376579UL) + ((uint64_t)op[6] * 7532343272095280190UL) + ((uint64_t)op[7] * 13843665365790255483UL) + ((uint64_t)op[8] * 10947683593237729900UL) + ((uint64_t)op[9] * 9704048954052723790UL) + ((uint64_t)op[10] * 16142903425487003975UL) + ((uint64_t)op[11] * 13707299581418419184UL);
	tmp_q[11] = ((uint64_t)op[0] * 11593094283000342024UL) + ((uint64_t)op[1] * 14371920072954873536UL) + ((uint64_t)op[2] * 3237228902217768122UL) + ((uint64_t)op[3] * 18215257407017648876UL) + ((uint64_t)op[4] * 16718196318350402300UL) + ((uint64_t)op[5] * 6711872113147943998UL) + ((uint64_t)op[6] * 17800426605072376579UL) + ((uint64_t)op[7] * 7532343272095280190UL) + ((uint64_t)op[8] * 13843665365790255483UL) + ((uint64_t)op[9] * 10947683593237729900UL) + ((uint64_t)op[10] * 9704048954052723790UL) + ((uint64_t)op[11] * 16142903425487003975UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4748544264975L) - ((((int128)tmp_q[1] * 4697466926975L) - ((int128)tmp_q[2] * 1960869871011L) - ((int128)tmp_q[3] * 1076186807945L) - ((int128)tmp_q[4] * 675743137264L) + ((int128)tmp_q[5] * 3777676356630L) - ((int128)tmp_q[6] * 5712553252371L) - ((int128)tmp_q[7] * 1015990806911L) + ((int128)tmp_q[8] * 5652726656622L) - ((int128)tmp_q[9] * 776218139309L) + ((int128)tmp_q[10] * 4096280334680L) + ((int128)tmp_q[11] * 3084521883096L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 3084521883096L) + ((int128)tmp_q[1] * 4748544264975L) - ((((int128)tmp_q[2] * 4697466926975L) - ((int128)tmp_q[3] * 1960869871011L) - ((int128)tmp_q[4] * 1076186807945L) - ((int128)tmp_q[5] * 675743137264L) + ((int128)tmp_q[6] * 3777676356630L) - ((int128)tmp_q[7] * 5712553252371L) - ((int128)tmp_q[8] * 1015990806911L) + ((int128)tmp_q[9] * 5652726656622L) - ((int128)tmp_q[10] * 776218139309L) + ((int128)tmp_q[11] * 4096280334680L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 4096280334680L) + ((int128)tmp_q[1] * 3084521883096L) + ((int128)tmp_q[2] * 4748544264975L) - ((((int128)tmp_q[3] * 4697466926975L) - ((int128)tmp_q[4] * 1960869871011L) - ((int128)tmp_q[5] * 1076186807945L) - ((int128)tmp_q[6] * 675743137264L) + ((int128)tmp_q[7] * 3777676356630L) - ((int128)tmp_q[8] * 5712553252371L) - ((int128)tmp_q[9] * 1015990806911L) + ((int128)tmp_q[10] * 5652726656622L) - ((int128)tmp_q[11] * 776218139309L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 776218139309L) + ((int128)tmp_q[1] * 4096280334680L) + ((int128)tmp_q[2] * 3084521883096L) + ((int128)tmp_q[3] * 4748544264975L) - ((((int128)tmp_q[4] * 4697466926975L) - ((int128)tmp_q[5] * 1960869871011L) - ((int128)tmp_q[6] * 1076186807945L) - ((int128)tmp_q[7] * 675743137264L) + ((int128)tmp_q[8] * 3777676356630L) - ((int128)tmp_q[9] * 5712553252371L) - ((int128)tmp_q[10] * 1015990806911L) + ((int128)tmp_q[11] * 5652726656622L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 5652726656622L) - ((int128)tmp_q[1] * 776218139309L) + ((int128)tmp_q[2] * 4096280334680L) + ((int128)tmp_q[3] * 3084521883096L) + ((int128)tmp_q[4] * 4748544264975L) - ((((int128)tmp_q[5] * 4697466926975L) - ((int128)tmp_q[6] * 1960869871011L) - ((int128)tmp_q[7] * 1076186807945L) - ((int128)tmp_q[8] * 675743137264L) + ((int128)tmp_q[9] * 3777676356630L) - ((int128)tmp_q[10] * 5712553252371L) - ((int128)tmp_q[11] * 1015990806911L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 1015990806911L) + ((int128)tmp_q[1] * 5652726656622L) - ((int128)tmp_q[2] * 776218139309L) + ((int128)tmp_q[3] * 4096280334680L) + ((int128)tmp_q[4] * 3084521883096L) + ((int128)tmp_q[5] * 4748544264975L) - ((((int128)tmp_q[6] * 4697466926975L) - ((int128)tmp_q[7] * 1960869871011L) - ((int128)tmp_q[8] * 1076186807945L) - ((int128)tmp_q[9] * 675743137264L) + ((int128)tmp_q[10] * 3777676356630L) - ((int128)tmp_q[11] * 5712553252371L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 5712553252371L) - ((int128)tmp_q[1] * 1015990806911L) + ((int128)tmp_q[2] * 5652726656622L) - ((int128)tmp_q[3] * 776218139309L) + ((int128)tmp_q[4] * 4096280334680L) + ((int128)tmp_q[5] * 3084521883096L) + ((int128)tmp_q[6] * 4748544264975L) - ((((int128)tmp_q[7] * 4697466926975L) - ((int128)tmp_q[8] * 1960869871011L) - ((int128)tmp_q[9] * 1076186807945L) - ((int128)tmp_q[10] * 675743137264L) + ((int128)tmp_q[11] * 3777676356630L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 3777676356630L) - ((int128)tmp_q[1] * 5712553252371L) - ((int128)tmp_q[2] * 1015990806911L) + ((int128)tmp_q[3] * 5652726656622L) - ((int128)tmp_q[4] * 776218139309L) + ((int128)tmp_q[5] * 4096280334680L) + ((int128)tmp_q[6] * 3084521883096L) + ((int128)tmp_q[7] * 4748544264975L) - ((((int128)tmp_q[8] * 4697466926975L) - ((int128)tmp_q[9] * 1960869871011L) - ((int128)tmp_q[10] * 1076186807945L) - ((int128)tmp_q[11] * 675743137264L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 675743137264L) + ((int128)tmp_q[1] * 3777676356630L) - ((int128)tmp_q[2] * 5712553252371L) - ((int128)tmp_q[3] * 1015990806911L) + ((int128)tmp_q[4] * 5652726656622L) - ((int128)tmp_q[5] * 776218139309L) + ((int128)tmp_q[6] * 4096280334680L) + ((int128)tmp_q[7] * 3084521883096L) + ((int128)tmp_q[8] * 4748544264975L) - ((((int128)tmp_q[9] * 4697466926975L) - ((int128)tmp_q[10] * 1960869871011L) - ((int128)tmp_q[11] * 1076186807945L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 1076186807945L) - ((int128)tmp_q[1] * 675743137264L) + ((int128)tmp_q[2] * 3777676356630L) - ((int128)tmp_q[3] * 5712553252371L) - ((int128)tmp_q[4] * 1015990806911L) + ((int128)tmp_q[5] * 5652726656622L) - ((int128)tmp_q[6] * 776218139309L) + ((int128)tmp_q[7] * 4096280334680L) + ((int128)tmp_q[8] * 3084521883096L) + ((int128)tmp_q[9] * 4748544264975L) - ((((int128)tmp_q[10] * 4697466926975L) - ((int128)tmp_q[11] * 1960869871011L)) * 2);
	tmp_zero[10] = -((int128)tmp_q[0] * 1960869871011L) - ((int128)tmp_q[1] * 1076186807945L) - ((int128)tmp_q[2] * 675743137264L) + ((int128)tmp_q[3] * 3777676356630L) - ((int128)tmp_q[4] * 5712553252371L) - ((int128)tmp_q[5] * 1015990806911L) + ((int128)tmp_q[6] * 5652726656622L) - ((int128)tmp_q[7] * 776218139309L) + ((int128)tmp_q[8] * 4096280334680L) + ((int128)tmp_q[9] * 3084521883096L) + ((int128)tmp_q[10] * 4748544264975L) - ((int128)tmp_q[11] * 9394933853950L);
	tmp_zero[11] = ((int128)tmp_q[0] * 4697466926975L) - ((int128)tmp_q[1] * 1960869871011L) - ((int128)tmp_q[2] * 1076186807945L) - ((int128)tmp_q[3] * 675743137264L) + ((int128)tmp_q[4] * 3777676356630L) - ((int128)tmp_q[5] * 5712553252371L) - ((int128)tmp_q[6] * 1015990806911L) + ((int128)tmp_q[7] * 5652726656622L) - ((int128)tmp_q[8] * 776218139309L) + ((int128)tmp_q[9] * 4096280334680L) + ((int128)tmp_q[10] * 3084521883096L) + ((int128)tmp_q[11] * 4748544264975L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

