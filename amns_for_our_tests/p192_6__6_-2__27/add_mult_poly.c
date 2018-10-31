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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3233446226836327973UL) + ((((uint64_t)op[1] * 7276767163156960643UL) + ((uint64_t)op[2] * 12320858394832778652UL) + ((uint64_t)op[3] * 4406020414241301444UL) + ((uint64_t)op[4] * 1150376364356985771UL) + ((uint64_t)op[5] * 4737894646318670415UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 4737894646318670415UL) + ((uint64_t)op[1] * 3233446226836327973UL) + ((((uint64_t)op[2] * 7276767163156960643UL) + ((uint64_t)op[3] * 12320858394832778652UL) + ((uint64_t)op[4] * 4406020414241301444UL) + ((uint64_t)op[5] * 1150376364356985771UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 1150376364356985771UL) + ((uint64_t)op[1] * 4737894646318670415UL) + ((uint64_t)op[2] * 3233446226836327973UL) + ((((uint64_t)op[3] * 7276767163156960643UL) + ((uint64_t)op[4] * 12320858394832778652UL) + ((uint64_t)op[5] * 4406020414241301444UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 4406020414241301444UL) + ((uint64_t)op[1] * 1150376364356985771UL) + ((uint64_t)op[2] * 4737894646318670415UL) + ((uint64_t)op[3] * 3233446226836327973UL) + ((((uint64_t)op[4] * 7276767163156960643UL) + ((uint64_t)op[5] * 12320858394832778652UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 12320858394832778652UL) + ((uint64_t)op[1] * 4406020414241301444UL) + ((uint64_t)op[2] * 1150376364356985771UL) + ((uint64_t)op[3] * 4737894646318670415UL) + ((uint64_t)op[4] * 3233446226836327973UL) + ((uint64_t)op[5] * 3893209747395630330UL);
	tmp_q[5] = ((uint64_t)op[0] * 7276767163156960643UL) + ((uint64_t)op[1] * 12320858394832778652UL) + ((uint64_t)op[2] * 4406020414241301444UL) + ((uint64_t)op[3] * 1150376364356985771UL) + ((uint64_t)op[4] * 4737894646318670415UL) + ((uint64_t)op[5] * 3233446226836327973UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2309391897L) - ((-((int128)tmp_q[1] * 610241737L) - ((int128)tmp_q[2] * 2164955317L) - ((int128)tmp_q[3] * 2559873693L) + ((int128)tmp_q[4] * 850852826L) - ((int128)tmp_q[5] * 1395892941L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1395892941L) + ((int128)tmp_q[1] * 2309391897L) - ((-((int128)tmp_q[2] * 610241737L) - ((int128)tmp_q[3] * 2164955317L) - ((int128)tmp_q[4] * 2559873693L) + ((int128)tmp_q[5] * 850852826L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 850852826L) - ((int128)tmp_q[1] * 1395892941L) + ((int128)tmp_q[2] * 2309391897L) - ((-((int128)tmp_q[3] * 610241737L) - ((int128)tmp_q[4] * 2164955317L) - ((int128)tmp_q[5] * 2559873693L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2559873693L) + ((int128)tmp_q[1] * 850852826L) - ((int128)tmp_q[2] * 1395892941L) + ((int128)tmp_q[3] * 2309391897L) - ((-((int128)tmp_q[4] * 610241737L) - ((int128)tmp_q[5] * 2164955317L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 2164955317L) - ((int128)tmp_q[1] * 2559873693L) + ((int128)tmp_q[2] * 850852826L) - ((int128)tmp_q[3] * 1395892941L) + ((int128)tmp_q[4] * 2309391897L) + ((int128)tmp_q[5] * 1220483474L);
	tmp_zero[5] = -((int128)tmp_q[0] * 610241737L) - ((int128)tmp_q[1] * 2164955317L) - ((int128)tmp_q[2] * 2559873693L) + ((int128)tmp_q[3] * 850852826L) - ((int128)tmp_q[4] * 1395892941L) + ((int128)tmp_q[5] * 2309391897L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

