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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5175912124953719232UL) + ((((uint64_t)op[1] * 14283265683245511563UL) + ((uint64_t)op[2] * 2050961555256759040UL) + ((uint64_t)op[3] * 10181126412894597434UL) + ((uint64_t)op[4] * 1525749306751904417UL) + ((uint64_t)op[5] * 56327873731988784UL) + ((uint64_t)op[6] * 12909558677773959221UL) + ((uint64_t)op[7] * 12895539616513512332UL) + ((uint64_t)op[8] * 3477391849392435152UL) + ((uint64_t)op[9] * 17320098306381788466UL) + ((uint64_t)op[10] * 14601107684561910288UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 14601107684561910288UL) + ((uint64_t)op[1] * 5175912124953719232UL) + ((((uint64_t)op[2] * 14283265683245511563UL) + ((uint64_t)op[3] * 2050961555256759040UL) + ((uint64_t)op[4] * 10181126412894597434UL) + ((uint64_t)op[5] * 1525749306751904417UL) + ((uint64_t)op[6] * 56327873731988784UL) + ((uint64_t)op[7] * 12909558677773959221UL) + ((uint64_t)op[8] * 12895539616513512332UL) + ((uint64_t)op[9] * 3477391849392435152UL) + ((uint64_t)op[10] * 17320098306381788466UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 17320098306381788466UL) + ((uint64_t)op[1] * 14601107684561910288UL) + ((uint64_t)op[2] * 5175912124953719232UL) + ((((uint64_t)op[3] * 14283265683245511563UL) + ((uint64_t)op[4] * 2050961555256759040UL) + ((uint64_t)op[5] * 10181126412894597434UL) + ((uint64_t)op[6] * 1525749306751904417UL) + ((uint64_t)op[7] * 56327873731988784UL) + ((uint64_t)op[8] * 12909558677773959221UL) + ((uint64_t)op[9] * 12895539616513512332UL) + ((uint64_t)op[10] * 3477391849392435152UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 3477391849392435152UL) + ((uint64_t)op[1] * 17320098306381788466UL) + ((uint64_t)op[2] * 14601107684561910288UL) + ((uint64_t)op[3] * 5175912124953719232UL) + ((((uint64_t)op[4] * 14283265683245511563UL) + ((uint64_t)op[5] * 2050961555256759040UL) + ((uint64_t)op[6] * 10181126412894597434UL) + ((uint64_t)op[7] * 1525749306751904417UL) + ((uint64_t)op[8] * 56327873731988784UL) + ((uint64_t)op[9] * 12909558677773959221UL) + ((uint64_t)op[10] * 12895539616513512332UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 12895539616513512332UL) + ((uint64_t)op[1] * 3477391849392435152UL) + ((uint64_t)op[2] * 17320098306381788466UL) + ((uint64_t)op[3] * 14601107684561910288UL) + ((uint64_t)op[4] * 5175912124953719232UL) + ((((uint64_t)op[5] * 14283265683245511563UL) + ((uint64_t)op[6] * 2050961555256759040UL) + ((uint64_t)op[7] * 10181126412894597434UL) + ((uint64_t)op[8] * 1525749306751904417UL) + ((uint64_t)op[9] * 56327873731988784UL) + ((uint64_t)op[10] * 12909558677773959221UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 12909558677773959221UL) + ((uint64_t)op[1] * 12895539616513512332UL) + ((uint64_t)op[2] * 3477391849392435152UL) + ((uint64_t)op[3] * 17320098306381788466UL) + ((uint64_t)op[4] * 14601107684561910288UL) + ((uint64_t)op[5] * 5175912124953719232UL) + ((((uint64_t)op[6] * 14283265683245511563UL) + ((uint64_t)op[7] * 2050961555256759040UL) + ((uint64_t)op[8] * 10181126412894597434UL) + ((uint64_t)op[9] * 1525749306751904417UL) + ((uint64_t)op[10] * 56327873731988784UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 56327873731988784UL) + ((uint64_t)op[1] * 12909558677773959221UL) + ((uint64_t)op[2] * 12895539616513512332UL) + ((uint64_t)op[3] * 3477391849392435152UL) + ((uint64_t)op[4] * 17320098306381788466UL) + ((uint64_t)op[5] * 14601107684561910288UL) + ((uint64_t)op[6] * 5175912124953719232UL) + ((((uint64_t)op[7] * 14283265683245511563UL) + ((uint64_t)op[8] * 2050961555256759040UL) + ((uint64_t)op[9] * 10181126412894597434UL) + ((uint64_t)op[10] * 1525749306751904417UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 1525749306751904417UL) + ((uint64_t)op[1] * 56327873731988784UL) + ((uint64_t)op[2] * 12909558677773959221UL) + ((uint64_t)op[3] * 12895539616513512332UL) + ((uint64_t)op[4] * 3477391849392435152UL) + ((uint64_t)op[5] * 17320098306381788466UL) + ((uint64_t)op[6] * 14601107684561910288UL) + ((uint64_t)op[7] * 5175912124953719232UL) + ((((uint64_t)op[8] * 14283265683245511563UL) + ((uint64_t)op[9] * 2050961555256759040UL) + ((uint64_t)op[10] * 10181126412894597434UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 10181126412894597434UL) + ((uint64_t)op[1] * 1525749306751904417UL) + ((uint64_t)op[2] * 56327873731988784UL) + ((uint64_t)op[3] * 12909558677773959221UL) + ((uint64_t)op[4] * 12895539616513512332UL) + ((uint64_t)op[5] * 3477391849392435152UL) + ((uint64_t)op[6] * 17320098306381788466UL) + ((uint64_t)op[7] * 14601107684561910288UL) + ((uint64_t)op[8] * 5175912124953719232UL) + ((((uint64_t)op[9] * 14283265683245511563UL) + ((uint64_t)op[10] * 2050961555256759040UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 2050961555256759040UL) + ((uint64_t)op[1] * 10181126412894597434UL) + ((uint64_t)op[2] * 1525749306751904417UL) + ((uint64_t)op[3] * 56327873731988784UL) + ((uint64_t)op[4] * 12909558677773959221UL) + ((uint64_t)op[5] * 12895539616513512332UL) + ((uint64_t)op[6] * 3477391849392435152UL) + ((uint64_t)op[7] * 17320098306381788466UL) + ((uint64_t)op[8] * 14601107684561910288UL) + ((uint64_t)op[9] * 5175912124953719232UL) + ((uint64_t)op[10] * 7749139414170822861UL);
	tmp_q[10] = ((uint64_t)op[0] * 14283265683245511563UL) + ((uint64_t)op[1] * 2050961555256759040UL) + ((uint64_t)op[2] * 10181126412894597434UL) + ((uint64_t)op[3] * 1525749306751904417UL) + ((uint64_t)op[4] * 56327873731988784UL) + ((uint64_t)op[5] * 12909558677773959221UL) + ((uint64_t)op[6] * 12895539616513512332UL) + ((uint64_t)op[7] * 3477391849392435152UL) + ((uint64_t)op[8] * 17320098306381788466UL) + ((uint64_t)op[9] * 14601107684561910288UL) + ((uint64_t)op[10] * 5175912124953719232UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 31513857712031L) + ((((int128)tmp_q[1] * 51362200382012L) - ((int128)tmp_q[2] * 10041363917776L) + ((int128)tmp_q[3] * 98501407319336L) - ((int128)tmp_q[4] * 30765262908377L) + ((int128)tmp_q[5] * 70159091510571L) - ((int128)tmp_q[6] * 23176772405676L) - ((int128)tmp_q[7] * 105957889171997L) - ((int128)tmp_q[8] * 4454133386583L) - ((int128)tmp_q[9] * 76223330256775L) + ((int128)tmp_q[10] * 85639519913667L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 85639519913667L) + ((int128)tmp_q[1] * 31513857712031L) + ((((int128)tmp_q[2] * 51362200382012L) - ((int128)tmp_q[3] * 10041363917776L) + ((int128)tmp_q[4] * 98501407319336L) - ((int128)tmp_q[5] * 30765262908377L) + ((int128)tmp_q[6] * 70159091510571L) - ((int128)tmp_q[7] * 23176772405676L) - ((int128)tmp_q[8] * 105957889171997L) - ((int128)tmp_q[9] * 4454133386583L) - ((int128)tmp_q[10] * 76223330256775L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 76223330256775L) + ((int128)tmp_q[1] * 85639519913667L) + ((int128)tmp_q[2] * 31513857712031L) + ((((int128)tmp_q[3] * 51362200382012L) - ((int128)tmp_q[4] * 10041363917776L) + ((int128)tmp_q[5] * 98501407319336L) - ((int128)tmp_q[6] * 30765262908377L) + ((int128)tmp_q[7] * 70159091510571L) - ((int128)tmp_q[8] * 23176772405676L) - ((int128)tmp_q[9] * 105957889171997L) - ((int128)tmp_q[10] * 4454133386583L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 4454133386583L) - ((int128)tmp_q[1] * 76223330256775L) + ((int128)tmp_q[2] * 85639519913667L) + ((int128)tmp_q[3] * 31513857712031L) + ((((int128)tmp_q[4] * 51362200382012L) - ((int128)tmp_q[5] * 10041363917776L) + ((int128)tmp_q[6] * 98501407319336L) - ((int128)tmp_q[7] * 30765262908377L) + ((int128)tmp_q[8] * 70159091510571L) - ((int128)tmp_q[9] * 23176772405676L) - ((int128)tmp_q[10] * 105957889171997L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 105957889171997L) - ((int128)tmp_q[1] * 4454133386583L) - ((int128)tmp_q[2] * 76223330256775L) + ((int128)tmp_q[3] * 85639519913667L) + ((int128)tmp_q[4] * 31513857712031L) + ((((int128)tmp_q[5] * 51362200382012L) - ((int128)tmp_q[6] * 10041363917776L) + ((int128)tmp_q[7] * 98501407319336L) - ((int128)tmp_q[8] * 30765262908377L) + ((int128)tmp_q[9] * 70159091510571L) - ((int128)tmp_q[10] * 23176772405676L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 23176772405676L) - ((int128)tmp_q[1] * 105957889171997L) - ((int128)tmp_q[2] * 4454133386583L) - ((int128)tmp_q[3] * 76223330256775L) + ((int128)tmp_q[4] * 85639519913667L) + ((int128)tmp_q[5] * 31513857712031L) + ((((int128)tmp_q[6] * 51362200382012L) - ((int128)tmp_q[7] * 10041363917776L) + ((int128)tmp_q[8] * 98501407319336L) - ((int128)tmp_q[9] * 30765262908377L) + ((int128)tmp_q[10] * 70159091510571L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 70159091510571L) - ((int128)tmp_q[1] * 23176772405676L) - ((int128)tmp_q[2] * 105957889171997L) - ((int128)tmp_q[3] * 4454133386583L) - ((int128)tmp_q[4] * 76223330256775L) + ((int128)tmp_q[5] * 85639519913667L) + ((int128)tmp_q[6] * 31513857712031L) + ((((int128)tmp_q[7] * 51362200382012L) - ((int128)tmp_q[8] * 10041363917776L) + ((int128)tmp_q[9] * 98501407319336L) - ((int128)tmp_q[10] * 30765262908377L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 30765262908377L) + ((int128)tmp_q[1] * 70159091510571L) - ((int128)tmp_q[2] * 23176772405676L) - ((int128)tmp_q[3] * 105957889171997L) - ((int128)tmp_q[4] * 4454133386583L) - ((int128)tmp_q[5] * 76223330256775L) + ((int128)tmp_q[6] * 85639519913667L) + ((int128)tmp_q[7] * 31513857712031L) + ((((int128)tmp_q[8] * 51362200382012L) - ((int128)tmp_q[9] * 10041363917776L) + ((int128)tmp_q[10] * 98501407319336L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 98501407319336L) - ((int128)tmp_q[1] * 30765262908377L) + ((int128)tmp_q[2] * 70159091510571L) - ((int128)tmp_q[3] * 23176772405676L) - ((int128)tmp_q[4] * 105957889171997L) - ((int128)tmp_q[5] * 4454133386583L) - ((int128)tmp_q[6] * 76223330256775L) + ((int128)tmp_q[7] * 85639519913667L) + ((int128)tmp_q[8] * 31513857712031L) + ((((int128)tmp_q[9] * 51362200382012L) - ((int128)tmp_q[10] * 10041363917776L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 10041363917776L) + ((int128)tmp_q[1] * 98501407319336L) - ((int128)tmp_q[2] * 30765262908377L) + ((int128)tmp_q[3] * 70159091510571L) - ((int128)tmp_q[4] * 23176772405676L) - ((int128)tmp_q[5] * 105957889171997L) - ((int128)tmp_q[6] * 4454133386583L) - ((int128)tmp_q[7] * 76223330256775L) + ((int128)tmp_q[8] * 85639519913667L) + ((int128)tmp_q[9] * 31513857712031L) + ((int128)tmp_q[10] * 359535402674084L);
	tmp_zero[10] = ((int128)tmp_q[0] * 51362200382012L) - ((int128)tmp_q[1] * 10041363917776L) + ((int128)tmp_q[2] * 98501407319336L) - ((int128)tmp_q[3] * 30765262908377L) + ((int128)tmp_q[4] * 70159091510571L) - ((int128)tmp_q[5] * 23176772405676L) - ((int128)tmp_q[6] * 105957889171997L) - ((int128)tmp_q[7] * 4454133386583L) - ((int128)tmp_q[8] * 76223330256775L) + ((int128)tmp_q[9] * 85639519913667L) + ((int128)tmp_q[10] * 31513857712031L);

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
}

