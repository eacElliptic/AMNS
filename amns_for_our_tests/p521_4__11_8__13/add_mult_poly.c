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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15317309942447422085UL) + ((((uint64_t)op[1] * 7896789175154933191UL) + ((uint64_t)op[2] * 18077693477038659959UL) + ((uint64_t)op[3] * 9115319693145717106UL) + ((uint64_t)op[4] * 14877187536627791829UL) + ((uint64_t)op[5] * 3192455549729999471UL) + ((uint64_t)op[6] * 9424587271780295328UL) + ((uint64_t)op[7] * 9354509106931889212UL) + ((uint64_t)op[8] * 15382615613571648623UL) + ((uint64_t)op[9] * 1983579656587025256UL) + ((uint64_t)op[10] * 12755306781996716162UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 12755306781996716162UL) + ((uint64_t)op[1] * 15317309942447422085UL) + ((((uint64_t)op[2] * 7896789175154933191UL) + ((uint64_t)op[3] * 18077693477038659959UL) + ((uint64_t)op[4] * 9115319693145717106UL) + ((uint64_t)op[5] * 14877187536627791829UL) + ((uint64_t)op[6] * 3192455549729999471UL) + ((uint64_t)op[7] * 9424587271780295328UL) + ((uint64_t)op[8] * 9354509106931889212UL) + ((uint64_t)op[9] * 15382615613571648623UL) + ((uint64_t)op[10] * 1983579656587025256UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 1983579656587025256UL) + ((uint64_t)op[1] * 12755306781996716162UL) + ((uint64_t)op[2] * 15317309942447422085UL) + ((((uint64_t)op[3] * 7896789175154933191UL) + ((uint64_t)op[4] * 18077693477038659959UL) + ((uint64_t)op[5] * 9115319693145717106UL) + ((uint64_t)op[6] * 14877187536627791829UL) + ((uint64_t)op[7] * 3192455549729999471UL) + ((uint64_t)op[8] * 9424587271780295328UL) + ((uint64_t)op[9] * 9354509106931889212UL) + ((uint64_t)op[10] * 15382615613571648623UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 15382615613571648623UL) + ((uint64_t)op[1] * 1983579656587025256UL) + ((uint64_t)op[2] * 12755306781996716162UL) + ((uint64_t)op[3] * 15317309942447422085UL) + ((((uint64_t)op[4] * 7896789175154933191UL) + ((uint64_t)op[5] * 18077693477038659959UL) + ((uint64_t)op[6] * 9115319693145717106UL) + ((uint64_t)op[7] * 14877187536627791829UL) + ((uint64_t)op[8] * 3192455549729999471UL) + ((uint64_t)op[9] * 9424587271780295328UL) + ((uint64_t)op[10] * 9354509106931889212UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 9354509106931889212UL) + ((uint64_t)op[1] * 15382615613571648623UL) + ((uint64_t)op[2] * 1983579656587025256UL) + ((uint64_t)op[3] * 12755306781996716162UL) + ((uint64_t)op[4] * 15317309942447422085UL) + ((((uint64_t)op[5] * 7896789175154933191UL) + ((uint64_t)op[6] * 18077693477038659959UL) + ((uint64_t)op[7] * 9115319693145717106UL) + ((uint64_t)op[8] * 14877187536627791829UL) + ((uint64_t)op[9] * 3192455549729999471UL) + ((uint64_t)op[10] * 9424587271780295328UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 9424587271780295328UL) + ((uint64_t)op[1] * 9354509106931889212UL) + ((uint64_t)op[2] * 15382615613571648623UL) + ((uint64_t)op[3] * 1983579656587025256UL) + ((uint64_t)op[4] * 12755306781996716162UL) + ((uint64_t)op[5] * 15317309942447422085UL) + ((((uint64_t)op[6] * 7896789175154933191UL) + ((uint64_t)op[7] * 18077693477038659959UL) + ((uint64_t)op[8] * 9115319693145717106UL) + ((uint64_t)op[9] * 14877187536627791829UL) + ((uint64_t)op[10] * 3192455549729999471UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 3192455549729999471UL) + ((uint64_t)op[1] * 9424587271780295328UL) + ((uint64_t)op[2] * 9354509106931889212UL) + ((uint64_t)op[3] * 15382615613571648623UL) + ((uint64_t)op[4] * 1983579656587025256UL) + ((uint64_t)op[5] * 12755306781996716162UL) + ((uint64_t)op[6] * 15317309942447422085UL) + ((((uint64_t)op[7] * 7896789175154933191UL) + ((uint64_t)op[8] * 18077693477038659959UL) + ((uint64_t)op[9] * 9115319693145717106UL) + ((uint64_t)op[10] * 14877187536627791829UL)) * 8);
	tmp_q[7] = ((uint64_t)op[0] * 14877187536627791829UL) + ((uint64_t)op[1] * 3192455549729999471UL) + ((uint64_t)op[2] * 9424587271780295328UL) + ((uint64_t)op[3] * 9354509106931889212UL) + ((uint64_t)op[4] * 15382615613571648623UL) + ((uint64_t)op[5] * 1983579656587025256UL) + ((uint64_t)op[6] * 12755306781996716162UL) + ((uint64_t)op[7] * 15317309942447422085UL) + ((((uint64_t)op[8] * 7896789175154933191UL) + ((uint64_t)op[9] * 18077693477038659959UL) + ((uint64_t)op[10] * 9115319693145717106UL)) * 8);
	tmp_q[8] = ((uint64_t)op[0] * 9115319693145717106UL) + ((uint64_t)op[1] * 14877187536627791829UL) + ((uint64_t)op[2] * 3192455549729999471UL) + ((uint64_t)op[3] * 9424587271780295328UL) + ((uint64_t)op[4] * 9354509106931889212UL) + ((uint64_t)op[5] * 15382615613571648623UL) + ((uint64_t)op[6] * 1983579656587025256UL) + ((uint64_t)op[7] * 12755306781996716162UL) + ((uint64_t)op[8] * 15317309942447422085UL) + ((((uint64_t)op[9] * 7896789175154933191UL) + ((uint64_t)op[10] * 18077693477038659959UL)) * 8);
	tmp_q[9] = ((uint64_t)op[0] * 18077693477038659959UL) + ((uint64_t)op[1] * 9115319693145717106UL) + ((uint64_t)op[2] * 14877187536627791829UL) + ((uint64_t)op[3] * 3192455549729999471UL) + ((uint64_t)op[4] * 9424587271780295328UL) + ((uint64_t)op[5] * 9354509106931889212UL) + ((uint64_t)op[6] * 15382615613571648623UL) + ((uint64_t)op[7] * 1983579656587025256UL) + ((uint64_t)op[8] * 12755306781996716162UL) + ((uint64_t)op[9] * 15317309942447422085UL) + ((uint64_t)op[10] * 7834081180110810680UL);
	tmp_q[10] = ((uint64_t)op[0] * 7896789175154933191UL) + ((uint64_t)op[1] * 18077693477038659959UL) + ((uint64_t)op[2] * 9115319693145717106UL) + ((uint64_t)op[3] * 14877187536627791829UL) + ((uint64_t)op[4] * 3192455549729999471UL) + ((uint64_t)op[5] * 9424587271780295328UL) + ((uint64_t)op[6] * 9354509106931889212UL) + ((uint64_t)op[7] * 15382615613571648623UL) + ((uint64_t)op[8] * 1983579656587025256UL) + ((uint64_t)op[9] * 12755306781996716162UL) + ((uint64_t)op[10] * 15317309942447422085UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 17091129130419L) + ((((int128)tmp_q[1] * 126794139226581L) - ((int128)tmp_q[2] * 60715814863912L) - ((int128)tmp_q[3] * 12257740293774L) - ((int128)tmp_q[4] * 38433240718113L) - ((int128)tmp_q[5] * 53203279154758L) + ((int128)tmp_q[6] * 27847835381356L) - ((int128)tmp_q[7] * 71421553164248L) + ((int128)tmp_q[8] * 42790867621879L) + ((int128)tmp_q[9] * 83546602093356L) - ((int128)tmp_q[10] * 984226853830L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 984226853830L) + ((int128)tmp_q[1] * 17091129130419L) + ((((int128)tmp_q[2] * 126794139226581L) - ((int128)tmp_q[3] * 60715814863912L) - ((int128)tmp_q[4] * 12257740293774L) - ((int128)tmp_q[5] * 38433240718113L) - ((int128)tmp_q[6] * 53203279154758L) + ((int128)tmp_q[7] * 27847835381356L) - ((int128)tmp_q[8] * 71421553164248L) + ((int128)tmp_q[9] * 42790867621879L) + ((int128)tmp_q[10] * 83546602093356L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 83546602093356L) - ((int128)tmp_q[1] * 984226853830L) + ((int128)tmp_q[2] * 17091129130419L) + ((((int128)tmp_q[3] * 126794139226581L) - ((int128)tmp_q[4] * 60715814863912L) - ((int128)tmp_q[5] * 12257740293774L) - ((int128)tmp_q[6] * 38433240718113L) - ((int128)tmp_q[7] * 53203279154758L) + ((int128)tmp_q[8] * 27847835381356L) - ((int128)tmp_q[9] * 71421553164248L) + ((int128)tmp_q[10] * 42790867621879L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 42790867621879L) + ((int128)tmp_q[1] * 83546602093356L) - ((int128)tmp_q[2] * 984226853830L) + ((int128)tmp_q[3] * 17091129130419L) + ((((int128)tmp_q[4] * 126794139226581L) - ((int128)tmp_q[5] * 60715814863912L) - ((int128)tmp_q[6] * 12257740293774L) - ((int128)tmp_q[7] * 38433240718113L) - ((int128)tmp_q[8] * 53203279154758L) + ((int128)tmp_q[9] * 27847835381356L) - ((int128)tmp_q[10] * 71421553164248L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 71421553164248L) + ((int128)tmp_q[1] * 42790867621879L) + ((int128)tmp_q[2] * 83546602093356L) - ((int128)tmp_q[3] * 984226853830L) + ((int128)tmp_q[4] * 17091129130419L) + ((((int128)tmp_q[5] * 126794139226581L) - ((int128)tmp_q[6] * 60715814863912L) - ((int128)tmp_q[7] * 12257740293774L) - ((int128)tmp_q[8] * 38433240718113L) - ((int128)tmp_q[9] * 53203279154758L) + ((int128)tmp_q[10] * 27847835381356L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 27847835381356L) - ((int128)tmp_q[1] * 71421553164248L) + ((int128)tmp_q[2] * 42790867621879L) + ((int128)tmp_q[3] * 83546602093356L) - ((int128)tmp_q[4] * 984226853830L) + ((int128)tmp_q[5] * 17091129130419L) + ((((int128)tmp_q[6] * 126794139226581L) - ((int128)tmp_q[7] * 60715814863912L) - ((int128)tmp_q[8] * 12257740293774L) - ((int128)tmp_q[9] * 38433240718113L) - ((int128)tmp_q[10] * 53203279154758L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 53203279154758L) + ((int128)tmp_q[1] * 27847835381356L) - ((int128)tmp_q[2] * 71421553164248L) + ((int128)tmp_q[3] * 42790867621879L) + ((int128)tmp_q[4] * 83546602093356L) - ((int128)tmp_q[5] * 984226853830L) + ((int128)tmp_q[6] * 17091129130419L) + ((((int128)tmp_q[7] * 126794139226581L) - ((int128)tmp_q[8] * 60715814863912L) - ((int128)tmp_q[9] * 12257740293774L) - ((int128)tmp_q[10] * 38433240718113L)) * 8);
	tmp_zero[7] = -((int128)tmp_q[0] * 38433240718113L) - ((int128)tmp_q[1] * 53203279154758L) + ((int128)tmp_q[2] * 27847835381356L) - ((int128)tmp_q[3] * 71421553164248L) + ((int128)tmp_q[4] * 42790867621879L) + ((int128)tmp_q[5] * 83546602093356L) - ((int128)tmp_q[6] * 984226853830L) + ((int128)tmp_q[7] * 17091129130419L) + ((((int128)tmp_q[8] * 126794139226581L) - ((int128)tmp_q[9] * 60715814863912L) - ((int128)tmp_q[10] * 12257740293774L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 12257740293774L) - ((int128)tmp_q[1] * 38433240718113L) - ((int128)tmp_q[2] * 53203279154758L) + ((int128)tmp_q[3] * 27847835381356L) - ((int128)tmp_q[4] * 71421553164248L) + ((int128)tmp_q[5] * 42790867621879L) + ((int128)tmp_q[6] * 83546602093356L) - ((int128)tmp_q[7] * 984226853830L) + ((int128)tmp_q[8] * 17091129130419L) + ((((int128)tmp_q[9] * 126794139226581L) - ((int128)tmp_q[10] * 60715814863912L)) * 8);
	tmp_zero[9] = -((int128)tmp_q[0] * 60715814863912L) - ((int128)tmp_q[1] * 12257740293774L) - ((int128)tmp_q[2] * 38433240718113L) - ((int128)tmp_q[3] * 53203279154758L) + ((int128)tmp_q[4] * 27847835381356L) - ((int128)tmp_q[5] * 71421553164248L) + ((int128)tmp_q[6] * 42790867621879L) + ((int128)tmp_q[7] * 83546602093356L) - ((int128)tmp_q[8] * 984226853830L) + ((int128)tmp_q[9] * 17091129130419L) + ((int128)tmp_q[10] * 1014353113812648L);
	tmp_zero[10] = ((int128)tmp_q[0] * 126794139226581L) - ((int128)tmp_q[1] * 60715814863912L) - ((int128)tmp_q[2] * 12257740293774L) - ((int128)tmp_q[3] * 38433240718113L) - ((int128)tmp_q[4] * 53203279154758L) + ((int128)tmp_q[5] * 27847835381356L) - ((int128)tmp_q[6] * 71421553164248L) + ((int128)tmp_q[7] * 42790867621879L) + ((int128)tmp_q[8] * 83546602093356L) - ((int128)tmp_q[9] * 984226853830L) + ((int128)tmp_q[10] * 17091129130419L);

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

