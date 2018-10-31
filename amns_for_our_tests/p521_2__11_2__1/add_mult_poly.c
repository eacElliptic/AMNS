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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3924592545441062793UL) + ((((uint64_t)op[1] * 9294611761905679284UL) + ((uint64_t)op[2] * 8376748358183247107UL) + ((uint64_t)op[3] * 4973102979265842654UL) + ((uint64_t)op[4] * 3631503543367703537UL) + ((uint64_t)op[5] * 6717984201616099862UL) + ((uint64_t)op[6] * 15607321781765079862UL) + ((uint64_t)op[7] * 16389354892969465965UL) + ((uint64_t)op[8] * 7512270091149622274UL) + ((uint64_t)op[9] * 12666450409243599205UL) + ((uint64_t)op[10] * 3963014979436180229UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 3963014979436180229UL) + ((uint64_t)op[1] * 3924592545441062793UL) + ((((uint64_t)op[2] * 9294611761905679284UL) + ((uint64_t)op[3] * 8376748358183247107UL) + ((uint64_t)op[4] * 4973102979265842654UL) + ((uint64_t)op[5] * 3631503543367703537UL) + ((uint64_t)op[6] * 6717984201616099862UL) + ((uint64_t)op[7] * 15607321781765079862UL) + ((uint64_t)op[8] * 16389354892969465965UL) + ((uint64_t)op[9] * 7512270091149622274UL) + ((uint64_t)op[10] * 12666450409243599205UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 12666450409243599205UL) + ((uint64_t)op[1] * 3963014979436180229UL) + ((uint64_t)op[2] * 3924592545441062793UL) + ((((uint64_t)op[3] * 9294611761905679284UL) + ((uint64_t)op[4] * 8376748358183247107UL) + ((uint64_t)op[5] * 4973102979265842654UL) + ((uint64_t)op[6] * 3631503543367703537UL) + ((uint64_t)op[7] * 6717984201616099862UL) + ((uint64_t)op[8] * 15607321781765079862UL) + ((uint64_t)op[9] * 16389354892969465965UL) + ((uint64_t)op[10] * 7512270091149622274UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 7512270091149622274UL) + ((uint64_t)op[1] * 12666450409243599205UL) + ((uint64_t)op[2] * 3963014979436180229UL) + ((uint64_t)op[3] * 3924592545441062793UL) + ((((uint64_t)op[4] * 9294611761905679284UL) + ((uint64_t)op[5] * 8376748358183247107UL) + ((uint64_t)op[6] * 4973102979265842654UL) + ((uint64_t)op[7] * 3631503543367703537UL) + ((uint64_t)op[8] * 6717984201616099862UL) + ((uint64_t)op[9] * 15607321781765079862UL) + ((uint64_t)op[10] * 16389354892969465965UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 16389354892969465965UL) + ((uint64_t)op[1] * 7512270091149622274UL) + ((uint64_t)op[2] * 12666450409243599205UL) + ((uint64_t)op[3] * 3963014979436180229UL) + ((uint64_t)op[4] * 3924592545441062793UL) + ((((uint64_t)op[5] * 9294611761905679284UL) + ((uint64_t)op[6] * 8376748358183247107UL) + ((uint64_t)op[7] * 4973102979265842654UL) + ((uint64_t)op[8] * 3631503543367703537UL) + ((uint64_t)op[9] * 6717984201616099862UL) + ((uint64_t)op[10] * 15607321781765079862UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 15607321781765079862UL) + ((uint64_t)op[1] * 16389354892969465965UL) + ((uint64_t)op[2] * 7512270091149622274UL) + ((uint64_t)op[3] * 12666450409243599205UL) + ((uint64_t)op[4] * 3963014979436180229UL) + ((uint64_t)op[5] * 3924592545441062793UL) + ((((uint64_t)op[6] * 9294611761905679284UL) + ((uint64_t)op[7] * 8376748358183247107UL) + ((uint64_t)op[8] * 4973102979265842654UL) + ((uint64_t)op[9] * 3631503543367703537UL) + ((uint64_t)op[10] * 6717984201616099862UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 6717984201616099862UL) + ((uint64_t)op[1] * 15607321781765079862UL) + ((uint64_t)op[2] * 16389354892969465965UL) + ((uint64_t)op[3] * 7512270091149622274UL) + ((uint64_t)op[4] * 12666450409243599205UL) + ((uint64_t)op[5] * 3963014979436180229UL) + ((uint64_t)op[6] * 3924592545441062793UL) + ((((uint64_t)op[7] * 9294611761905679284UL) + ((uint64_t)op[8] * 8376748358183247107UL) + ((uint64_t)op[9] * 4973102979265842654UL) + ((uint64_t)op[10] * 3631503543367703537UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 3631503543367703537UL) + ((uint64_t)op[1] * 6717984201616099862UL) + ((uint64_t)op[2] * 15607321781765079862UL) + ((uint64_t)op[3] * 16389354892969465965UL) + ((uint64_t)op[4] * 7512270091149622274UL) + ((uint64_t)op[5] * 12666450409243599205UL) + ((uint64_t)op[6] * 3963014979436180229UL) + ((uint64_t)op[7] * 3924592545441062793UL) + ((((uint64_t)op[8] * 9294611761905679284UL) + ((uint64_t)op[9] * 8376748358183247107UL) + ((uint64_t)op[10] * 4973102979265842654UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 4973102979265842654UL) + ((uint64_t)op[1] * 3631503543367703537UL) + ((uint64_t)op[2] * 6717984201616099862UL) + ((uint64_t)op[3] * 15607321781765079862UL) + ((uint64_t)op[4] * 16389354892969465965UL) + ((uint64_t)op[5] * 7512270091149622274UL) + ((uint64_t)op[6] * 12666450409243599205UL) + ((uint64_t)op[7] * 3963014979436180229UL) + ((uint64_t)op[8] * 3924592545441062793UL) + ((((uint64_t)op[9] * 9294611761905679284UL) + ((uint64_t)op[10] * 8376748358183247107UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 8376748358183247107UL) + ((uint64_t)op[1] * 4973102979265842654UL) + ((uint64_t)op[2] * 3631503543367703537UL) + ((uint64_t)op[3] * 6717984201616099862UL) + ((uint64_t)op[4] * 15607321781765079862UL) + ((uint64_t)op[5] * 16389354892969465965UL) + ((uint64_t)op[6] * 7512270091149622274UL) + ((uint64_t)op[7] * 12666450409243599205UL) + ((uint64_t)op[8] * 3963014979436180229UL) + ((uint64_t)op[9] * 3924592545441062793UL) + ((uint64_t)op[10] * 142479450101806952UL);
	tmp_q[10] = ((uint64_t)op[0] * 9294611761905679284UL) + ((uint64_t)op[1] * 8376748358183247107UL) + ((uint64_t)op[2] * 4973102979265842654UL) + ((uint64_t)op[3] * 3631503543367703537UL) + ((uint64_t)op[4] * 6717984201616099862UL) + ((uint64_t)op[5] * 15607321781765079862UL) + ((uint64_t)op[6] * 16389354892969465965UL) + ((uint64_t)op[7] * 7512270091149622274UL) + ((uint64_t)op[8] * 12666450409243599205UL) + ((uint64_t)op[9] * 3963014979436180229UL) + ((uint64_t)op[10] * 3924592545441062793UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 22027651108917L) + ((-((int128)tmp_q[1] * 87975868049191L) - ((int128)tmp_q[2] * 6785694917982L) + ((int128)tmp_q[3] * 95417023519301L) + ((int128)tmp_q[4] * 4148608318308L) - ((int128)tmp_q[5] * 82996923405520L) + ((int128)tmp_q[6] * 34749022693828L) + ((int128)tmp_q[7] * 13233824289516L) + ((int128)tmp_q[8] * 27932696071653L) - ((int128)tmp_q[9] * 81591714821414L) + ((int128)tmp_q[10] * 1758724138133L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1758724138133L) + ((int128)tmp_q[1] * 22027651108917L) + ((-((int128)tmp_q[2] * 87975868049191L) - ((int128)tmp_q[3] * 6785694917982L) + ((int128)tmp_q[4] * 95417023519301L) + ((int128)tmp_q[5] * 4148608318308L) - ((int128)tmp_q[6] * 82996923405520L) + ((int128)tmp_q[7] * 34749022693828L) + ((int128)tmp_q[8] * 13233824289516L) + ((int128)tmp_q[9] * 27932696071653L) - ((int128)tmp_q[10] * 81591714821414L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 81591714821414L) + ((int128)tmp_q[1] * 1758724138133L) + ((int128)tmp_q[2] * 22027651108917L) + ((-((int128)tmp_q[3] * 87975868049191L) - ((int128)tmp_q[4] * 6785694917982L) + ((int128)tmp_q[5] * 95417023519301L) + ((int128)tmp_q[6] * 4148608318308L) - ((int128)tmp_q[7] * 82996923405520L) + ((int128)tmp_q[8] * 34749022693828L) + ((int128)tmp_q[9] * 13233824289516L) + ((int128)tmp_q[10] * 27932696071653L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 27932696071653L) - ((int128)tmp_q[1] * 81591714821414L) + ((int128)tmp_q[2] * 1758724138133L) + ((int128)tmp_q[3] * 22027651108917L) + ((-((int128)tmp_q[4] * 87975868049191L) - ((int128)tmp_q[5] * 6785694917982L) + ((int128)tmp_q[6] * 95417023519301L) + ((int128)tmp_q[7] * 4148608318308L) - ((int128)tmp_q[8] * 82996923405520L) + ((int128)tmp_q[9] * 34749022693828L) + ((int128)tmp_q[10] * 13233824289516L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 13233824289516L) + ((int128)tmp_q[1] * 27932696071653L) - ((int128)tmp_q[2] * 81591714821414L) + ((int128)tmp_q[3] * 1758724138133L) + ((int128)tmp_q[4] * 22027651108917L) + ((-((int128)tmp_q[5] * 87975868049191L) - ((int128)tmp_q[6] * 6785694917982L) + ((int128)tmp_q[7] * 95417023519301L) + ((int128)tmp_q[8] * 4148608318308L) - ((int128)tmp_q[9] * 82996923405520L) + ((int128)tmp_q[10] * 34749022693828L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 34749022693828L) + ((int128)tmp_q[1] * 13233824289516L) + ((int128)tmp_q[2] * 27932696071653L) - ((int128)tmp_q[3] * 81591714821414L) + ((int128)tmp_q[4] * 1758724138133L) + ((int128)tmp_q[5] * 22027651108917L) + ((-((int128)tmp_q[6] * 87975868049191L) - ((int128)tmp_q[7] * 6785694917982L) + ((int128)tmp_q[8] * 95417023519301L) + ((int128)tmp_q[9] * 4148608318308L) - ((int128)tmp_q[10] * 82996923405520L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 82996923405520L) + ((int128)tmp_q[1] * 34749022693828L) + ((int128)tmp_q[2] * 13233824289516L) + ((int128)tmp_q[3] * 27932696071653L) - ((int128)tmp_q[4] * 81591714821414L) + ((int128)tmp_q[5] * 1758724138133L) + ((int128)tmp_q[6] * 22027651108917L) + ((-((int128)tmp_q[7] * 87975868049191L) - ((int128)tmp_q[8] * 6785694917982L) + ((int128)tmp_q[9] * 95417023519301L) + ((int128)tmp_q[10] * 4148608318308L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 4148608318308L) - ((int128)tmp_q[1] * 82996923405520L) + ((int128)tmp_q[2] * 34749022693828L) + ((int128)tmp_q[3] * 13233824289516L) + ((int128)tmp_q[4] * 27932696071653L) - ((int128)tmp_q[5] * 81591714821414L) + ((int128)tmp_q[6] * 1758724138133L) + ((int128)tmp_q[7] * 22027651108917L) + ((-((int128)tmp_q[8] * 87975868049191L) - ((int128)tmp_q[9] * 6785694917982L) + ((int128)tmp_q[10] * 95417023519301L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 95417023519301L) + ((int128)tmp_q[1] * 4148608318308L) - ((int128)tmp_q[2] * 82996923405520L) + ((int128)tmp_q[3] * 34749022693828L) + ((int128)tmp_q[4] * 13233824289516L) + ((int128)tmp_q[5] * 27932696071653L) - ((int128)tmp_q[6] * 81591714821414L) + ((int128)tmp_q[7] * 1758724138133L) + ((int128)tmp_q[8] * 22027651108917L) + ((-((int128)tmp_q[9] * 87975868049191L) - ((int128)tmp_q[10] * 6785694917982L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 6785694917982L) + ((int128)tmp_q[1] * 95417023519301L) + ((int128)tmp_q[2] * 4148608318308L) - ((int128)tmp_q[3] * 82996923405520L) + ((int128)tmp_q[4] * 34749022693828L) + ((int128)tmp_q[5] * 13233824289516L) + ((int128)tmp_q[6] * 27932696071653L) - ((int128)tmp_q[7] * 81591714821414L) + ((int128)tmp_q[8] * 1758724138133L) + ((int128)tmp_q[9] * 22027651108917L) - ((int128)tmp_q[10] * 175951736098382L);
	tmp_zero[10] = -((int128)tmp_q[0] * 87975868049191L) - ((int128)tmp_q[1] * 6785694917982L) + ((int128)tmp_q[2] * 95417023519301L) + ((int128)tmp_q[3] * 4148608318308L) - ((int128)tmp_q[4] * 82996923405520L) + ((int128)tmp_q[5] * 34749022693828L) + ((int128)tmp_q[6] * 13233824289516L) + ((int128)tmp_q[7] * 27932696071653L) - ((int128)tmp_q[8] * 81591714821414L) + ((int128)tmp_q[9] * 1758724138133L) + ((int128)tmp_q[10] * 22027651108917L);

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

