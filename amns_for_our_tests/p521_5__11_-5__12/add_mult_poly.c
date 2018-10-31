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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3467738515399735100UL) + ((((uint64_t)op[1] * 5481740930395249026UL) + ((uint64_t)op[2] * 4414832429634507167UL) + ((uint64_t)op[3] * 6297156649723115876UL) + ((uint64_t)op[4] * 3267188526617048291UL) + ((uint64_t)op[5] * 14142017263886853389UL) + ((uint64_t)op[6] * 10065014271304783847UL) + ((uint64_t)op[7] * 243034440342256713UL) + ((uint64_t)op[8] * 16782230240695485653UL) + ((uint64_t)op[9] * 8138829103953240050UL) + ((uint64_t)op[10] * 5599251018076713UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 5599251018076713UL) + ((uint64_t)op[1] * 3467738515399735100UL) + ((((uint64_t)op[2] * 5481740930395249026UL) + ((uint64_t)op[3] * 4414832429634507167UL) + ((uint64_t)op[4] * 6297156649723115876UL) + ((uint64_t)op[5] * 3267188526617048291UL) + ((uint64_t)op[6] * 14142017263886853389UL) + ((uint64_t)op[7] * 10065014271304783847UL) + ((uint64_t)op[8] * 243034440342256713UL) + ((uint64_t)op[9] * 16782230240695485653UL) + ((uint64_t)op[10] * 8138829103953240050UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 8138829103953240050UL) + ((uint64_t)op[1] * 5599251018076713UL) + ((uint64_t)op[2] * 3467738515399735100UL) + ((((uint64_t)op[3] * 5481740930395249026UL) + ((uint64_t)op[4] * 4414832429634507167UL) + ((uint64_t)op[5] * 6297156649723115876UL) + ((uint64_t)op[6] * 3267188526617048291UL) + ((uint64_t)op[7] * 14142017263886853389UL) + ((uint64_t)op[8] * 10065014271304783847UL) + ((uint64_t)op[9] * 243034440342256713UL) + ((uint64_t)op[10] * 16782230240695485653UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16782230240695485653UL) + ((uint64_t)op[1] * 8138829103953240050UL) + ((uint64_t)op[2] * 5599251018076713UL) + ((uint64_t)op[3] * 3467738515399735100UL) + ((((uint64_t)op[4] * 5481740930395249026UL) + ((uint64_t)op[5] * 4414832429634507167UL) + ((uint64_t)op[6] * 6297156649723115876UL) + ((uint64_t)op[7] * 3267188526617048291UL) + ((uint64_t)op[8] * 14142017263886853389UL) + ((uint64_t)op[9] * 10065014271304783847UL) + ((uint64_t)op[10] * 243034440342256713UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 243034440342256713UL) + ((uint64_t)op[1] * 16782230240695485653UL) + ((uint64_t)op[2] * 8138829103953240050UL) + ((uint64_t)op[3] * 5599251018076713UL) + ((uint64_t)op[4] * 3467738515399735100UL) + ((((uint64_t)op[5] * 5481740930395249026UL) + ((uint64_t)op[6] * 4414832429634507167UL) + ((uint64_t)op[7] * 6297156649723115876UL) + ((uint64_t)op[8] * 3267188526617048291UL) + ((uint64_t)op[9] * 14142017263886853389UL) + ((uint64_t)op[10] * 10065014271304783847UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 10065014271304783847UL) + ((uint64_t)op[1] * 243034440342256713UL) + ((uint64_t)op[2] * 16782230240695485653UL) + ((uint64_t)op[3] * 8138829103953240050UL) + ((uint64_t)op[4] * 5599251018076713UL) + ((uint64_t)op[5] * 3467738515399735100UL) + ((((uint64_t)op[6] * 5481740930395249026UL) + ((uint64_t)op[7] * 4414832429634507167UL) + ((uint64_t)op[8] * 6297156649723115876UL) + ((uint64_t)op[9] * 3267188526617048291UL) + ((uint64_t)op[10] * 14142017263886853389UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 14142017263886853389UL) + ((uint64_t)op[1] * 10065014271304783847UL) + ((uint64_t)op[2] * 243034440342256713UL) + ((uint64_t)op[3] * 16782230240695485653UL) + ((uint64_t)op[4] * 8138829103953240050UL) + ((uint64_t)op[5] * 5599251018076713UL) + ((uint64_t)op[6] * 3467738515399735100UL) + ((((uint64_t)op[7] * 5481740930395249026UL) + ((uint64_t)op[8] * 4414832429634507167UL) + ((uint64_t)op[9] * 6297156649723115876UL) + ((uint64_t)op[10] * 3267188526617048291UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 3267188526617048291UL) + ((uint64_t)op[1] * 14142017263886853389UL) + ((uint64_t)op[2] * 10065014271304783847UL) + ((uint64_t)op[3] * 243034440342256713UL) + ((uint64_t)op[4] * 16782230240695485653UL) + ((uint64_t)op[5] * 8138829103953240050UL) + ((uint64_t)op[6] * 5599251018076713UL) + ((uint64_t)op[7] * 3467738515399735100UL) + ((((uint64_t)op[8] * 5481740930395249026UL) + ((uint64_t)op[9] * 4414832429634507167UL) + ((uint64_t)op[10] * 6297156649723115876UL)) * 18446744073709551611);
	tmp_q[8] = ((uint64_t)op[0] * 6297156649723115876UL) + ((uint64_t)op[1] * 3267188526617048291UL) + ((uint64_t)op[2] * 14142017263886853389UL) + ((uint64_t)op[3] * 10065014271304783847UL) + ((uint64_t)op[4] * 243034440342256713UL) + ((uint64_t)op[5] * 16782230240695485653UL) + ((uint64_t)op[6] * 8138829103953240050UL) + ((uint64_t)op[7] * 5599251018076713UL) + ((uint64_t)op[8] * 3467738515399735100UL) + ((((uint64_t)op[9] * 5481740930395249026UL) + ((uint64_t)op[10] * 4414832429634507167UL)) * 18446744073709551611);
	tmp_q[9] = ((uint64_t)op[0] * 4414832429634507167UL) + ((uint64_t)op[1] * 6297156649723115876UL) + ((uint64_t)op[2] * 3267188526617048291UL) + ((uint64_t)op[3] * 14142017263886853389UL) + ((uint64_t)op[4] * 10065014271304783847UL) + ((uint64_t)op[5] * 243034440342256713UL) + ((uint64_t)op[6] * 16782230240695485653UL) + ((uint64_t)op[7] * 8138829103953240050UL) + ((uint64_t)op[8] * 5599251018076713UL) + ((uint64_t)op[9] * 3467738515399735100UL) + ((uint64_t)op[10] * 9484783495442858102UL);
	tmp_q[10] = ((uint64_t)op[0] * 5481740930395249026UL) + ((uint64_t)op[1] * 4414832429634507167UL) + ((uint64_t)op[2] * 6297156649723115876UL) + ((uint64_t)op[3] * 3267188526617048291UL) + ((uint64_t)op[4] * 14142017263886853389UL) + ((uint64_t)op[5] * 10065014271304783847UL) + ((uint64_t)op[6] * 243034440342256713UL) + ((uint64_t)op[7] * 16782230240695485653UL) + ((uint64_t)op[8] * 8138829103953240050UL) + ((uint64_t)op[9] * 5599251018076713UL) + ((uint64_t)op[10] * 3467738515399735100UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 110742357262282L) - ((((int128)tmp_q[1] * 306486274754L) - ((int128)tmp_q[2] * 4971964600529L) + ((int128)tmp_q[3] * 27254721890921L) - ((int128)tmp_q[4] * 71739492194924L) - ((int128)tmp_q[5] * 51895018647351L) + ((int128)tmp_q[6] * 57893300744108L) + ((int128)tmp_q[7] * 73481362701343L) - ((int128)tmp_q[8] * 54232925047373L) + ((int128)tmp_q[9] * 91080668153764L) - ((int128)tmp_q[10] * 69128026284814L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 69128026284814L) - ((int128)tmp_q[1] * 110742357262282L) - ((((int128)tmp_q[2] * 306486274754L) - ((int128)tmp_q[3] * 4971964600529L) + ((int128)tmp_q[4] * 27254721890921L) - ((int128)tmp_q[5] * 71739492194924L) - ((int128)tmp_q[6] * 51895018647351L) + ((int128)tmp_q[7] * 57893300744108L) + ((int128)tmp_q[8] * 73481362701343L) - ((int128)tmp_q[9] * 54232925047373L) + ((int128)tmp_q[10] * 91080668153764L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 91080668153764L) - ((int128)tmp_q[1] * 69128026284814L) - ((int128)tmp_q[2] * 110742357262282L) - ((((int128)tmp_q[3] * 306486274754L) - ((int128)tmp_q[4] * 4971964600529L) + ((int128)tmp_q[5] * 27254721890921L) - ((int128)tmp_q[6] * 71739492194924L) - ((int128)tmp_q[7] * 51895018647351L) + ((int128)tmp_q[8] * 57893300744108L) + ((int128)tmp_q[9] * 73481362701343L) - ((int128)tmp_q[10] * 54232925047373L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 54232925047373L) + ((int128)tmp_q[1] * 91080668153764L) - ((int128)tmp_q[2] * 69128026284814L) - ((int128)tmp_q[3] * 110742357262282L) - ((((int128)tmp_q[4] * 306486274754L) - ((int128)tmp_q[5] * 4971964600529L) + ((int128)tmp_q[6] * 27254721890921L) - ((int128)tmp_q[7] * 71739492194924L) - ((int128)tmp_q[8] * 51895018647351L) + ((int128)tmp_q[9] * 57893300744108L) + ((int128)tmp_q[10] * 73481362701343L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 73481362701343L) - ((int128)tmp_q[1] * 54232925047373L) + ((int128)tmp_q[2] * 91080668153764L) - ((int128)tmp_q[3] * 69128026284814L) - ((int128)tmp_q[4] * 110742357262282L) - ((((int128)tmp_q[5] * 306486274754L) - ((int128)tmp_q[6] * 4971964600529L) + ((int128)tmp_q[7] * 27254721890921L) - ((int128)tmp_q[8] * 71739492194924L) - ((int128)tmp_q[9] * 51895018647351L) + ((int128)tmp_q[10] * 57893300744108L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 57893300744108L) + ((int128)tmp_q[1] * 73481362701343L) - ((int128)tmp_q[2] * 54232925047373L) + ((int128)tmp_q[3] * 91080668153764L) - ((int128)tmp_q[4] * 69128026284814L) - ((int128)tmp_q[5] * 110742357262282L) - ((((int128)tmp_q[6] * 306486274754L) - ((int128)tmp_q[7] * 4971964600529L) + ((int128)tmp_q[8] * 27254721890921L) - ((int128)tmp_q[9] * 71739492194924L) - ((int128)tmp_q[10] * 51895018647351L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 51895018647351L) + ((int128)tmp_q[1] * 57893300744108L) + ((int128)tmp_q[2] * 73481362701343L) - ((int128)tmp_q[3] * 54232925047373L) + ((int128)tmp_q[4] * 91080668153764L) - ((int128)tmp_q[5] * 69128026284814L) - ((int128)tmp_q[6] * 110742357262282L) - ((((int128)tmp_q[7] * 306486274754L) - ((int128)tmp_q[8] * 4971964600529L) + ((int128)tmp_q[9] * 27254721890921L) - ((int128)tmp_q[10] * 71739492194924L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 71739492194924L) - ((int128)tmp_q[1] * 51895018647351L) + ((int128)tmp_q[2] * 57893300744108L) + ((int128)tmp_q[3] * 73481362701343L) - ((int128)tmp_q[4] * 54232925047373L) + ((int128)tmp_q[5] * 91080668153764L) - ((int128)tmp_q[6] * 69128026284814L) - ((int128)tmp_q[7] * 110742357262282L) - ((((int128)tmp_q[8] * 306486274754L) - ((int128)tmp_q[9] * 4971964600529L) + ((int128)tmp_q[10] * 27254721890921L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 27254721890921L) - ((int128)tmp_q[1] * 71739492194924L) - ((int128)tmp_q[2] * 51895018647351L) + ((int128)tmp_q[3] * 57893300744108L) + ((int128)tmp_q[4] * 73481362701343L) - ((int128)tmp_q[5] * 54232925047373L) + ((int128)tmp_q[6] * 91080668153764L) - ((int128)tmp_q[7] * 69128026284814L) - ((int128)tmp_q[8] * 110742357262282L) - ((((int128)tmp_q[9] * 306486274754L) - ((int128)tmp_q[10] * 4971964600529L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 4971964600529L) + ((int128)tmp_q[1] * 27254721890921L) - ((int128)tmp_q[2] * 71739492194924L) - ((int128)tmp_q[3] * 51895018647351L) + ((int128)tmp_q[4] * 57893300744108L) + ((int128)tmp_q[5] * 73481362701343L) - ((int128)tmp_q[6] * 54232925047373L) + ((int128)tmp_q[7] * 91080668153764L) - ((int128)tmp_q[8] * 69128026284814L) - ((int128)tmp_q[9] * 110742357262282L) - ((int128)tmp_q[10] * 1532431373770L);
	tmp_zero[10] = ((int128)tmp_q[0] * 306486274754L) - ((int128)tmp_q[1] * 4971964600529L) + ((int128)tmp_q[2] * 27254721890921L) - ((int128)tmp_q[3] * 71739492194924L) - ((int128)tmp_q[4] * 51895018647351L) + ((int128)tmp_q[5] * 57893300744108L) + ((int128)tmp_q[6] * 73481362701343L) - ((int128)tmp_q[7] * 54232925047373L) + ((int128)tmp_q[8] * 91080668153764L) - ((int128)tmp_q[9] * 69128026284814L) - ((int128)tmp_q[10] * 110742357262282L);

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

