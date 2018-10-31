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
	tmp_q[0] = ((uint64_t)op[0] * 14046666055595018931UL) + ((((uint64_t)op[1] * 13727134387113831382UL) + ((uint64_t)op[2] * 15760254885456009892UL) + ((uint64_t)op[3] * 7413246680987626094UL) + ((uint64_t)op[4] * 6412178840472283244UL) + ((uint64_t)op[5] * 10925349792597239858UL) + ((uint64_t)op[6] * 12882844154144869295UL) + ((uint64_t)op[7] * 10531053262461838584UL) + ((uint64_t)op[8] * 3688889386642748936UL) + ((uint64_t)op[9] * 7447675040221524173UL) + ((uint64_t)op[10] * 8508338676159000086UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 8508338676159000086UL) + ((uint64_t)op[1] * 14046666055595018931UL) + ((((uint64_t)op[2] * 13727134387113831382UL) + ((uint64_t)op[3] * 15760254885456009892UL) + ((uint64_t)op[4] * 7413246680987626094UL) + ((uint64_t)op[5] * 6412178840472283244UL) + ((uint64_t)op[6] * 10925349792597239858UL) + ((uint64_t)op[7] * 12882844154144869295UL) + ((uint64_t)op[8] * 10531053262461838584UL) + ((uint64_t)op[9] * 3688889386642748936UL) + ((uint64_t)op[10] * 7447675040221524173UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 7447675040221524173UL) + ((uint64_t)op[1] * 8508338676159000086UL) + ((uint64_t)op[2] * 14046666055595018931UL) + ((((uint64_t)op[3] * 13727134387113831382UL) + ((uint64_t)op[4] * 15760254885456009892UL) + ((uint64_t)op[5] * 7413246680987626094UL) + ((uint64_t)op[6] * 6412178840472283244UL) + ((uint64_t)op[7] * 10925349792597239858UL) + ((uint64_t)op[8] * 12882844154144869295UL) + ((uint64_t)op[9] * 10531053262461838584UL) + ((uint64_t)op[10] * 3688889386642748936UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 3688889386642748936UL) + ((uint64_t)op[1] * 7447675040221524173UL) + ((uint64_t)op[2] * 8508338676159000086UL) + ((uint64_t)op[3] * 14046666055595018931UL) + ((((uint64_t)op[4] * 13727134387113831382UL) + ((uint64_t)op[5] * 15760254885456009892UL) + ((uint64_t)op[6] * 7413246680987626094UL) + ((uint64_t)op[7] * 6412178840472283244UL) + ((uint64_t)op[8] * 10925349792597239858UL) + ((uint64_t)op[9] * 12882844154144869295UL) + ((uint64_t)op[10] * 10531053262461838584UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 10531053262461838584UL) + ((uint64_t)op[1] * 3688889386642748936UL) + ((uint64_t)op[2] * 7447675040221524173UL) + ((uint64_t)op[3] * 8508338676159000086UL) + ((uint64_t)op[4] * 14046666055595018931UL) + ((((uint64_t)op[5] * 13727134387113831382UL) + ((uint64_t)op[6] * 15760254885456009892UL) + ((uint64_t)op[7] * 7413246680987626094UL) + ((uint64_t)op[8] * 6412178840472283244UL) + ((uint64_t)op[9] * 10925349792597239858UL) + ((uint64_t)op[10] * 12882844154144869295UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 12882844154144869295UL) + ((uint64_t)op[1] * 10531053262461838584UL) + ((uint64_t)op[2] * 3688889386642748936UL) + ((uint64_t)op[3] * 7447675040221524173UL) + ((uint64_t)op[4] * 8508338676159000086UL) + ((uint64_t)op[5] * 14046666055595018931UL) + ((((uint64_t)op[6] * 13727134387113831382UL) + ((uint64_t)op[7] * 15760254885456009892UL) + ((uint64_t)op[8] * 7413246680987626094UL) + ((uint64_t)op[9] * 6412178840472283244UL) + ((uint64_t)op[10] * 10925349792597239858UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 10925349792597239858UL) + ((uint64_t)op[1] * 12882844154144869295UL) + ((uint64_t)op[2] * 10531053262461838584UL) + ((uint64_t)op[3] * 3688889386642748936UL) + ((uint64_t)op[4] * 7447675040221524173UL) + ((uint64_t)op[5] * 8508338676159000086UL) + ((uint64_t)op[6] * 14046666055595018931UL) + ((((uint64_t)op[7] * 13727134387113831382UL) + ((uint64_t)op[8] * 15760254885456009892UL) + ((uint64_t)op[9] * 7413246680987626094UL) + ((uint64_t)op[10] * 6412178840472283244UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 6412178840472283244UL) + ((uint64_t)op[1] * 10925349792597239858UL) + ((uint64_t)op[2] * 12882844154144869295UL) + ((uint64_t)op[3] * 10531053262461838584UL) + ((uint64_t)op[4] * 3688889386642748936UL) + ((uint64_t)op[5] * 7447675040221524173UL) + ((uint64_t)op[6] * 8508338676159000086UL) + ((uint64_t)op[7] * 14046666055595018931UL) + ((((uint64_t)op[8] * 13727134387113831382UL) + ((uint64_t)op[9] * 15760254885456009892UL) + ((uint64_t)op[10] * 7413246680987626094UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 7413246680987626094UL) + ((uint64_t)op[1] * 6412178840472283244UL) + ((uint64_t)op[2] * 10925349792597239858UL) + ((uint64_t)op[3] * 12882844154144869295UL) + ((uint64_t)op[4] * 10531053262461838584UL) + ((uint64_t)op[5] * 3688889386642748936UL) + ((uint64_t)op[6] * 7447675040221524173UL) + ((uint64_t)op[7] * 8508338676159000086UL) + ((uint64_t)op[8] * 14046666055595018931UL) + ((((uint64_t)op[9] * 13727134387113831382UL) + ((uint64_t)op[10] * 15760254885456009892UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 15760254885456009892UL) + ((uint64_t)op[1] * 7413246680987626094UL) + ((uint64_t)op[2] * 6412178840472283244UL) + ((uint64_t)op[3] * 10925349792597239858UL) + ((uint64_t)op[4] * 12882844154144869295UL) + ((uint64_t)op[5] * 10531053262461838584UL) + ((uint64_t)op[6] * 3688889386642748936UL) + ((uint64_t)op[7] * 7447675040221524173UL) + ((uint64_t)op[8] * 8508338676159000086UL) + ((uint64_t)op[9] * 14046666055595018931UL) + ((uint64_t)op[10] * 3856220341249061594UL);
	tmp_q[10] = ((uint64_t)op[0] * 13727134387113831382UL) + ((uint64_t)op[1] * 15760254885456009892UL) + ((uint64_t)op[2] * 7413246680987626094UL) + ((uint64_t)op[3] * 6412178840472283244UL) + ((uint64_t)op[4] * 10925349792597239858UL) + ((uint64_t)op[5] * 12882844154144869295UL) + ((uint64_t)op[6] * 10531053262461838584UL) + ((uint64_t)op[7] * 3688889386642748936UL) + ((uint64_t)op[8] * 7447675040221524173UL) + ((uint64_t)op[9] * 8508338676159000086UL) + ((uint64_t)op[10] * 14046666055595018931UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 21001225812547L) + ((((int128)tmp_q[1] * 44493669441358L) - ((int128)tmp_q[2] * 77003492284469L) + ((int128)tmp_q[3] * 69353763775211L) - ((int128)tmp_q[4] * 31680320407519L) + ((int128)tmp_q[5] * 35154864097235L) + ((int128)tmp_q[6] * 83127124860157L) + ((int128)tmp_q[7] * 29835390787654L) - ((int128)tmp_q[8] * 42941823588526L) + ((int128)tmp_q[9] * 50606594198398L) - ((int128)tmp_q[10] * 58979888709301L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 58979888709301L) - ((int128)tmp_q[1] * 21001225812547L) + ((((int128)tmp_q[2] * 44493669441358L) - ((int128)tmp_q[3] * 77003492284469L) + ((int128)tmp_q[4] * 69353763775211L) - ((int128)tmp_q[5] * 31680320407519L) + ((int128)tmp_q[6] * 35154864097235L) + ((int128)tmp_q[7] * 83127124860157L) + ((int128)tmp_q[8] * 29835390787654L) - ((int128)tmp_q[9] * 42941823588526L) + ((int128)tmp_q[10] * 50606594198398L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 50606594198398L) - ((int128)tmp_q[1] * 58979888709301L) - ((int128)tmp_q[2] * 21001225812547L) + ((((int128)tmp_q[3] * 44493669441358L) - ((int128)tmp_q[4] * 77003492284469L) + ((int128)tmp_q[5] * 69353763775211L) - ((int128)tmp_q[6] * 31680320407519L) + ((int128)tmp_q[7] * 35154864097235L) + ((int128)tmp_q[8] * 83127124860157L) + ((int128)tmp_q[9] * 29835390787654L) - ((int128)tmp_q[10] * 42941823588526L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 42941823588526L) + ((int128)tmp_q[1] * 50606594198398L) - ((int128)tmp_q[2] * 58979888709301L) - ((int128)tmp_q[3] * 21001225812547L) + ((((int128)tmp_q[4] * 44493669441358L) - ((int128)tmp_q[5] * 77003492284469L) + ((int128)tmp_q[6] * 69353763775211L) - ((int128)tmp_q[7] * 31680320407519L) + ((int128)tmp_q[8] * 35154864097235L) + ((int128)tmp_q[9] * 83127124860157L) + ((int128)tmp_q[10] * 29835390787654L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 29835390787654L) - ((int128)tmp_q[1] * 42941823588526L) + ((int128)tmp_q[2] * 50606594198398L) - ((int128)tmp_q[3] * 58979888709301L) - ((int128)tmp_q[4] * 21001225812547L) + ((((int128)tmp_q[5] * 44493669441358L) - ((int128)tmp_q[6] * 77003492284469L) + ((int128)tmp_q[7] * 69353763775211L) - ((int128)tmp_q[8] * 31680320407519L) + ((int128)tmp_q[9] * 35154864097235L) + ((int128)tmp_q[10] * 83127124860157L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 83127124860157L) + ((int128)tmp_q[1] * 29835390787654L) - ((int128)tmp_q[2] * 42941823588526L) + ((int128)tmp_q[3] * 50606594198398L) - ((int128)tmp_q[4] * 58979888709301L) - ((int128)tmp_q[5] * 21001225812547L) + ((((int128)tmp_q[6] * 44493669441358L) - ((int128)tmp_q[7] * 77003492284469L) + ((int128)tmp_q[8] * 69353763775211L) - ((int128)tmp_q[9] * 31680320407519L) + ((int128)tmp_q[10] * 35154864097235L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 35154864097235L) + ((int128)tmp_q[1] * 83127124860157L) + ((int128)tmp_q[2] * 29835390787654L) - ((int128)tmp_q[3] * 42941823588526L) + ((int128)tmp_q[4] * 50606594198398L) - ((int128)tmp_q[5] * 58979888709301L) - ((int128)tmp_q[6] * 21001225812547L) + ((((int128)tmp_q[7] * 44493669441358L) - ((int128)tmp_q[8] * 77003492284469L) + ((int128)tmp_q[9] * 69353763775211L) - ((int128)tmp_q[10] * 31680320407519L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 31680320407519L) + ((int128)tmp_q[1] * 35154864097235L) + ((int128)tmp_q[2] * 83127124860157L) + ((int128)tmp_q[3] * 29835390787654L) - ((int128)tmp_q[4] * 42941823588526L) + ((int128)tmp_q[5] * 50606594198398L) - ((int128)tmp_q[6] * 58979888709301L) - ((int128)tmp_q[7] * 21001225812547L) + ((((int128)tmp_q[8] * 44493669441358L) - ((int128)tmp_q[9] * 77003492284469L) + ((int128)tmp_q[10] * 69353763775211L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 69353763775211L) - ((int128)tmp_q[1] * 31680320407519L) + ((int128)tmp_q[2] * 35154864097235L) + ((int128)tmp_q[3] * 83127124860157L) + ((int128)tmp_q[4] * 29835390787654L) - ((int128)tmp_q[5] * 42941823588526L) + ((int128)tmp_q[6] * 50606594198398L) - ((int128)tmp_q[7] * 58979888709301L) - ((int128)tmp_q[8] * 21001225812547L) + ((((int128)tmp_q[9] * 44493669441358L) - ((int128)tmp_q[10] * 77003492284469L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 77003492284469L) + ((int128)tmp_q[1] * 69353763775211L) - ((int128)tmp_q[2] * 31680320407519L) + ((int128)tmp_q[3] * 35154864097235L) + ((int128)tmp_q[4] * 83127124860157L) + ((int128)tmp_q[5] * 29835390787654L) - ((int128)tmp_q[6] * 42941823588526L) + ((int128)tmp_q[7] * 50606594198398L) - ((int128)tmp_q[8] * 58979888709301L) - ((int128)tmp_q[9] * 21001225812547L) + ((int128)tmp_q[10] * 311455686089506L);
	tmp_zero[10] = ((int128)tmp_q[0] * 44493669441358L) - ((int128)tmp_q[1] * 77003492284469L) + ((int128)tmp_q[2] * 69353763775211L) - ((int128)tmp_q[3] * 31680320407519L) + ((int128)tmp_q[4] * 35154864097235L) + ((int128)tmp_q[5] * 83127124860157L) + ((int128)tmp_q[6] * 29835390787654L) - ((int128)tmp_q[7] * 42941823588526L) + ((int128)tmp_q[8] * 50606594198398L) - ((int128)tmp_q[9] * 58979888709301L) - ((int128)tmp_q[10] * 21001225812547L);

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

