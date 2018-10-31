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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 5);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 10);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 10);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 5);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11807281361146987837UL) + ((((uint64_t)op[1] * 15769946635396677728UL) + ((uint64_t)op[2] * 11806021937643102228UL) + ((uint64_t)op[3] * 12663329286194755152UL) + ((uint64_t)op[4] * 440206687894763904UL) + ((uint64_t)op[5] * 16745027192440930508UL) + ((uint64_t)op[6] * 2731250561113367518UL) + ((uint64_t)op[7] * 4460059072069131818UL) + ((uint64_t)op[8] * 16406034343609594324UL) + ((uint64_t)op[9] * 10131181015433319287UL) + ((uint64_t)op[10] * 15256327758115961151UL) + ((uint64_t)op[11] * 3719345613807219042UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 3719345613807219042UL) + ((uint64_t)op[1] * 11807281361146987837UL) + ((((uint64_t)op[2] * 15769946635396677728UL) + ((uint64_t)op[3] * 11806021937643102228UL) + ((uint64_t)op[4] * 12663329286194755152UL) + ((uint64_t)op[5] * 440206687894763904UL) + ((uint64_t)op[6] * 16745027192440930508UL) + ((uint64_t)op[7] * 2731250561113367518UL) + ((uint64_t)op[8] * 4460059072069131818UL) + ((uint64_t)op[9] * 16406034343609594324UL) + ((uint64_t)op[10] * 10131181015433319287UL) + ((uint64_t)op[11] * 15256327758115961151UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 15256327758115961151UL) + ((uint64_t)op[1] * 3719345613807219042UL) + ((uint64_t)op[2] * 11807281361146987837UL) + ((((uint64_t)op[3] * 15769946635396677728UL) + ((uint64_t)op[4] * 11806021937643102228UL) + ((uint64_t)op[5] * 12663329286194755152UL) + ((uint64_t)op[6] * 440206687894763904UL) + ((uint64_t)op[7] * 16745027192440930508UL) + ((uint64_t)op[8] * 2731250561113367518UL) + ((uint64_t)op[9] * 4460059072069131818UL) + ((uint64_t)op[10] * 16406034343609594324UL) + ((uint64_t)op[11] * 10131181015433319287UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 10131181015433319287UL) + ((uint64_t)op[1] * 15256327758115961151UL) + ((uint64_t)op[2] * 3719345613807219042UL) + ((uint64_t)op[3] * 11807281361146987837UL) + ((((uint64_t)op[4] * 15769946635396677728UL) + ((uint64_t)op[5] * 11806021937643102228UL) + ((uint64_t)op[6] * 12663329286194755152UL) + ((uint64_t)op[7] * 440206687894763904UL) + ((uint64_t)op[8] * 16745027192440930508UL) + ((uint64_t)op[9] * 2731250561113367518UL) + ((uint64_t)op[10] * 4460059072069131818UL) + ((uint64_t)op[11] * 16406034343609594324UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 16406034343609594324UL) + ((uint64_t)op[1] * 10131181015433319287UL) + ((uint64_t)op[2] * 15256327758115961151UL) + ((uint64_t)op[3] * 3719345613807219042UL) + ((uint64_t)op[4] * 11807281361146987837UL) + ((((uint64_t)op[5] * 15769946635396677728UL) + ((uint64_t)op[6] * 11806021937643102228UL) + ((uint64_t)op[7] * 12663329286194755152UL) + ((uint64_t)op[8] * 440206687894763904UL) + ((uint64_t)op[9] * 16745027192440930508UL) + ((uint64_t)op[10] * 2731250561113367518UL) + ((uint64_t)op[11] * 4460059072069131818UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 4460059072069131818UL) + ((uint64_t)op[1] * 16406034343609594324UL) + ((uint64_t)op[2] * 10131181015433319287UL) + ((uint64_t)op[3] * 15256327758115961151UL) + ((uint64_t)op[4] * 3719345613807219042UL) + ((uint64_t)op[5] * 11807281361146987837UL) + ((((uint64_t)op[6] * 15769946635396677728UL) + ((uint64_t)op[7] * 11806021937643102228UL) + ((uint64_t)op[8] * 12663329286194755152UL) + ((uint64_t)op[9] * 440206687894763904UL) + ((uint64_t)op[10] * 16745027192440930508UL) + ((uint64_t)op[11] * 2731250561113367518UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 2731250561113367518UL) + ((uint64_t)op[1] * 4460059072069131818UL) + ((uint64_t)op[2] * 16406034343609594324UL) + ((uint64_t)op[3] * 10131181015433319287UL) + ((uint64_t)op[4] * 15256327758115961151UL) + ((uint64_t)op[5] * 3719345613807219042UL) + ((uint64_t)op[6] * 11807281361146987837UL) + ((((uint64_t)op[7] * 15769946635396677728UL) + ((uint64_t)op[8] * 11806021937643102228UL) + ((uint64_t)op[9] * 12663329286194755152UL) + ((uint64_t)op[10] * 440206687894763904UL) + ((uint64_t)op[11] * 16745027192440930508UL)) * 5);
	tmp_q[7] = ((uint64_t)op[0] * 16745027192440930508UL) + ((uint64_t)op[1] * 2731250561113367518UL) + ((uint64_t)op[2] * 4460059072069131818UL) + ((uint64_t)op[3] * 16406034343609594324UL) + ((uint64_t)op[4] * 10131181015433319287UL) + ((uint64_t)op[5] * 15256327758115961151UL) + ((uint64_t)op[6] * 3719345613807219042UL) + ((uint64_t)op[7] * 11807281361146987837UL) + ((((uint64_t)op[8] * 15769946635396677728UL) + ((uint64_t)op[9] * 11806021937643102228UL) + ((uint64_t)op[10] * 12663329286194755152UL) + ((uint64_t)op[11] * 440206687894763904UL)) * 5);
	tmp_q[8] = ((uint64_t)op[0] * 440206687894763904UL) + ((uint64_t)op[1] * 16745027192440930508UL) + ((uint64_t)op[2] * 2731250561113367518UL) + ((uint64_t)op[3] * 4460059072069131818UL) + ((uint64_t)op[4] * 16406034343609594324UL) + ((uint64_t)op[5] * 10131181015433319287UL) + ((uint64_t)op[6] * 15256327758115961151UL) + ((uint64_t)op[7] * 3719345613807219042UL) + ((uint64_t)op[8] * 11807281361146987837UL) + ((((uint64_t)op[9] * 15769946635396677728UL) + ((uint64_t)op[10] * 11806021937643102228UL) + ((uint64_t)op[11] * 12663329286194755152UL)) * 5);
	tmp_q[9] = ((uint64_t)op[0] * 12663329286194755152UL) + ((uint64_t)op[1] * 440206687894763904UL) + ((uint64_t)op[2] * 16745027192440930508UL) + ((uint64_t)op[3] * 2731250561113367518UL) + ((uint64_t)op[4] * 4460059072069131818UL) + ((uint64_t)op[5] * 16406034343609594324UL) + ((uint64_t)op[6] * 10131181015433319287UL) + ((uint64_t)op[7] * 15256327758115961151UL) + ((uint64_t)op[8] * 3719345613807219042UL) + ((uint64_t)op[9] * 11807281361146987837UL) + ((((uint64_t)op[10] * 15769946635396677728UL) + ((uint64_t)op[11] * 11806021937643102228UL)) * 5);
	tmp_q[10] = ((uint64_t)op[0] * 11806021937643102228UL) + ((uint64_t)op[1] * 12663329286194755152UL) + ((uint64_t)op[2] * 440206687894763904UL) + ((uint64_t)op[3] * 16745027192440930508UL) + ((uint64_t)op[4] * 2731250561113367518UL) + ((uint64_t)op[5] * 4460059072069131818UL) + ((uint64_t)op[6] * 16406034343609594324UL) + ((uint64_t)op[7] * 10131181015433319287UL) + ((uint64_t)op[8] * 15256327758115961151UL) + ((uint64_t)op[9] * 3719345613807219042UL) + ((uint64_t)op[10] * 11807281361146987837UL) + ((uint64_t)op[11] * 5062756882145182176UL);
	tmp_q[11] = ((uint64_t)op[0] * 15769946635396677728UL) + ((uint64_t)op[1] * 11806021937643102228UL) + ((uint64_t)op[2] * 12663329286194755152UL) + ((uint64_t)op[3] * 440206687894763904UL) + ((uint64_t)op[4] * 16745027192440930508UL) + ((uint64_t)op[5] * 2731250561113367518UL) + ((uint64_t)op[6] * 4460059072069131818UL) + ((uint64_t)op[7] * 16406034343609594324UL) + ((uint64_t)op[8] * 10131181015433319287UL) + ((uint64_t)op[9] * 15256327758115961151UL) + ((uint64_t)op[10] * 3719345613807219042UL) + ((uint64_t)op[11] * 11807281361146987837UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 644053800861L) + ((((int128)tmp_q[1] * 5732031081337L) + ((int128)tmp_q[2] * 1229492444530L) - ((int128)tmp_q[3] * 4169766352232L) - ((int128)tmp_q[4] * 891133241793L) - ((int128)tmp_q[5] * 1285180556663L) + ((int128)tmp_q[6] * 3301406012735L) + ((int128)tmp_q[7] * 1992425239574L) + ((int128)tmp_q[8] * 3410169732919L) - ((int128)tmp_q[9] * 500176616740L) + ((int128)tmp_q[10] * 7734654581252L) + ((int128)tmp_q[11] * 757732270717L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 757732270717L) - ((int128)tmp_q[1] * 644053800861L) + ((((int128)tmp_q[2] * 5732031081337L) + ((int128)tmp_q[3] * 1229492444530L) - ((int128)tmp_q[4] * 4169766352232L) - ((int128)tmp_q[5] * 891133241793L) - ((int128)tmp_q[6] * 1285180556663L) + ((int128)tmp_q[7] * 3301406012735L) + ((int128)tmp_q[8] * 1992425239574L) + ((int128)tmp_q[9] * 3410169732919L) - ((int128)tmp_q[10] * 500176616740L) + ((int128)tmp_q[11] * 7734654581252L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 7734654581252L) + ((int128)tmp_q[1] * 757732270717L) - ((int128)tmp_q[2] * 644053800861L) + ((((int128)tmp_q[3] * 5732031081337L) + ((int128)tmp_q[4] * 1229492444530L) - ((int128)tmp_q[5] * 4169766352232L) - ((int128)tmp_q[6] * 891133241793L) - ((int128)tmp_q[7] * 1285180556663L) + ((int128)tmp_q[8] * 3301406012735L) + ((int128)tmp_q[9] * 1992425239574L) + ((int128)tmp_q[10] * 3410169732919L) - ((int128)tmp_q[11] * 500176616740L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 500176616740L) + ((int128)tmp_q[1] * 7734654581252L) + ((int128)tmp_q[2] * 757732270717L) - ((int128)tmp_q[3] * 644053800861L) + ((((int128)tmp_q[4] * 5732031081337L) + ((int128)tmp_q[5] * 1229492444530L) - ((int128)tmp_q[6] * 4169766352232L) - ((int128)tmp_q[7] * 891133241793L) - ((int128)tmp_q[8] * 1285180556663L) + ((int128)tmp_q[9] * 3301406012735L) + ((int128)tmp_q[10] * 1992425239574L) + ((int128)tmp_q[11] * 3410169732919L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 3410169732919L) - ((int128)tmp_q[1] * 500176616740L) + ((int128)tmp_q[2] * 7734654581252L) + ((int128)tmp_q[3] * 757732270717L) - ((int128)tmp_q[4] * 644053800861L) + ((((int128)tmp_q[5] * 5732031081337L) + ((int128)tmp_q[6] * 1229492444530L) - ((int128)tmp_q[7] * 4169766352232L) - ((int128)tmp_q[8] * 891133241793L) - ((int128)tmp_q[9] * 1285180556663L) + ((int128)tmp_q[10] * 3301406012735L) + ((int128)tmp_q[11] * 1992425239574L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 1992425239574L) + ((int128)tmp_q[1] * 3410169732919L) - ((int128)tmp_q[2] * 500176616740L) + ((int128)tmp_q[3] * 7734654581252L) + ((int128)tmp_q[4] * 757732270717L) - ((int128)tmp_q[5] * 644053800861L) + ((((int128)tmp_q[6] * 5732031081337L) + ((int128)tmp_q[7] * 1229492444530L) - ((int128)tmp_q[8] * 4169766352232L) - ((int128)tmp_q[9] * 891133241793L) - ((int128)tmp_q[10] * 1285180556663L) + ((int128)tmp_q[11] * 3301406012735L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 3301406012735L) + ((int128)tmp_q[1] * 1992425239574L) + ((int128)tmp_q[2] * 3410169732919L) - ((int128)tmp_q[3] * 500176616740L) + ((int128)tmp_q[4] * 7734654581252L) + ((int128)tmp_q[5] * 757732270717L) - ((int128)tmp_q[6] * 644053800861L) + ((((int128)tmp_q[7] * 5732031081337L) + ((int128)tmp_q[8] * 1229492444530L) - ((int128)tmp_q[9] * 4169766352232L) - ((int128)tmp_q[10] * 891133241793L) - ((int128)tmp_q[11] * 1285180556663L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 1285180556663L) + ((int128)tmp_q[1] * 3301406012735L) + ((int128)tmp_q[2] * 1992425239574L) + ((int128)tmp_q[3] * 3410169732919L) - ((int128)tmp_q[4] * 500176616740L) + ((int128)tmp_q[5] * 7734654581252L) + ((int128)tmp_q[6] * 757732270717L) - ((int128)tmp_q[7] * 644053800861L) + ((((int128)tmp_q[8] * 5732031081337L) + ((int128)tmp_q[9] * 1229492444530L) - ((int128)tmp_q[10] * 4169766352232L) - ((int128)tmp_q[11] * 891133241793L)) * 5);
	tmp_zero[8] = -((int128)tmp_q[0] * 891133241793L) - ((int128)tmp_q[1] * 1285180556663L) + ((int128)tmp_q[2] * 3301406012735L) + ((int128)tmp_q[3] * 1992425239574L) + ((int128)tmp_q[4] * 3410169732919L) - ((int128)tmp_q[5] * 500176616740L) + ((int128)tmp_q[6] * 7734654581252L) + ((int128)tmp_q[7] * 757732270717L) - ((int128)tmp_q[8] * 644053800861L) + ((((int128)tmp_q[9] * 5732031081337L) + ((int128)tmp_q[10] * 1229492444530L) - ((int128)tmp_q[11] * 4169766352232L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 4169766352232L) - ((int128)tmp_q[1] * 891133241793L) - ((int128)tmp_q[2] * 1285180556663L) + ((int128)tmp_q[3] * 3301406012735L) + ((int128)tmp_q[4] * 1992425239574L) + ((int128)tmp_q[5] * 3410169732919L) - ((int128)tmp_q[6] * 500176616740L) + ((int128)tmp_q[7] * 7734654581252L) + ((int128)tmp_q[8] * 757732270717L) - ((int128)tmp_q[9] * 644053800861L) + ((((int128)tmp_q[10] * 5732031081337L) + ((int128)tmp_q[11] * 1229492444530L)) * 5);
	tmp_zero[10] = ((int128)tmp_q[0] * 1229492444530L) - ((int128)tmp_q[1] * 4169766352232L) - ((int128)tmp_q[2] * 891133241793L) - ((int128)tmp_q[3] * 1285180556663L) + ((int128)tmp_q[4] * 3301406012735L) + ((int128)tmp_q[5] * 1992425239574L) + ((int128)tmp_q[6] * 3410169732919L) - ((int128)tmp_q[7] * 500176616740L) + ((int128)tmp_q[8] * 7734654581252L) + ((int128)tmp_q[9] * 757732270717L) - ((int128)tmp_q[10] * 644053800861L) + ((int128)tmp_q[11] * 28660155406685L);
	tmp_zero[11] = ((int128)tmp_q[0] * 5732031081337L) + ((int128)tmp_q[1] * 1229492444530L) - ((int128)tmp_q[2] * 4169766352232L) - ((int128)tmp_q[3] * 891133241793L) - ((int128)tmp_q[4] * 1285180556663L) + ((int128)tmp_q[5] * 3301406012735L) + ((int128)tmp_q[6] * 1992425239574L) + ((int128)tmp_q[7] * 3410169732919L) - ((int128)tmp_q[8] * 500176616740L) + ((int128)tmp_q[9] * 7734654581252L) + ((int128)tmp_q[10] * 757732270717L) - ((int128)tmp_q[11] * 644053800861L);

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

