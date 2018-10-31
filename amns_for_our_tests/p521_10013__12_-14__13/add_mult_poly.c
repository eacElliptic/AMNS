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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6503687679131644735UL) + ((((uint64_t)op[1] * 8203184972811538585UL) + ((uint64_t)op[2] * 14274373566720472531UL) + ((uint64_t)op[3] * 14331108768534316570UL) + ((uint64_t)op[4] * 16813954280956205044UL) + ((uint64_t)op[5] * 3461274974002850526UL) + ((uint64_t)op[6] * 1498261524989162869UL) + ((uint64_t)op[7] * 8174790301846847154UL) + ((uint64_t)op[8] * 12323962059799690213UL) + ((uint64_t)op[9] * 10853367982675705306UL) + ((uint64_t)op[10] * 3959950905931286400UL) + ((uint64_t)op[11] * 14574864607800013350UL)) * 18446744073709551602);
	tmp_q[1] = ((uint64_t)op[0] * 14574864607800013350UL) + ((uint64_t)op[1] * 6503687679131644735UL) + ((((uint64_t)op[2] * 8203184972811538585UL) + ((uint64_t)op[3] * 14274373566720472531UL) + ((uint64_t)op[4] * 14331108768534316570UL) + ((uint64_t)op[5] * 16813954280956205044UL) + ((uint64_t)op[6] * 3461274974002850526UL) + ((uint64_t)op[7] * 1498261524989162869UL) + ((uint64_t)op[8] * 8174790301846847154UL) + ((uint64_t)op[9] * 12323962059799690213UL) + ((uint64_t)op[10] * 10853367982675705306UL) + ((uint64_t)op[11] * 3959950905931286400UL)) * 18446744073709551602);
	tmp_q[2] = ((uint64_t)op[0] * 3959950905931286400UL) + ((uint64_t)op[1] * 14574864607800013350UL) + ((uint64_t)op[2] * 6503687679131644735UL) + ((((uint64_t)op[3] * 8203184972811538585UL) + ((uint64_t)op[4] * 14274373566720472531UL) + ((uint64_t)op[5] * 14331108768534316570UL) + ((uint64_t)op[6] * 16813954280956205044UL) + ((uint64_t)op[7] * 3461274974002850526UL) + ((uint64_t)op[8] * 1498261524989162869UL) + ((uint64_t)op[9] * 8174790301846847154UL) + ((uint64_t)op[10] * 12323962059799690213UL) + ((uint64_t)op[11] * 10853367982675705306UL)) * 18446744073709551602);
	tmp_q[3] = ((uint64_t)op[0] * 10853367982675705306UL) + ((uint64_t)op[1] * 3959950905931286400UL) + ((uint64_t)op[2] * 14574864607800013350UL) + ((uint64_t)op[3] * 6503687679131644735UL) + ((((uint64_t)op[4] * 8203184972811538585UL) + ((uint64_t)op[5] * 14274373566720472531UL) + ((uint64_t)op[6] * 14331108768534316570UL) + ((uint64_t)op[7] * 16813954280956205044UL) + ((uint64_t)op[8] * 3461274974002850526UL) + ((uint64_t)op[9] * 1498261524989162869UL) + ((uint64_t)op[10] * 8174790301846847154UL) + ((uint64_t)op[11] * 12323962059799690213UL)) * 18446744073709551602);
	tmp_q[4] = ((uint64_t)op[0] * 12323962059799690213UL) + ((uint64_t)op[1] * 10853367982675705306UL) + ((uint64_t)op[2] * 3959950905931286400UL) + ((uint64_t)op[3] * 14574864607800013350UL) + ((uint64_t)op[4] * 6503687679131644735UL) + ((((uint64_t)op[5] * 8203184972811538585UL) + ((uint64_t)op[6] * 14274373566720472531UL) + ((uint64_t)op[7] * 14331108768534316570UL) + ((uint64_t)op[8] * 16813954280956205044UL) + ((uint64_t)op[9] * 3461274974002850526UL) + ((uint64_t)op[10] * 1498261524989162869UL) + ((uint64_t)op[11] * 8174790301846847154UL)) * 18446744073709551602);
	tmp_q[5] = ((uint64_t)op[0] * 8174790301846847154UL) + ((uint64_t)op[1] * 12323962059799690213UL) + ((uint64_t)op[2] * 10853367982675705306UL) + ((uint64_t)op[3] * 3959950905931286400UL) + ((uint64_t)op[4] * 14574864607800013350UL) + ((uint64_t)op[5] * 6503687679131644735UL) + ((((uint64_t)op[6] * 8203184972811538585UL) + ((uint64_t)op[7] * 14274373566720472531UL) + ((uint64_t)op[8] * 14331108768534316570UL) + ((uint64_t)op[9] * 16813954280956205044UL) + ((uint64_t)op[10] * 3461274974002850526UL) + ((uint64_t)op[11] * 1498261524989162869UL)) * 18446744073709551602);
	tmp_q[6] = ((uint64_t)op[0] * 1498261524989162869UL) + ((uint64_t)op[1] * 8174790301846847154UL) + ((uint64_t)op[2] * 12323962059799690213UL) + ((uint64_t)op[3] * 10853367982675705306UL) + ((uint64_t)op[4] * 3959950905931286400UL) + ((uint64_t)op[5] * 14574864607800013350UL) + ((uint64_t)op[6] * 6503687679131644735UL) + ((((uint64_t)op[7] * 8203184972811538585UL) + ((uint64_t)op[8] * 14274373566720472531UL) + ((uint64_t)op[9] * 14331108768534316570UL) + ((uint64_t)op[10] * 16813954280956205044UL) + ((uint64_t)op[11] * 3461274974002850526UL)) * 18446744073709551602);
	tmp_q[7] = ((uint64_t)op[0] * 3461274974002850526UL) + ((uint64_t)op[1] * 1498261524989162869UL) + ((uint64_t)op[2] * 8174790301846847154UL) + ((uint64_t)op[3] * 12323962059799690213UL) + ((uint64_t)op[4] * 10853367982675705306UL) + ((uint64_t)op[5] * 3959950905931286400UL) + ((uint64_t)op[6] * 14574864607800013350UL) + ((uint64_t)op[7] * 6503687679131644735UL) + ((((uint64_t)op[8] * 8203184972811538585UL) + ((uint64_t)op[9] * 14274373566720472531UL) + ((uint64_t)op[10] * 14331108768534316570UL) + ((uint64_t)op[11] * 16813954280956205044UL)) * 18446744073709551602);
	tmp_q[8] = ((uint64_t)op[0] * 16813954280956205044UL) + ((uint64_t)op[1] * 3461274974002850526UL) + ((uint64_t)op[2] * 1498261524989162869UL) + ((uint64_t)op[3] * 8174790301846847154UL) + ((uint64_t)op[4] * 12323962059799690213UL) + ((uint64_t)op[5] * 10853367982675705306UL) + ((uint64_t)op[6] * 3959950905931286400UL) + ((uint64_t)op[7] * 14574864607800013350UL) + ((uint64_t)op[8] * 6503687679131644735UL) + ((((uint64_t)op[9] * 8203184972811538585UL) + ((uint64_t)op[10] * 14274373566720472531UL) + ((uint64_t)op[11] * 14331108768534316570UL)) * 18446744073709551602);
	tmp_q[9] = ((uint64_t)op[0] * 14331108768534316570UL) + ((uint64_t)op[1] * 16813954280956205044UL) + ((uint64_t)op[2] * 3461274974002850526UL) + ((uint64_t)op[3] * 1498261524989162869UL) + ((uint64_t)op[4] * 8174790301846847154UL) + ((uint64_t)op[5] * 12323962059799690213UL) + ((uint64_t)op[6] * 10853367982675705306UL) + ((uint64_t)op[7] * 3959950905931286400UL) + ((uint64_t)op[8] * 14574864607800013350UL) + ((uint64_t)op[9] * 6503687679131644735UL) + ((((uint64_t)op[10] * 8203184972811538585UL) + ((uint64_t)op[11] * 14274373566720472531UL)) * 18446744073709551602);
	tmp_q[10] = ((uint64_t)op[0] * 14274373566720472531UL) + ((uint64_t)op[1] * 14331108768534316570UL) + ((uint64_t)op[2] * 16813954280956205044UL) + ((uint64_t)op[3] * 3461274974002850526UL) + ((uint64_t)op[4] * 1498261524989162869UL) + ((uint64_t)op[5] * 8174790301846847154UL) + ((uint64_t)op[6] * 12323962059799690213UL) + ((uint64_t)op[7] * 10853367982675705306UL) + ((uint64_t)op[8] * 3959950905931286400UL) + ((uint64_t)op[9] * 14574864607800013350UL) + ((uint64_t)op[10] * 6503687679131644735UL) + ((uint64_t)op[11] * 14282618896605321122UL);
	tmp_q[11] = ((uint64_t)op[0] * 8203184972811538585UL) + ((uint64_t)op[1] * 14274373566720472531UL) + ((uint64_t)op[2] * 14331108768534316570UL) + ((uint64_t)op[3] * 16813954280956205044UL) + ((uint64_t)op[4] * 3461274974002850526UL) + ((uint64_t)op[5] * 1498261524989162869UL) + ((uint64_t)op[6] * 8174790301846847154UL) + ((uint64_t)op[7] * 12323962059799690213UL) + ((uint64_t)op[8] * 10853367982675705306UL) + ((uint64_t)op[9] * 3959950905931286400UL) + ((uint64_t)op[10] * 14574864607800013350UL) + ((uint64_t)op[11] * 6503687679131644735UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2245198196175L) - ((((int128)tmp_q[1] * 3787097022577L) + ((int128)tmp_q[2] * 3689777337111L) - ((int128)tmp_q[3] * 1799711524612L) + ((int128)tmp_q[4] * 3892794976093L) - ((int128)tmp_q[5] * 2909872507096L) - ((int128)tmp_q[6] * 6924463652147L) - ((int128)tmp_q[7] * 89716629698L) - ((int128)tmp_q[8] * 5870664435607L) + ((int128)tmp_q[9] * 3354785155162L) + ((int128)tmp_q[10] * 928054669338L) + ((int128)tmp_q[11] * 4585266827502L)) * 14);
	tmp_zero[1] = ((int128)tmp_q[0] * 4585266827502L) - ((int128)tmp_q[1] * 2245198196175L) - ((((int128)tmp_q[2] * 3787097022577L) + ((int128)tmp_q[3] * 3689777337111L) - ((int128)tmp_q[4] * 1799711524612L) + ((int128)tmp_q[5] * 3892794976093L) - ((int128)tmp_q[6] * 2909872507096L) - ((int128)tmp_q[7] * 6924463652147L) - ((int128)tmp_q[8] * 89716629698L) - ((int128)tmp_q[9] * 5870664435607L) + ((int128)tmp_q[10] * 3354785155162L) + ((int128)tmp_q[11] * 928054669338L)) * 14);
	tmp_zero[2] = ((int128)tmp_q[0] * 928054669338L) + ((int128)tmp_q[1] * 4585266827502L) - ((int128)tmp_q[2] * 2245198196175L) - ((((int128)tmp_q[3] * 3787097022577L) + ((int128)tmp_q[4] * 3689777337111L) - ((int128)tmp_q[5] * 1799711524612L) + ((int128)tmp_q[6] * 3892794976093L) - ((int128)tmp_q[7] * 2909872507096L) - ((int128)tmp_q[8] * 6924463652147L) - ((int128)tmp_q[9] * 89716629698L) - ((int128)tmp_q[10] * 5870664435607L) + ((int128)tmp_q[11] * 3354785155162L)) * 14);
	tmp_zero[3] = ((int128)tmp_q[0] * 3354785155162L) + ((int128)tmp_q[1] * 928054669338L) + ((int128)tmp_q[2] * 4585266827502L) - ((int128)tmp_q[3] * 2245198196175L) - ((((int128)tmp_q[4] * 3787097022577L) + ((int128)tmp_q[5] * 3689777337111L) - ((int128)tmp_q[6] * 1799711524612L) + ((int128)tmp_q[7] * 3892794976093L) - ((int128)tmp_q[8] * 2909872507096L) - ((int128)tmp_q[9] * 6924463652147L) - ((int128)tmp_q[10] * 89716629698L) - ((int128)tmp_q[11] * 5870664435607L)) * 14);
	tmp_zero[4] = -((int128)tmp_q[0] * 5870664435607L) + ((int128)tmp_q[1] * 3354785155162L) + ((int128)tmp_q[2] * 928054669338L) + ((int128)tmp_q[3] * 4585266827502L) - ((int128)tmp_q[4] * 2245198196175L) - ((((int128)tmp_q[5] * 3787097022577L) + ((int128)tmp_q[6] * 3689777337111L) - ((int128)tmp_q[7] * 1799711524612L) + ((int128)tmp_q[8] * 3892794976093L) - ((int128)tmp_q[9] * 2909872507096L) - ((int128)tmp_q[10] * 6924463652147L) - ((int128)tmp_q[11] * 89716629698L)) * 14);
	tmp_zero[5] = -((int128)tmp_q[0] * 89716629698L) - ((int128)tmp_q[1] * 5870664435607L) + ((int128)tmp_q[2] * 3354785155162L) + ((int128)tmp_q[3] * 928054669338L) + ((int128)tmp_q[4] * 4585266827502L) - ((int128)tmp_q[5] * 2245198196175L) - ((((int128)tmp_q[6] * 3787097022577L) + ((int128)tmp_q[7] * 3689777337111L) - ((int128)tmp_q[8] * 1799711524612L) + ((int128)tmp_q[9] * 3892794976093L) - ((int128)tmp_q[10] * 2909872507096L) - ((int128)tmp_q[11] * 6924463652147L)) * 14);
	tmp_zero[6] = -((int128)tmp_q[0] * 6924463652147L) - ((int128)tmp_q[1] * 89716629698L) - ((int128)tmp_q[2] * 5870664435607L) + ((int128)tmp_q[3] * 3354785155162L) + ((int128)tmp_q[4] * 928054669338L) + ((int128)tmp_q[5] * 4585266827502L) - ((int128)tmp_q[6] * 2245198196175L) - ((((int128)tmp_q[7] * 3787097022577L) + ((int128)tmp_q[8] * 3689777337111L) - ((int128)tmp_q[9] * 1799711524612L) + ((int128)tmp_q[10] * 3892794976093L) - ((int128)tmp_q[11] * 2909872507096L)) * 14);
	tmp_zero[7] = -((int128)tmp_q[0] * 2909872507096L) - ((int128)tmp_q[1] * 6924463652147L) - ((int128)tmp_q[2] * 89716629698L) - ((int128)tmp_q[3] * 5870664435607L) + ((int128)tmp_q[4] * 3354785155162L) + ((int128)tmp_q[5] * 928054669338L) + ((int128)tmp_q[6] * 4585266827502L) - ((int128)tmp_q[7] * 2245198196175L) - ((((int128)tmp_q[8] * 3787097022577L) + ((int128)tmp_q[9] * 3689777337111L) - ((int128)tmp_q[10] * 1799711524612L) + ((int128)tmp_q[11] * 3892794976093L)) * 14);
	tmp_zero[8] = ((int128)tmp_q[0] * 3892794976093L) - ((int128)tmp_q[1] * 2909872507096L) - ((int128)tmp_q[2] * 6924463652147L) - ((int128)tmp_q[3] * 89716629698L) - ((int128)tmp_q[4] * 5870664435607L) + ((int128)tmp_q[5] * 3354785155162L) + ((int128)tmp_q[6] * 928054669338L) + ((int128)tmp_q[7] * 4585266827502L) - ((int128)tmp_q[8] * 2245198196175L) - ((((int128)tmp_q[9] * 3787097022577L) + ((int128)tmp_q[10] * 3689777337111L) - ((int128)tmp_q[11] * 1799711524612L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 1799711524612L) + ((int128)tmp_q[1] * 3892794976093L) - ((int128)tmp_q[2] * 2909872507096L) - ((int128)tmp_q[3] * 6924463652147L) - ((int128)tmp_q[4] * 89716629698L) - ((int128)tmp_q[5] * 5870664435607L) + ((int128)tmp_q[6] * 3354785155162L) + ((int128)tmp_q[7] * 928054669338L) + ((int128)tmp_q[8] * 4585266827502L) - ((int128)tmp_q[9] * 2245198196175L) - ((((int128)tmp_q[10] * 3787097022577L) + ((int128)tmp_q[11] * 3689777337111L)) * 14);
	tmp_zero[10] = ((int128)tmp_q[0] * 3689777337111L) - ((int128)tmp_q[1] * 1799711524612L) + ((int128)tmp_q[2] * 3892794976093L) - ((int128)tmp_q[3] * 2909872507096L) - ((int128)tmp_q[4] * 6924463652147L) - ((int128)tmp_q[5] * 89716629698L) - ((int128)tmp_q[6] * 5870664435607L) + ((int128)tmp_q[7] * 3354785155162L) + ((int128)tmp_q[8] * 928054669338L) + ((int128)tmp_q[9] * 4585266827502L) - ((int128)tmp_q[10] * 2245198196175L) - ((int128)tmp_q[11] * 53019358316078L);
	tmp_zero[11] = ((int128)tmp_q[0] * 3787097022577L) + ((int128)tmp_q[1] * 3689777337111L) - ((int128)tmp_q[2] * 1799711524612L) + ((int128)tmp_q[3] * 3892794976093L) - ((int128)tmp_q[4] * 2909872507096L) - ((int128)tmp_q[5] * 6924463652147L) - ((int128)tmp_q[6] * 89716629698L) - ((int128)tmp_q[7] * 5870664435607L) + ((int128)tmp_q[8] * 3354785155162L) + ((int128)tmp_q[9] * 928054669338L) + ((int128)tmp_q[10] * 4585266827502L) - ((int128)tmp_q[11] * 2245198196175L);

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

