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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11118434160020890875UL) + ((((uint64_t)op[1] * 14370409012785167890UL) + ((uint64_t)op[2] * 7166158576714992454UL) + ((uint64_t)op[3] * 12204392295438114394UL) + ((uint64_t)op[4] * 3215697565427119172UL) + ((uint64_t)op[5] * 12167477915755237047UL) + ((uint64_t)op[6] * 15130946489549850276UL) + ((uint64_t)op[7] * 14947172906696856941UL) + ((uint64_t)op[8] * 13898040695765672688UL) + ((uint64_t)op[9] * 15208459982896985771UL) + ((uint64_t)op[10] * 7092585686045929206UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 7092585686045929206UL) + ((uint64_t)op[1] * 11118434160020890875UL) + ((((uint64_t)op[2] * 14370409012785167890UL) + ((uint64_t)op[3] * 7166158576714992454UL) + ((uint64_t)op[4] * 12204392295438114394UL) + ((uint64_t)op[5] * 3215697565427119172UL) + ((uint64_t)op[6] * 12167477915755237047UL) + ((uint64_t)op[7] * 15130946489549850276UL) + ((uint64_t)op[8] * 14947172906696856941UL) + ((uint64_t)op[9] * 13898040695765672688UL) + ((uint64_t)op[10] * 15208459982896985771UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 15208459982896985771UL) + ((uint64_t)op[1] * 7092585686045929206UL) + ((uint64_t)op[2] * 11118434160020890875UL) + ((((uint64_t)op[3] * 14370409012785167890UL) + ((uint64_t)op[4] * 7166158576714992454UL) + ((uint64_t)op[5] * 12204392295438114394UL) + ((uint64_t)op[6] * 3215697565427119172UL) + ((uint64_t)op[7] * 12167477915755237047UL) + ((uint64_t)op[8] * 15130946489549850276UL) + ((uint64_t)op[9] * 14947172906696856941UL) + ((uint64_t)op[10] * 13898040695765672688UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 13898040695765672688UL) + ((uint64_t)op[1] * 15208459982896985771UL) + ((uint64_t)op[2] * 7092585686045929206UL) + ((uint64_t)op[3] * 11118434160020890875UL) + ((((uint64_t)op[4] * 14370409012785167890UL) + ((uint64_t)op[5] * 7166158576714992454UL) + ((uint64_t)op[6] * 12204392295438114394UL) + ((uint64_t)op[7] * 3215697565427119172UL) + ((uint64_t)op[8] * 12167477915755237047UL) + ((uint64_t)op[9] * 15130946489549850276UL) + ((uint64_t)op[10] * 14947172906696856941UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 14947172906696856941UL) + ((uint64_t)op[1] * 13898040695765672688UL) + ((uint64_t)op[2] * 15208459982896985771UL) + ((uint64_t)op[3] * 7092585686045929206UL) + ((uint64_t)op[4] * 11118434160020890875UL) + ((((uint64_t)op[5] * 14370409012785167890UL) + ((uint64_t)op[6] * 7166158576714992454UL) + ((uint64_t)op[7] * 12204392295438114394UL) + ((uint64_t)op[8] * 3215697565427119172UL) + ((uint64_t)op[9] * 12167477915755237047UL) + ((uint64_t)op[10] * 15130946489549850276UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 15130946489549850276UL) + ((uint64_t)op[1] * 14947172906696856941UL) + ((uint64_t)op[2] * 13898040695765672688UL) + ((uint64_t)op[3] * 15208459982896985771UL) + ((uint64_t)op[4] * 7092585686045929206UL) + ((uint64_t)op[5] * 11118434160020890875UL) + ((((uint64_t)op[6] * 14370409012785167890UL) + ((uint64_t)op[7] * 7166158576714992454UL) + ((uint64_t)op[8] * 12204392295438114394UL) + ((uint64_t)op[9] * 3215697565427119172UL) + ((uint64_t)op[10] * 12167477915755237047UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 12167477915755237047UL) + ((uint64_t)op[1] * 15130946489549850276UL) + ((uint64_t)op[2] * 14947172906696856941UL) + ((uint64_t)op[3] * 13898040695765672688UL) + ((uint64_t)op[4] * 15208459982896985771UL) + ((uint64_t)op[5] * 7092585686045929206UL) + ((uint64_t)op[6] * 11118434160020890875UL) + ((((uint64_t)op[7] * 14370409012785167890UL) + ((uint64_t)op[8] * 7166158576714992454UL) + ((uint64_t)op[9] * 12204392295438114394UL) + ((uint64_t)op[10] * 3215697565427119172UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 3215697565427119172UL) + ((uint64_t)op[1] * 12167477915755237047UL) + ((uint64_t)op[2] * 15130946489549850276UL) + ((uint64_t)op[3] * 14947172906696856941UL) + ((uint64_t)op[4] * 13898040695765672688UL) + ((uint64_t)op[5] * 15208459982896985771UL) + ((uint64_t)op[6] * 7092585686045929206UL) + ((uint64_t)op[7] * 11118434160020890875UL) + ((((uint64_t)op[8] * 14370409012785167890UL) + ((uint64_t)op[9] * 7166158576714992454UL) + ((uint64_t)op[10] * 12204392295438114394UL)) * 18446744073709551610);
	tmp_q[8] = ((uint64_t)op[0] * 12204392295438114394UL) + ((uint64_t)op[1] * 3215697565427119172UL) + ((uint64_t)op[2] * 12167477915755237047UL) + ((uint64_t)op[3] * 15130946489549850276UL) + ((uint64_t)op[4] * 14947172906696856941UL) + ((uint64_t)op[5] * 13898040695765672688UL) + ((uint64_t)op[6] * 15208459982896985771UL) + ((uint64_t)op[7] * 7092585686045929206UL) + ((uint64_t)op[8] * 11118434160020890875UL) + ((((uint64_t)op[9] * 14370409012785167890UL) + ((uint64_t)op[10] * 7166158576714992454UL)) * 18446744073709551610);
	tmp_q[9] = ((uint64_t)op[0] * 7166158576714992454UL) + ((uint64_t)op[1] * 12204392295438114394UL) + ((uint64_t)op[2] * 3215697565427119172UL) + ((uint64_t)op[3] * 12167477915755237047UL) + ((uint64_t)op[4] * 15130946489549850276UL) + ((uint64_t)op[5] * 14947172906696856941UL) + ((uint64_t)op[6] * 13898040695765672688UL) + ((uint64_t)op[7] * 15208459982896985771UL) + ((uint64_t)op[8] * 7092585686045929206UL) + ((uint64_t)op[9] * 11118434160020890875UL) + ((uint64_t)op[10] * 6011266291836750740UL);
	tmp_q[10] = ((uint64_t)op[0] * 14370409012785167890UL) + ((uint64_t)op[1] * 7166158576714992454UL) + ((uint64_t)op[2] * 12204392295438114394UL) + ((uint64_t)op[3] * 3215697565427119172UL) + ((uint64_t)op[4] * 12167477915755237047UL) + ((uint64_t)op[5] * 15130946489549850276UL) + ((uint64_t)op[6] * 14947172906696856941UL) + ((uint64_t)op[7] * 13898040695765672688UL) + ((uint64_t)op[8] * 15208459982896985771UL) + ((uint64_t)op[9] * 7092585686045929206UL) + ((uint64_t)op[10] * 11118434160020890875UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 89333904220963L) - ((((int128)tmp_q[1] * 8085010669349L) - ((int128)tmp_q[2] * 77843782700158L) - ((int128)tmp_q[3] * 27099995537727L) + ((int128)tmp_q[4] * 61951466755538L) + ((int128)tmp_q[5] * 76473458177408L) + ((int128)tmp_q[6] * 45057682085804L) + ((int128)tmp_q[7] * 14338499133338L) - ((int128)tmp_q[8] * 86594860238056L) + ((int128)tmp_q[9] * 35051242980527L) + ((int128)tmp_q[10] * 53293311810018L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 53293311810018L) - ((int128)tmp_q[1] * 89333904220963L) - ((((int128)tmp_q[2] * 8085010669349L) - ((int128)tmp_q[3] * 77843782700158L) - ((int128)tmp_q[4] * 27099995537727L) + ((int128)tmp_q[5] * 61951466755538L) + ((int128)tmp_q[6] * 76473458177408L) + ((int128)tmp_q[7] * 45057682085804L) + ((int128)tmp_q[8] * 14338499133338L) - ((int128)tmp_q[9] * 86594860238056L) + ((int128)tmp_q[10] * 35051242980527L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 35051242980527L) + ((int128)tmp_q[1] * 53293311810018L) - ((int128)tmp_q[2] * 89333904220963L) - ((((int128)tmp_q[3] * 8085010669349L) - ((int128)tmp_q[4] * 77843782700158L) - ((int128)tmp_q[5] * 27099995537727L) + ((int128)tmp_q[6] * 61951466755538L) + ((int128)tmp_q[7] * 76473458177408L) + ((int128)tmp_q[8] * 45057682085804L) + ((int128)tmp_q[9] * 14338499133338L) - ((int128)tmp_q[10] * 86594860238056L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 86594860238056L) + ((int128)tmp_q[1] * 35051242980527L) + ((int128)tmp_q[2] * 53293311810018L) - ((int128)tmp_q[3] * 89333904220963L) - ((((int128)tmp_q[4] * 8085010669349L) - ((int128)tmp_q[5] * 77843782700158L) - ((int128)tmp_q[6] * 27099995537727L) + ((int128)tmp_q[7] * 61951466755538L) + ((int128)tmp_q[8] * 76473458177408L) + ((int128)tmp_q[9] * 45057682085804L) + ((int128)tmp_q[10] * 14338499133338L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 14338499133338L) - ((int128)tmp_q[1] * 86594860238056L) + ((int128)tmp_q[2] * 35051242980527L) + ((int128)tmp_q[3] * 53293311810018L) - ((int128)tmp_q[4] * 89333904220963L) - ((((int128)tmp_q[5] * 8085010669349L) - ((int128)tmp_q[6] * 77843782700158L) - ((int128)tmp_q[7] * 27099995537727L) + ((int128)tmp_q[8] * 61951466755538L) + ((int128)tmp_q[9] * 76473458177408L) + ((int128)tmp_q[10] * 45057682085804L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 45057682085804L) + ((int128)tmp_q[1] * 14338499133338L) - ((int128)tmp_q[2] * 86594860238056L) + ((int128)tmp_q[3] * 35051242980527L) + ((int128)tmp_q[4] * 53293311810018L) - ((int128)tmp_q[5] * 89333904220963L) - ((((int128)tmp_q[6] * 8085010669349L) - ((int128)tmp_q[7] * 77843782700158L) - ((int128)tmp_q[8] * 27099995537727L) + ((int128)tmp_q[9] * 61951466755538L) + ((int128)tmp_q[10] * 76473458177408L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 76473458177408L) + ((int128)tmp_q[1] * 45057682085804L) + ((int128)tmp_q[2] * 14338499133338L) - ((int128)tmp_q[3] * 86594860238056L) + ((int128)tmp_q[4] * 35051242980527L) + ((int128)tmp_q[5] * 53293311810018L) - ((int128)tmp_q[6] * 89333904220963L) - ((((int128)tmp_q[7] * 8085010669349L) - ((int128)tmp_q[8] * 77843782700158L) - ((int128)tmp_q[9] * 27099995537727L) + ((int128)tmp_q[10] * 61951466755538L)) * 6);
	tmp_zero[7] = ((int128)tmp_q[0] * 61951466755538L) + ((int128)tmp_q[1] * 76473458177408L) + ((int128)tmp_q[2] * 45057682085804L) + ((int128)tmp_q[3] * 14338499133338L) - ((int128)tmp_q[4] * 86594860238056L) + ((int128)tmp_q[5] * 35051242980527L) + ((int128)tmp_q[6] * 53293311810018L) - ((int128)tmp_q[7] * 89333904220963L) - ((((int128)tmp_q[8] * 8085010669349L) - ((int128)tmp_q[9] * 77843782700158L) - ((int128)tmp_q[10] * 27099995537727L)) * 6);
	tmp_zero[8] = -((int128)tmp_q[0] * 27099995537727L) + ((int128)tmp_q[1] * 61951466755538L) + ((int128)tmp_q[2] * 76473458177408L) + ((int128)tmp_q[3] * 45057682085804L) + ((int128)tmp_q[4] * 14338499133338L) - ((int128)tmp_q[5] * 86594860238056L) + ((int128)tmp_q[6] * 35051242980527L) + ((int128)tmp_q[7] * 53293311810018L) - ((int128)tmp_q[8] * 89333904220963L) - ((((int128)tmp_q[9] * 8085010669349L) - ((int128)tmp_q[10] * 77843782700158L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 77843782700158L) - ((int128)tmp_q[1] * 27099995537727L) + ((int128)tmp_q[2] * 61951466755538L) + ((int128)tmp_q[3] * 76473458177408L) + ((int128)tmp_q[4] * 45057682085804L) + ((int128)tmp_q[5] * 14338499133338L) - ((int128)tmp_q[6] * 86594860238056L) + ((int128)tmp_q[7] * 35051242980527L) + ((int128)tmp_q[8] * 53293311810018L) - ((int128)tmp_q[9] * 89333904220963L) - ((int128)tmp_q[10] * 48510064016094L);
	tmp_zero[10] = ((int128)tmp_q[0] * 8085010669349L) - ((int128)tmp_q[1] * 77843782700158L) - ((int128)tmp_q[2] * 27099995537727L) + ((int128)tmp_q[3] * 61951466755538L) + ((int128)tmp_q[4] * 76473458177408L) + ((int128)tmp_q[5] * 45057682085804L) + ((int128)tmp_q[6] * 14338499133338L) - ((int128)tmp_q[7] * 86594860238056L) + ((int128)tmp_q[8] * 35051242980527L) + ((int128)tmp_q[9] * 53293311810018L) - ((int128)tmp_q[10] * 89333904220963L);

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

