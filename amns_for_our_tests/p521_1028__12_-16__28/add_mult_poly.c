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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) << 4);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) << 4);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) << 4);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) << 4);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) << 4);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) << 4);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) << 4);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) << 4);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) << 4);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) << 4);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) << 4);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) << 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) << 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) << 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) << 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) << 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) << 4);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10060942121681104165UL) + ((((uint64_t)op[1] * 6030321762404648907UL) + ((uint64_t)op[2] * 13565664018539748236UL) + ((uint64_t)op[3] * 14382504041555547534UL) + ((uint64_t)op[4] * 8639521636822833954UL) + ((uint64_t)op[5] * 8579763471416540542UL) + ((uint64_t)op[6] * 2282533160865985468UL) + ((uint64_t)op[7] * 8644937727222260509UL) + ((uint64_t)op[8] * 3701840624430744017UL) + ((uint64_t)op[9] * 9337075829208288084UL) + ((uint64_t)op[10] * 17694152710284881901UL) + ((uint64_t)op[11] * 428995828988955876UL)) * 18446744073709551600);
	tmp_q[1] = ((uint64_t)op[0] * 428995828988955876UL) + ((uint64_t)op[1] * 10060942121681104165UL) + ((((uint64_t)op[2] * 6030321762404648907UL) + ((uint64_t)op[3] * 13565664018539748236UL) + ((uint64_t)op[4] * 14382504041555547534UL) + ((uint64_t)op[5] * 8639521636822833954UL) + ((uint64_t)op[6] * 8579763471416540542UL) + ((uint64_t)op[7] * 2282533160865985468UL) + ((uint64_t)op[8] * 8644937727222260509UL) + ((uint64_t)op[9] * 3701840624430744017UL) + ((uint64_t)op[10] * 9337075829208288084UL) + ((uint64_t)op[11] * 17694152710284881901UL)) * 18446744073709551600);
	tmp_q[2] = ((uint64_t)op[0] * 17694152710284881901UL) + ((uint64_t)op[1] * 428995828988955876UL) + ((uint64_t)op[2] * 10060942121681104165UL) + ((((uint64_t)op[3] * 6030321762404648907UL) + ((uint64_t)op[4] * 13565664018539748236UL) + ((uint64_t)op[5] * 14382504041555547534UL) + ((uint64_t)op[6] * 8639521636822833954UL) + ((uint64_t)op[7] * 8579763471416540542UL) + ((uint64_t)op[8] * 2282533160865985468UL) + ((uint64_t)op[9] * 8644937727222260509UL) + ((uint64_t)op[10] * 3701840624430744017UL) + ((uint64_t)op[11] * 9337075829208288084UL)) * 18446744073709551600);
	tmp_q[3] = ((uint64_t)op[0] * 9337075829208288084UL) + ((uint64_t)op[1] * 17694152710284881901UL) + ((uint64_t)op[2] * 428995828988955876UL) + ((uint64_t)op[3] * 10060942121681104165UL) + ((((uint64_t)op[4] * 6030321762404648907UL) + ((uint64_t)op[5] * 13565664018539748236UL) + ((uint64_t)op[6] * 14382504041555547534UL) + ((uint64_t)op[7] * 8639521636822833954UL) + ((uint64_t)op[8] * 8579763471416540542UL) + ((uint64_t)op[9] * 2282533160865985468UL) + ((uint64_t)op[10] * 8644937727222260509UL) + ((uint64_t)op[11] * 3701840624430744017UL)) * 18446744073709551600);
	tmp_q[4] = ((uint64_t)op[0] * 3701840624430744017UL) + ((uint64_t)op[1] * 9337075829208288084UL) + ((uint64_t)op[2] * 17694152710284881901UL) + ((uint64_t)op[3] * 428995828988955876UL) + ((uint64_t)op[4] * 10060942121681104165UL) + ((((uint64_t)op[5] * 6030321762404648907UL) + ((uint64_t)op[6] * 13565664018539748236UL) + ((uint64_t)op[7] * 14382504041555547534UL) + ((uint64_t)op[8] * 8639521636822833954UL) + ((uint64_t)op[9] * 8579763471416540542UL) + ((uint64_t)op[10] * 2282533160865985468UL) + ((uint64_t)op[11] * 8644937727222260509UL)) * 18446744073709551600);
	tmp_q[5] = ((uint64_t)op[0] * 8644937727222260509UL) + ((uint64_t)op[1] * 3701840624430744017UL) + ((uint64_t)op[2] * 9337075829208288084UL) + ((uint64_t)op[3] * 17694152710284881901UL) + ((uint64_t)op[4] * 428995828988955876UL) + ((uint64_t)op[5] * 10060942121681104165UL) + ((((uint64_t)op[6] * 6030321762404648907UL) + ((uint64_t)op[7] * 13565664018539748236UL) + ((uint64_t)op[8] * 14382504041555547534UL) + ((uint64_t)op[9] * 8639521636822833954UL) + ((uint64_t)op[10] * 8579763471416540542UL) + ((uint64_t)op[11] * 2282533160865985468UL)) * 18446744073709551600);
	tmp_q[6] = ((uint64_t)op[0] * 2282533160865985468UL) + ((uint64_t)op[1] * 8644937727222260509UL) + ((uint64_t)op[2] * 3701840624430744017UL) + ((uint64_t)op[3] * 9337075829208288084UL) + ((uint64_t)op[4] * 17694152710284881901UL) + ((uint64_t)op[5] * 428995828988955876UL) + ((uint64_t)op[6] * 10060942121681104165UL) + ((((uint64_t)op[7] * 6030321762404648907UL) + ((uint64_t)op[8] * 13565664018539748236UL) + ((uint64_t)op[9] * 14382504041555547534UL) + ((uint64_t)op[10] * 8639521636822833954UL) + ((uint64_t)op[11] * 8579763471416540542UL)) * 18446744073709551600);
	tmp_q[7] = ((uint64_t)op[0] * 8579763471416540542UL) + ((uint64_t)op[1] * 2282533160865985468UL) + ((uint64_t)op[2] * 8644937727222260509UL) + ((uint64_t)op[3] * 3701840624430744017UL) + ((uint64_t)op[4] * 9337075829208288084UL) + ((uint64_t)op[5] * 17694152710284881901UL) + ((uint64_t)op[6] * 428995828988955876UL) + ((uint64_t)op[7] * 10060942121681104165UL) + ((((uint64_t)op[8] * 6030321762404648907UL) + ((uint64_t)op[9] * 13565664018539748236UL) + ((uint64_t)op[10] * 14382504041555547534UL) + ((uint64_t)op[11] * 8639521636822833954UL)) * 18446744073709551600);
	tmp_q[8] = ((uint64_t)op[0] * 8639521636822833954UL) + ((uint64_t)op[1] * 8579763471416540542UL) + ((uint64_t)op[2] * 2282533160865985468UL) + ((uint64_t)op[3] * 8644937727222260509UL) + ((uint64_t)op[4] * 3701840624430744017UL) + ((uint64_t)op[5] * 9337075829208288084UL) + ((uint64_t)op[6] * 17694152710284881901UL) + ((uint64_t)op[7] * 428995828988955876UL) + ((uint64_t)op[8] * 10060942121681104165UL) + ((((uint64_t)op[9] * 6030321762404648907UL) + ((uint64_t)op[10] * 13565664018539748236UL) + ((uint64_t)op[11] * 14382504041555547534UL)) * 18446744073709551600);
	tmp_q[9] = ((uint64_t)op[0] * 14382504041555547534UL) + ((uint64_t)op[1] * 8639521636822833954UL) + ((uint64_t)op[2] * 8579763471416540542UL) + ((uint64_t)op[3] * 2282533160865985468UL) + ((uint64_t)op[4] * 8644937727222260509UL) + ((uint64_t)op[5] * 3701840624430744017UL) + ((uint64_t)op[6] * 9337075829208288084UL) + ((uint64_t)op[7] * 17694152710284881901UL) + ((uint64_t)op[8] * 428995828988955876UL) + ((uint64_t)op[9] * 10060942121681104165UL) + ((((uint64_t)op[10] * 6030321762404648907UL) + ((uint64_t)op[11] * 13565664018539748236UL)) * 18446744073709551600);
	tmp_q[10] = ((uint64_t)op[0] * 13565664018539748236UL) + ((uint64_t)op[1] * 14382504041555547534UL) + ((uint64_t)op[2] * 8639521636822833954UL) + ((uint64_t)op[3] * 8579763471416540542UL) + ((uint64_t)op[4] * 2282533160865985468UL) + ((uint64_t)op[5] * 8644937727222260509UL) + ((uint64_t)op[6] * 3701840624430744017UL) + ((uint64_t)op[7] * 9337075829208288084UL) + ((uint64_t)op[8] * 17694152710284881901UL) + ((uint64_t)op[9] * 428995828988955876UL) + ((uint64_t)op[10] * 10060942121681104165UL) + ((uint64_t)op[11] * 14195316243782927184UL);
	tmp_q[11] = ((uint64_t)op[0] * 6030321762404648907UL) + ((uint64_t)op[1] * 13565664018539748236UL) + ((uint64_t)op[2] * 14382504041555547534UL) + ((uint64_t)op[3] * 8639521636822833954UL) + ((uint64_t)op[4] * 8579763471416540542UL) + ((uint64_t)op[5] * 2282533160865985468UL) + ((uint64_t)op[6] * 8644937727222260509UL) + ((uint64_t)op[7] * 3701840624430744017UL) + ((uint64_t)op[8] * 9337075829208288084UL) + ((uint64_t)op[9] * 17694152710284881901UL) + ((uint64_t)op[10] * 428995828988955876UL) + ((uint64_t)op[11] * 10060942121681104165UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1107724485587L) - ((-((int128)tmp_q[1] * 3290925942941L) + ((int128)tmp_q[2] * 1785981438415L) - ((int128)tmp_q[3] * 1795927443961L) + ((int128)tmp_q[4] * 1167417851443L) + ((int128)tmp_q[5] * 3988619881376L) - ((int128)tmp_q[6] * 273717809225L) - ((int128)tmp_q[7] * 6842206356495L) - ((int128)tmp_q[8] * 6620660170100L) - ((int128)tmp_q[9] * 2404246732660L) + ((int128)tmp_q[10] * 5370016043269L) + ((int128)tmp_q[11] * 1229503114692L)) * 16);
	tmp_zero[1] = ((int128)tmp_q[0] * 1229503114692L) + ((int128)tmp_q[1] * 1107724485587L) - ((-((int128)tmp_q[2] * 3290925942941L) + ((int128)tmp_q[3] * 1785981438415L) - ((int128)tmp_q[4] * 1795927443961L) + ((int128)tmp_q[5] * 1167417851443L) + ((int128)tmp_q[6] * 3988619881376L) - ((int128)tmp_q[7] * 273717809225L) - ((int128)tmp_q[8] * 6842206356495L) - ((int128)tmp_q[9] * 6620660170100L) - ((int128)tmp_q[10] * 2404246732660L) + ((int128)tmp_q[11] * 5370016043269L)) * 16);
	tmp_zero[2] = ((int128)tmp_q[0] * 5370016043269L) + ((int128)tmp_q[1] * 1229503114692L) + ((int128)tmp_q[2] * 1107724485587L) - ((-((int128)tmp_q[3] * 3290925942941L) + ((int128)tmp_q[4] * 1785981438415L) - ((int128)tmp_q[5] * 1795927443961L) + ((int128)tmp_q[6] * 1167417851443L) + ((int128)tmp_q[7] * 3988619881376L) - ((int128)tmp_q[8] * 273717809225L) - ((int128)tmp_q[9] * 6842206356495L) - ((int128)tmp_q[10] * 6620660170100L) - ((int128)tmp_q[11] * 2404246732660L)) * 16);
	tmp_zero[3] = -((int128)tmp_q[0] * 2404246732660L) + ((int128)tmp_q[1] * 5370016043269L) + ((int128)tmp_q[2] * 1229503114692L) + ((int128)tmp_q[3] * 1107724485587L) - ((-((int128)tmp_q[4] * 3290925942941L) + ((int128)tmp_q[5] * 1785981438415L) - ((int128)tmp_q[6] * 1795927443961L) + ((int128)tmp_q[7] * 1167417851443L) + ((int128)tmp_q[8] * 3988619881376L) - ((int128)tmp_q[9] * 273717809225L) - ((int128)tmp_q[10] * 6842206356495L) - ((int128)tmp_q[11] * 6620660170100L)) * 16);
	tmp_zero[4] = -((int128)tmp_q[0] * 6620660170100L) - ((int128)tmp_q[1] * 2404246732660L) + ((int128)tmp_q[2] * 5370016043269L) + ((int128)tmp_q[3] * 1229503114692L) + ((int128)tmp_q[4] * 1107724485587L) - ((-((int128)tmp_q[5] * 3290925942941L) + ((int128)tmp_q[6] * 1785981438415L) - ((int128)tmp_q[7] * 1795927443961L) + ((int128)tmp_q[8] * 1167417851443L) + ((int128)tmp_q[9] * 3988619881376L) - ((int128)tmp_q[10] * 273717809225L) - ((int128)tmp_q[11] * 6842206356495L)) * 16);
	tmp_zero[5] = -((int128)tmp_q[0] * 6842206356495L) - ((int128)tmp_q[1] * 6620660170100L) - ((int128)tmp_q[2] * 2404246732660L) + ((int128)tmp_q[3] * 5370016043269L) + ((int128)tmp_q[4] * 1229503114692L) + ((int128)tmp_q[5] * 1107724485587L) - ((-((int128)tmp_q[6] * 3290925942941L) + ((int128)tmp_q[7] * 1785981438415L) - ((int128)tmp_q[8] * 1795927443961L) + ((int128)tmp_q[9] * 1167417851443L) + ((int128)tmp_q[10] * 3988619881376L) - ((int128)tmp_q[11] * 273717809225L)) * 16);
	tmp_zero[6] = -((int128)tmp_q[0] * 273717809225L) - ((int128)tmp_q[1] * 6842206356495L) - ((int128)tmp_q[2] * 6620660170100L) - ((int128)tmp_q[3] * 2404246732660L) + ((int128)tmp_q[4] * 5370016043269L) + ((int128)tmp_q[5] * 1229503114692L) + ((int128)tmp_q[6] * 1107724485587L) - ((-((int128)tmp_q[7] * 3290925942941L) + ((int128)tmp_q[8] * 1785981438415L) - ((int128)tmp_q[9] * 1795927443961L) + ((int128)tmp_q[10] * 1167417851443L) + ((int128)tmp_q[11] * 3988619881376L)) * 16);
	tmp_zero[7] = ((int128)tmp_q[0] * 3988619881376L) - ((int128)tmp_q[1] * 273717809225L) - ((int128)tmp_q[2] * 6842206356495L) - ((int128)tmp_q[3] * 6620660170100L) - ((int128)tmp_q[4] * 2404246732660L) + ((int128)tmp_q[5] * 5370016043269L) + ((int128)tmp_q[6] * 1229503114692L) + ((int128)tmp_q[7] * 1107724485587L) - ((-((int128)tmp_q[8] * 3290925942941L) + ((int128)tmp_q[9] * 1785981438415L) - ((int128)tmp_q[10] * 1795927443961L) + ((int128)tmp_q[11] * 1167417851443L)) * 16);
	tmp_zero[8] = ((int128)tmp_q[0] * 1167417851443L) + ((int128)tmp_q[1] * 3988619881376L) - ((int128)tmp_q[2] * 273717809225L) - ((int128)tmp_q[3] * 6842206356495L) - ((int128)tmp_q[4] * 6620660170100L) - ((int128)tmp_q[5] * 2404246732660L) + ((int128)tmp_q[6] * 5370016043269L) + ((int128)tmp_q[7] * 1229503114692L) + ((int128)tmp_q[8] * 1107724485587L) - ((-((int128)tmp_q[9] * 3290925942941L) + ((int128)tmp_q[10] * 1785981438415L) - ((int128)tmp_q[11] * 1795927443961L)) * 16);
	tmp_zero[9] = -((int128)tmp_q[0] * 1795927443961L) + ((int128)tmp_q[1] * 1167417851443L) + ((int128)tmp_q[2] * 3988619881376L) - ((int128)tmp_q[3] * 273717809225L) - ((int128)tmp_q[4] * 6842206356495L) - ((int128)tmp_q[5] * 6620660170100L) - ((int128)tmp_q[6] * 2404246732660L) + ((int128)tmp_q[7] * 5370016043269L) + ((int128)tmp_q[8] * 1229503114692L) + ((int128)tmp_q[9] * 1107724485587L) - ((-((int128)tmp_q[10] * 3290925942941L) + ((int128)tmp_q[11] * 1785981438415L)) * 16);
	tmp_zero[10] = ((int128)tmp_q[0] * 1785981438415L) - ((int128)tmp_q[1] * 1795927443961L) + ((int128)tmp_q[2] * 1167417851443L) + ((int128)tmp_q[3] * 3988619881376L) - ((int128)tmp_q[4] * 273717809225L) - ((int128)tmp_q[5] * 6842206356495L) - ((int128)tmp_q[6] * 6620660170100L) - ((int128)tmp_q[7] * 2404246732660L) + ((int128)tmp_q[8] * 5370016043269L) + ((int128)tmp_q[9] * 1229503114692L) + ((int128)tmp_q[10] * 1107724485587L) + ((int128)tmp_q[11] * 52654815087056L);
	tmp_zero[11] = -((int128)tmp_q[0] * 3290925942941L) + ((int128)tmp_q[1] * 1785981438415L) - ((int128)tmp_q[2] * 1795927443961L) + ((int128)tmp_q[3] * 1167417851443L) + ((int128)tmp_q[4] * 3988619881376L) - ((int128)tmp_q[5] * 273717809225L) - ((int128)tmp_q[6] * 6842206356495L) - ((int128)tmp_q[7] * 6620660170100L) - ((int128)tmp_q[8] * 2404246732660L) + ((int128)tmp_q[9] * 5370016043269L) + ((int128)tmp_q[10] * 1229503114692L) + ((int128)tmp_q[11] * 1107724485587L);

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

