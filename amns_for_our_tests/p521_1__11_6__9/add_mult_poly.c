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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3814736251555248829UL) + ((((uint64_t)op[1] * 7063588192017465263UL) + ((uint64_t)op[2] * 13785887791050718948UL) + ((uint64_t)op[3] * 6489736271999419374UL) + ((uint64_t)op[4] * 7810650901096139443UL) + ((uint64_t)op[5] * 10992916032076318400UL) + ((uint64_t)op[6] * 17459363388283345339UL) + ((uint64_t)op[7] * 13180690705578554535UL) + ((uint64_t)op[8] * 15807971680528568648UL) + ((uint64_t)op[9] * 10679651914351631252UL) + ((uint64_t)op[10] * 12715087827704712706UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 12715087827704712706UL) + ((uint64_t)op[1] * 3814736251555248829UL) + ((((uint64_t)op[2] * 7063588192017465263UL) + ((uint64_t)op[3] * 13785887791050718948UL) + ((uint64_t)op[4] * 6489736271999419374UL) + ((uint64_t)op[5] * 7810650901096139443UL) + ((uint64_t)op[6] * 10992916032076318400UL) + ((uint64_t)op[7] * 17459363388283345339UL) + ((uint64_t)op[8] * 13180690705578554535UL) + ((uint64_t)op[9] * 15807971680528568648UL) + ((uint64_t)op[10] * 10679651914351631252UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 10679651914351631252UL) + ((uint64_t)op[1] * 12715087827704712706UL) + ((uint64_t)op[2] * 3814736251555248829UL) + ((((uint64_t)op[3] * 7063588192017465263UL) + ((uint64_t)op[4] * 13785887791050718948UL) + ((uint64_t)op[5] * 6489736271999419374UL) + ((uint64_t)op[6] * 7810650901096139443UL) + ((uint64_t)op[7] * 10992916032076318400UL) + ((uint64_t)op[8] * 17459363388283345339UL) + ((uint64_t)op[9] * 13180690705578554535UL) + ((uint64_t)op[10] * 15807971680528568648UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 15807971680528568648UL) + ((uint64_t)op[1] * 10679651914351631252UL) + ((uint64_t)op[2] * 12715087827704712706UL) + ((uint64_t)op[3] * 3814736251555248829UL) + ((((uint64_t)op[4] * 7063588192017465263UL) + ((uint64_t)op[5] * 13785887791050718948UL) + ((uint64_t)op[6] * 6489736271999419374UL) + ((uint64_t)op[7] * 7810650901096139443UL) + ((uint64_t)op[8] * 10992916032076318400UL) + ((uint64_t)op[9] * 17459363388283345339UL) + ((uint64_t)op[10] * 13180690705578554535UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 13180690705578554535UL) + ((uint64_t)op[1] * 15807971680528568648UL) + ((uint64_t)op[2] * 10679651914351631252UL) + ((uint64_t)op[3] * 12715087827704712706UL) + ((uint64_t)op[4] * 3814736251555248829UL) + ((((uint64_t)op[5] * 7063588192017465263UL) + ((uint64_t)op[6] * 13785887791050718948UL) + ((uint64_t)op[7] * 6489736271999419374UL) + ((uint64_t)op[8] * 7810650901096139443UL) + ((uint64_t)op[9] * 10992916032076318400UL) + ((uint64_t)op[10] * 17459363388283345339UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 17459363388283345339UL) + ((uint64_t)op[1] * 13180690705578554535UL) + ((uint64_t)op[2] * 15807971680528568648UL) + ((uint64_t)op[3] * 10679651914351631252UL) + ((uint64_t)op[4] * 12715087827704712706UL) + ((uint64_t)op[5] * 3814736251555248829UL) + ((((uint64_t)op[6] * 7063588192017465263UL) + ((uint64_t)op[7] * 13785887791050718948UL) + ((uint64_t)op[8] * 6489736271999419374UL) + ((uint64_t)op[9] * 7810650901096139443UL) + ((uint64_t)op[10] * 10992916032076318400UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 10992916032076318400UL) + ((uint64_t)op[1] * 17459363388283345339UL) + ((uint64_t)op[2] * 13180690705578554535UL) + ((uint64_t)op[3] * 15807971680528568648UL) + ((uint64_t)op[4] * 10679651914351631252UL) + ((uint64_t)op[5] * 12715087827704712706UL) + ((uint64_t)op[6] * 3814736251555248829UL) + ((((uint64_t)op[7] * 7063588192017465263UL) + ((uint64_t)op[8] * 13785887791050718948UL) + ((uint64_t)op[9] * 6489736271999419374UL) + ((uint64_t)op[10] * 7810650901096139443UL)) * 6);
	tmp_q[7] = ((uint64_t)op[0] * 7810650901096139443UL) + ((uint64_t)op[1] * 10992916032076318400UL) + ((uint64_t)op[2] * 17459363388283345339UL) + ((uint64_t)op[3] * 13180690705578554535UL) + ((uint64_t)op[4] * 15807971680528568648UL) + ((uint64_t)op[5] * 10679651914351631252UL) + ((uint64_t)op[6] * 12715087827704712706UL) + ((uint64_t)op[7] * 3814736251555248829UL) + ((((uint64_t)op[8] * 7063588192017465263UL) + ((uint64_t)op[9] * 13785887791050718948UL) + ((uint64_t)op[10] * 6489736271999419374UL)) * 6);
	tmp_q[8] = ((uint64_t)op[0] * 6489736271999419374UL) + ((uint64_t)op[1] * 7810650901096139443UL) + ((uint64_t)op[2] * 10992916032076318400UL) + ((uint64_t)op[3] * 17459363388283345339UL) + ((uint64_t)op[4] * 13180690705578554535UL) + ((uint64_t)op[5] * 15807971680528568648UL) + ((uint64_t)op[6] * 10679651914351631252UL) + ((uint64_t)op[7] * 12715087827704712706UL) + ((uint64_t)op[8] * 3814736251555248829UL) + ((((uint64_t)op[9] * 7063588192017465263UL) + ((uint64_t)op[10] * 13785887791050718948UL)) * 6);
	tmp_q[9] = ((uint64_t)op[0] * 13785887791050718948UL) + ((uint64_t)op[1] * 6489736271999419374UL) + ((uint64_t)op[2] * 7810650901096139443UL) + ((uint64_t)op[3] * 10992916032076318400UL) + ((uint64_t)op[4] * 17459363388283345339UL) + ((uint64_t)op[5] * 13180690705578554535UL) + ((uint64_t)op[6] * 15807971680528568648UL) + ((uint64_t)op[7] * 10679651914351631252UL) + ((uint64_t)op[8] * 12715087827704712706UL) + ((uint64_t)op[9] * 3814736251555248829UL) + ((uint64_t)op[10] * 5488041004685688346UL);
	tmp_q[10] = ((uint64_t)op[0] * 7063588192017465263UL) + ((uint64_t)op[1] * 13785887791050718948UL) + ((uint64_t)op[2] * 6489736271999419374UL) + ((uint64_t)op[3] * 7810650901096139443UL) + ((uint64_t)op[4] * 10992916032076318400UL) + ((uint64_t)op[5] * 17459363388283345339UL) + ((uint64_t)op[6] * 13180690705578554535UL) + ((uint64_t)op[7] * 15807971680528568648UL) + ((uint64_t)op[8] * 10679651914351631252UL) + ((uint64_t)op[9] * 12715087827704712706UL) + ((uint64_t)op[10] * 3814736251555248829UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 13733484044925L) + ((-((int128)tmp_q[1] * 91313630923954L) - ((int128)tmp_q[2] * 76280152600576L) + ((int128)tmp_q[3] * 44896512025963L) - ((int128)tmp_q[4] * 25217876202649L) + ((int128)tmp_q[5] * 76724082589754L) - ((int128)tmp_q[6] * 31124482196847L) - ((int128)tmp_q[7] * 29626453907169L) + ((int128)tmp_q[8] * 41245934410412L) + ((int128)tmp_q[9] * 19248129186934L) - ((int128)tmp_q[10] * 38224028833948L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 38224028833948L) - ((int128)tmp_q[1] * 13733484044925L) + ((-((int128)tmp_q[2] * 91313630923954L) - ((int128)tmp_q[3] * 76280152600576L) + ((int128)tmp_q[4] * 44896512025963L) - ((int128)tmp_q[5] * 25217876202649L) + ((int128)tmp_q[6] * 76724082589754L) - ((int128)tmp_q[7] * 31124482196847L) - ((int128)tmp_q[8] * 29626453907169L) + ((int128)tmp_q[9] * 41245934410412L) + ((int128)tmp_q[10] * 19248129186934L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 19248129186934L) - ((int128)tmp_q[1] * 38224028833948L) - ((int128)tmp_q[2] * 13733484044925L) + ((-((int128)tmp_q[3] * 91313630923954L) - ((int128)tmp_q[4] * 76280152600576L) + ((int128)tmp_q[5] * 44896512025963L) - ((int128)tmp_q[6] * 25217876202649L) + ((int128)tmp_q[7] * 76724082589754L) - ((int128)tmp_q[8] * 31124482196847L) - ((int128)tmp_q[9] * 29626453907169L) + ((int128)tmp_q[10] * 41245934410412L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 41245934410412L) + ((int128)tmp_q[1] * 19248129186934L) - ((int128)tmp_q[2] * 38224028833948L) - ((int128)tmp_q[3] * 13733484044925L) + ((-((int128)tmp_q[4] * 91313630923954L) - ((int128)tmp_q[5] * 76280152600576L) + ((int128)tmp_q[6] * 44896512025963L) - ((int128)tmp_q[7] * 25217876202649L) + ((int128)tmp_q[8] * 76724082589754L) - ((int128)tmp_q[9] * 31124482196847L) - ((int128)tmp_q[10] * 29626453907169L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 29626453907169L) + ((int128)tmp_q[1] * 41245934410412L) + ((int128)tmp_q[2] * 19248129186934L) - ((int128)tmp_q[3] * 38224028833948L) - ((int128)tmp_q[4] * 13733484044925L) + ((-((int128)tmp_q[5] * 91313630923954L) - ((int128)tmp_q[6] * 76280152600576L) + ((int128)tmp_q[7] * 44896512025963L) - ((int128)tmp_q[8] * 25217876202649L) + ((int128)tmp_q[9] * 76724082589754L) - ((int128)tmp_q[10] * 31124482196847L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 31124482196847L) - ((int128)tmp_q[1] * 29626453907169L) + ((int128)tmp_q[2] * 41245934410412L) + ((int128)tmp_q[3] * 19248129186934L) - ((int128)tmp_q[4] * 38224028833948L) - ((int128)tmp_q[5] * 13733484044925L) + ((-((int128)tmp_q[6] * 91313630923954L) - ((int128)tmp_q[7] * 76280152600576L) + ((int128)tmp_q[8] * 44896512025963L) - ((int128)tmp_q[9] * 25217876202649L) + ((int128)tmp_q[10] * 76724082589754L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 76724082589754L) - ((int128)tmp_q[1] * 31124482196847L) - ((int128)tmp_q[2] * 29626453907169L) + ((int128)tmp_q[3] * 41245934410412L) + ((int128)tmp_q[4] * 19248129186934L) - ((int128)tmp_q[5] * 38224028833948L) - ((int128)tmp_q[6] * 13733484044925L) + ((-((int128)tmp_q[7] * 91313630923954L) - ((int128)tmp_q[8] * 76280152600576L) + ((int128)tmp_q[9] * 44896512025963L) - ((int128)tmp_q[10] * 25217876202649L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 25217876202649L) + ((int128)tmp_q[1] * 76724082589754L) - ((int128)tmp_q[2] * 31124482196847L) - ((int128)tmp_q[3] * 29626453907169L) + ((int128)tmp_q[4] * 41245934410412L) + ((int128)tmp_q[5] * 19248129186934L) - ((int128)tmp_q[6] * 38224028833948L) - ((int128)tmp_q[7] * 13733484044925L) + ((-((int128)tmp_q[8] * 91313630923954L) - ((int128)tmp_q[9] * 76280152600576L) + ((int128)tmp_q[10] * 44896512025963L)) * 6);
	tmp_zero[8] = ((int128)tmp_q[0] * 44896512025963L) - ((int128)tmp_q[1] * 25217876202649L) + ((int128)tmp_q[2] * 76724082589754L) - ((int128)tmp_q[3] * 31124482196847L) - ((int128)tmp_q[4] * 29626453907169L) + ((int128)tmp_q[5] * 41245934410412L) + ((int128)tmp_q[6] * 19248129186934L) - ((int128)tmp_q[7] * 38224028833948L) - ((int128)tmp_q[8] * 13733484044925L) + ((-((int128)tmp_q[9] * 91313630923954L) - ((int128)tmp_q[10] * 76280152600576L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 76280152600576L) + ((int128)tmp_q[1] * 44896512025963L) - ((int128)tmp_q[2] * 25217876202649L) + ((int128)tmp_q[3] * 76724082589754L) - ((int128)tmp_q[4] * 31124482196847L) - ((int128)tmp_q[5] * 29626453907169L) + ((int128)tmp_q[6] * 41245934410412L) + ((int128)tmp_q[7] * 19248129186934L) - ((int128)tmp_q[8] * 38224028833948L) - ((int128)tmp_q[9] * 13733484044925L) - ((int128)tmp_q[10] * 547881785543724L);
	tmp_zero[10] = -((int128)tmp_q[0] * 91313630923954L) - ((int128)tmp_q[1] * 76280152600576L) + ((int128)tmp_q[2] * 44896512025963L) - ((int128)tmp_q[3] * 25217876202649L) + ((int128)tmp_q[4] * 76724082589754L) - ((int128)tmp_q[5] * 31124482196847L) - ((int128)tmp_q[6] * 29626453907169L) + ((int128)tmp_q[7] * 41245934410412L) + ((int128)tmp_q[8] * 19248129186934L) - ((int128)tmp_q[9] * 38224028833948L) - ((int128)tmp_q[10] * 13733484044925L);

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

