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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 3);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 3);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7707991230282528475UL) + ((((uint64_t)op[1] * 13423033651596950334UL) + ((uint64_t)op[2] * 14778165546131683497UL) + ((uint64_t)op[3] * 6351247225915032413UL) + ((uint64_t)op[4] * 1001678668802483749UL) + ((uint64_t)op[5] * 6661174754742074390UL) + ((uint64_t)op[6] * 1141185034377407385UL) + ((uint64_t)op[7] * 4185578370872289756UL) + ((uint64_t)op[8] * 14019192090982233556UL) + ((uint64_t)op[9] * 7766537469366477796UL) + ((uint64_t)op[10] * 3170690343034795015UL) + ((uint64_t)op[11] * 6040873488720305521UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 6040873488720305521UL) + ((uint64_t)op[1] * 7707991230282528475UL) + ((((uint64_t)op[2] * 13423033651596950334UL) + ((uint64_t)op[3] * 14778165546131683497UL) + ((uint64_t)op[4] * 6351247225915032413UL) + ((uint64_t)op[5] * 1001678668802483749UL) + ((uint64_t)op[6] * 6661174754742074390UL) + ((uint64_t)op[7] * 1141185034377407385UL) + ((uint64_t)op[8] * 4185578370872289756UL) + ((uint64_t)op[9] * 14019192090982233556UL) + ((uint64_t)op[10] * 7766537469366477796UL) + ((uint64_t)op[11] * 3170690343034795015UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 3170690343034795015UL) + ((uint64_t)op[1] * 6040873488720305521UL) + ((uint64_t)op[2] * 7707991230282528475UL) + ((((uint64_t)op[3] * 13423033651596950334UL) + ((uint64_t)op[4] * 14778165546131683497UL) + ((uint64_t)op[5] * 6351247225915032413UL) + ((uint64_t)op[6] * 1001678668802483749UL) + ((uint64_t)op[7] * 6661174754742074390UL) + ((uint64_t)op[8] * 1141185034377407385UL) + ((uint64_t)op[9] * 4185578370872289756UL) + ((uint64_t)op[10] * 14019192090982233556UL) + ((uint64_t)op[11] * 7766537469366477796UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7766537469366477796UL) + ((uint64_t)op[1] * 3170690343034795015UL) + ((uint64_t)op[2] * 6040873488720305521UL) + ((uint64_t)op[3] * 7707991230282528475UL) + ((((uint64_t)op[4] * 13423033651596950334UL) + ((uint64_t)op[5] * 14778165546131683497UL) + ((uint64_t)op[6] * 6351247225915032413UL) + ((uint64_t)op[7] * 1001678668802483749UL) + ((uint64_t)op[8] * 6661174754742074390UL) + ((uint64_t)op[9] * 1141185034377407385UL) + ((uint64_t)op[10] * 4185578370872289756UL) + ((uint64_t)op[11] * 14019192090982233556UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 14019192090982233556UL) + ((uint64_t)op[1] * 7766537469366477796UL) + ((uint64_t)op[2] * 3170690343034795015UL) + ((uint64_t)op[3] * 6040873488720305521UL) + ((uint64_t)op[4] * 7707991230282528475UL) + ((((uint64_t)op[5] * 13423033651596950334UL) + ((uint64_t)op[6] * 14778165546131683497UL) + ((uint64_t)op[7] * 6351247225915032413UL) + ((uint64_t)op[8] * 1001678668802483749UL) + ((uint64_t)op[9] * 6661174754742074390UL) + ((uint64_t)op[10] * 1141185034377407385UL) + ((uint64_t)op[11] * 4185578370872289756UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 4185578370872289756UL) + ((uint64_t)op[1] * 14019192090982233556UL) + ((uint64_t)op[2] * 7766537469366477796UL) + ((uint64_t)op[3] * 3170690343034795015UL) + ((uint64_t)op[4] * 6040873488720305521UL) + ((uint64_t)op[5] * 7707991230282528475UL) + ((((uint64_t)op[6] * 13423033651596950334UL) + ((uint64_t)op[7] * 14778165546131683497UL) + ((uint64_t)op[8] * 6351247225915032413UL) + ((uint64_t)op[9] * 1001678668802483749UL) + ((uint64_t)op[10] * 6661174754742074390UL) + ((uint64_t)op[11] * 1141185034377407385UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 1141185034377407385UL) + ((uint64_t)op[1] * 4185578370872289756UL) + ((uint64_t)op[2] * 14019192090982233556UL) + ((uint64_t)op[3] * 7766537469366477796UL) + ((uint64_t)op[4] * 3170690343034795015UL) + ((uint64_t)op[5] * 6040873488720305521UL) + ((uint64_t)op[6] * 7707991230282528475UL) + ((((uint64_t)op[7] * 13423033651596950334UL) + ((uint64_t)op[8] * 14778165546131683497UL) + ((uint64_t)op[9] * 6351247225915032413UL) + ((uint64_t)op[10] * 1001678668802483749UL) + ((uint64_t)op[11] * 6661174754742074390UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 6661174754742074390UL) + ((uint64_t)op[1] * 1141185034377407385UL) + ((uint64_t)op[2] * 4185578370872289756UL) + ((uint64_t)op[3] * 14019192090982233556UL) + ((uint64_t)op[4] * 7766537469366477796UL) + ((uint64_t)op[5] * 3170690343034795015UL) + ((uint64_t)op[6] * 6040873488720305521UL) + ((uint64_t)op[7] * 7707991230282528475UL) + ((((uint64_t)op[8] * 13423033651596950334UL) + ((uint64_t)op[9] * 14778165546131683497UL) + ((uint64_t)op[10] * 6351247225915032413UL) + ((uint64_t)op[11] * 1001678668802483749UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 1001678668802483749UL) + ((uint64_t)op[1] * 6661174754742074390UL) + ((uint64_t)op[2] * 1141185034377407385UL) + ((uint64_t)op[3] * 4185578370872289756UL) + ((uint64_t)op[4] * 14019192090982233556UL) + ((uint64_t)op[5] * 7766537469366477796UL) + ((uint64_t)op[6] * 3170690343034795015UL) + ((uint64_t)op[7] * 6040873488720305521UL) + ((uint64_t)op[8] * 7707991230282528475UL) + ((((uint64_t)op[9] * 13423033651596950334UL) + ((uint64_t)op[10] * 14778165546131683497UL) + ((uint64_t)op[11] * 6351247225915032413UL)) * 18446744073709551613);
	tmp_q[9] = ((uint64_t)op[0] * 6351247225915032413UL) + ((uint64_t)op[1] * 1001678668802483749UL) + ((uint64_t)op[2] * 6661174754742074390UL) + ((uint64_t)op[3] * 1141185034377407385UL) + ((uint64_t)op[4] * 4185578370872289756UL) + ((uint64_t)op[5] * 14019192090982233556UL) + ((uint64_t)op[6] * 7766537469366477796UL) + ((uint64_t)op[7] * 3170690343034795015UL) + ((uint64_t)op[8] * 6040873488720305521UL) + ((uint64_t)op[9] * 7707991230282528475UL) + ((((uint64_t)op[10] * 13423033651596950334UL) + ((uint64_t)op[11] * 14778165546131683497UL)) * 18446744073709551613);
	tmp_q[10] = ((uint64_t)op[0] * 14778165546131683497UL) + ((uint64_t)op[1] * 6351247225915032413UL) + ((uint64_t)op[2] * 1001678668802483749UL) + ((uint64_t)op[3] * 6661174754742074390UL) + ((uint64_t)op[4] * 1141185034377407385UL) + ((uint64_t)op[5] * 4185578370872289756UL) + ((uint64_t)op[6] * 14019192090982233556UL) + ((uint64_t)op[7] * 7766537469366477796UL) + ((uint64_t)op[8] * 3170690343034795015UL) + ((uint64_t)op[9] * 6040873488720305521UL) + ((uint64_t)op[10] * 7707991230282528475UL) + ((uint64_t)op[11] * 15071131266337803846UL);
	tmp_q[11] = ((uint64_t)op[0] * 13423033651596950334UL) + ((uint64_t)op[1] * 14778165546131683497UL) + ((uint64_t)op[2] * 6351247225915032413UL) + ((uint64_t)op[3] * 1001678668802483749UL) + ((uint64_t)op[4] * 6661174754742074390UL) + ((uint64_t)op[5] * 1141185034377407385UL) + ((uint64_t)op[6] * 4185578370872289756UL) + ((uint64_t)op[7] * 14019192090982233556UL) + ((uint64_t)op[8] * 7766537469366477796UL) + ((uint64_t)op[9] * 3170690343034795015UL) + ((uint64_t)op[10] * 6040873488720305521UL) + ((uint64_t)op[11] * 7707991230282528475UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2924914894690L) - ((-((int128)tmp_q[1] * 202435862179L) - ((int128)tmp_q[2] * 3370810424012L) + ((int128)tmp_q[3] * 788793267333L) + ((int128)tmp_q[4] * 2450650351909L) + ((int128)tmp_q[5] * 1794048204769L) - ((int128)tmp_q[6] * 5917939366890L) - ((int128)tmp_q[7] * 5935450709161L) + ((int128)tmp_q[8] * 4843110636457L) + ((int128)tmp_q[9] * 3362464491198L) - ((int128)tmp_q[10] * 5221456108027L) + ((int128)tmp_q[11] * 6768668106762L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 6768668106762L) - ((int128)tmp_q[1] * 2924914894690L) - ((-((int128)tmp_q[2] * 202435862179L) - ((int128)tmp_q[3] * 3370810424012L) + ((int128)tmp_q[4] * 788793267333L) + ((int128)tmp_q[5] * 2450650351909L) + ((int128)tmp_q[6] * 1794048204769L) - ((int128)tmp_q[7] * 5917939366890L) - ((int128)tmp_q[8] * 5935450709161L) + ((int128)tmp_q[9] * 4843110636457L) + ((int128)tmp_q[10] * 3362464491198L) - ((int128)tmp_q[11] * 5221456108027L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 5221456108027L) + ((int128)tmp_q[1] * 6768668106762L) - ((int128)tmp_q[2] * 2924914894690L) - ((-((int128)tmp_q[3] * 202435862179L) - ((int128)tmp_q[4] * 3370810424012L) + ((int128)tmp_q[5] * 788793267333L) + ((int128)tmp_q[6] * 2450650351909L) + ((int128)tmp_q[7] * 1794048204769L) - ((int128)tmp_q[8] * 5917939366890L) - ((int128)tmp_q[9] * 5935450709161L) + ((int128)tmp_q[10] * 4843110636457L) + ((int128)tmp_q[11] * 3362464491198L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 3362464491198L) - ((int128)tmp_q[1] * 5221456108027L) + ((int128)tmp_q[2] * 6768668106762L) - ((int128)tmp_q[3] * 2924914894690L) - ((-((int128)tmp_q[4] * 202435862179L) - ((int128)tmp_q[5] * 3370810424012L) + ((int128)tmp_q[6] * 788793267333L) + ((int128)tmp_q[7] * 2450650351909L) + ((int128)tmp_q[8] * 1794048204769L) - ((int128)tmp_q[9] * 5917939366890L) - ((int128)tmp_q[10] * 5935450709161L) + ((int128)tmp_q[11] * 4843110636457L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 4843110636457L) + ((int128)tmp_q[1] * 3362464491198L) - ((int128)tmp_q[2] * 5221456108027L) + ((int128)tmp_q[3] * 6768668106762L) - ((int128)tmp_q[4] * 2924914894690L) - ((-((int128)tmp_q[5] * 202435862179L) - ((int128)tmp_q[6] * 3370810424012L) + ((int128)tmp_q[7] * 788793267333L) + ((int128)tmp_q[8] * 2450650351909L) + ((int128)tmp_q[9] * 1794048204769L) - ((int128)tmp_q[10] * 5917939366890L) - ((int128)tmp_q[11] * 5935450709161L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 5935450709161L) + ((int128)tmp_q[1] * 4843110636457L) + ((int128)tmp_q[2] * 3362464491198L) - ((int128)tmp_q[3] * 5221456108027L) + ((int128)tmp_q[4] * 6768668106762L) - ((int128)tmp_q[5] * 2924914894690L) - ((-((int128)tmp_q[6] * 202435862179L) - ((int128)tmp_q[7] * 3370810424012L) + ((int128)tmp_q[8] * 788793267333L) + ((int128)tmp_q[9] * 2450650351909L) + ((int128)tmp_q[10] * 1794048204769L) - ((int128)tmp_q[11] * 5917939366890L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 5917939366890L) - ((int128)tmp_q[1] * 5935450709161L) + ((int128)tmp_q[2] * 4843110636457L) + ((int128)tmp_q[3] * 3362464491198L) - ((int128)tmp_q[4] * 5221456108027L) + ((int128)tmp_q[5] * 6768668106762L) - ((int128)tmp_q[6] * 2924914894690L) - ((-((int128)tmp_q[7] * 202435862179L) - ((int128)tmp_q[8] * 3370810424012L) + ((int128)tmp_q[9] * 788793267333L) + ((int128)tmp_q[10] * 2450650351909L) + ((int128)tmp_q[11] * 1794048204769L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 1794048204769L) - ((int128)tmp_q[1] * 5917939366890L) - ((int128)tmp_q[2] * 5935450709161L) + ((int128)tmp_q[3] * 4843110636457L) + ((int128)tmp_q[4] * 3362464491198L) - ((int128)tmp_q[5] * 5221456108027L) + ((int128)tmp_q[6] * 6768668106762L) - ((int128)tmp_q[7] * 2924914894690L) - ((-((int128)tmp_q[8] * 202435862179L) - ((int128)tmp_q[9] * 3370810424012L) + ((int128)tmp_q[10] * 788793267333L) + ((int128)tmp_q[11] * 2450650351909L)) * 3);
	tmp_zero[8] = ((int128)tmp_q[0] * 2450650351909L) + ((int128)tmp_q[1] * 1794048204769L) - ((int128)tmp_q[2] * 5917939366890L) - ((int128)tmp_q[3] * 5935450709161L) + ((int128)tmp_q[4] * 4843110636457L) + ((int128)tmp_q[5] * 3362464491198L) - ((int128)tmp_q[6] * 5221456108027L) + ((int128)tmp_q[7] * 6768668106762L) - ((int128)tmp_q[8] * 2924914894690L) - ((-((int128)tmp_q[9] * 202435862179L) - ((int128)tmp_q[10] * 3370810424012L) + ((int128)tmp_q[11] * 788793267333L)) * 3);
	tmp_zero[9] = ((int128)tmp_q[0] * 788793267333L) + ((int128)tmp_q[1] * 2450650351909L) + ((int128)tmp_q[2] * 1794048204769L) - ((int128)tmp_q[3] * 5917939366890L) - ((int128)tmp_q[4] * 5935450709161L) + ((int128)tmp_q[5] * 4843110636457L) + ((int128)tmp_q[6] * 3362464491198L) - ((int128)tmp_q[7] * 5221456108027L) + ((int128)tmp_q[8] * 6768668106762L) - ((int128)tmp_q[9] * 2924914894690L) - ((-((int128)tmp_q[10] * 202435862179L) - ((int128)tmp_q[11] * 3370810424012L)) * 3);
	tmp_zero[10] = -((int128)tmp_q[0] * 3370810424012L) + ((int128)tmp_q[1] * 788793267333L) + ((int128)tmp_q[2] * 2450650351909L) + ((int128)tmp_q[3] * 1794048204769L) - ((int128)tmp_q[4] * 5917939366890L) - ((int128)tmp_q[5] * 5935450709161L) + ((int128)tmp_q[6] * 4843110636457L) + ((int128)tmp_q[7] * 3362464491198L) - ((int128)tmp_q[8] * 5221456108027L) + ((int128)tmp_q[9] * 6768668106762L) - ((int128)tmp_q[10] * 2924914894690L) + ((int128)tmp_q[11] * 607307586537L);
	tmp_zero[11] = -((int128)tmp_q[0] * 202435862179L) - ((int128)tmp_q[1] * 3370810424012L) + ((int128)tmp_q[2] * 788793267333L) + ((int128)tmp_q[3] * 2450650351909L) + ((int128)tmp_q[4] * 1794048204769L) - ((int128)tmp_q[5] * 5917939366890L) - ((int128)tmp_q[6] * 5935450709161L) + ((int128)tmp_q[7] * 4843110636457L) + ((int128)tmp_q[8] * 3362464491198L) - ((int128)tmp_q[9] * 5221456108027L) + ((int128)tmp_q[10] * 6768668106762L) - ((int128)tmp_q[11] * 2924914894690L);

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

