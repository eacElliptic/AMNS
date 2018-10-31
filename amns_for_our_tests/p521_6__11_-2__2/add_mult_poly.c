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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14270945964062237637UL) + ((((uint64_t)op[1] * 2669329516822033678UL) + ((uint64_t)op[2] * 12018502344290750430UL) + ((uint64_t)op[3] * 13568889618045874363UL) + ((uint64_t)op[4] * 12180621306793228983UL) + ((uint64_t)op[5] * 11363757983100489948UL) + ((uint64_t)op[6] * 233692023464269553UL) + ((uint64_t)op[7] * 3119192727961219111UL) + ((uint64_t)op[8] * 7243046096440637249UL) + ((uint64_t)op[9] * 11861273293012357964UL) + ((uint64_t)op[10] * 10768416583599282701UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 10768416583599282701UL) + ((uint64_t)op[1] * 14270945964062237637UL) + ((((uint64_t)op[2] * 2669329516822033678UL) + ((uint64_t)op[3] * 12018502344290750430UL) + ((uint64_t)op[4] * 13568889618045874363UL) + ((uint64_t)op[5] * 12180621306793228983UL) + ((uint64_t)op[6] * 11363757983100489948UL) + ((uint64_t)op[7] * 233692023464269553UL) + ((uint64_t)op[8] * 3119192727961219111UL) + ((uint64_t)op[9] * 7243046096440637249UL) + ((uint64_t)op[10] * 11861273293012357964UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 11861273293012357964UL) + ((uint64_t)op[1] * 10768416583599282701UL) + ((uint64_t)op[2] * 14270945964062237637UL) + ((((uint64_t)op[3] * 2669329516822033678UL) + ((uint64_t)op[4] * 12018502344290750430UL) + ((uint64_t)op[5] * 13568889618045874363UL) + ((uint64_t)op[6] * 12180621306793228983UL) + ((uint64_t)op[7] * 11363757983100489948UL) + ((uint64_t)op[8] * 233692023464269553UL) + ((uint64_t)op[9] * 3119192727961219111UL) + ((uint64_t)op[10] * 7243046096440637249UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 7243046096440637249UL) + ((uint64_t)op[1] * 11861273293012357964UL) + ((uint64_t)op[2] * 10768416583599282701UL) + ((uint64_t)op[3] * 14270945964062237637UL) + ((((uint64_t)op[4] * 2669329516822033678UL) + ((uint64_t)op[5] * 12018502344290750430UL) + ((uint64_t)op[6] * 13568889618045874363UL) + ((uint64_t)op[7] * 12180621306793228983UL) + ((uint64_t)op[8] * 11363757983100489948UL) + ((uint64_t)op[9] * 233692023464269553UL) + ((uint64_t)op[10] * 3119192727961219111UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 3119192727961219111UL) + ((uint64_t)op[1] * 7243046096440637249UL) + ((uint64_t)op[2] * 11861273293012357964UL) + ((uint64_t)op[3] * 10768416583599282701UL) + ((uint64_t)op[4] * 14270945964062237637UL) + ((((uint64_t)op[5] * 2669329516822033678UL) + ((uint64_t)op[6] * 12018502344290750430UL) + ((uint64_t)op[7] * 13568889618045874363UL) + ((uint64_t)op[8] * 12180621306793228983UL) + ((uint64_t)op[9] * 11363757983100489948UL) + ((uint64_t)op[10] * 233692023464269553UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 233692023464269553UL) + ((uint64_t)op[1] * 3119192727961219111UL) + ((uint64_t)op[2] * 7243046096440637249UL) + ((uint64_t)op[3] * 11861273293012357964UL) + ((uint64_t)op[4] * 10768416583599282701UL) + ((uint64_t)op[5] * 14270945964062237637UL) + ((((uint64_t)op[6] * 2669329516822033678UL) + ((uint64_t)op[7] * 12018502344290750430UL) + ((uint64_t)op[8] * 13568889618045874363UL) + ((uint64_t)op[9] * 12180621306793228983UL) + ((uint64_t)op[10] * 11363757983100489948UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 11363757983100489948UL) + ((uint64_t)op[1] * 233692023464269553UL) + ((uint64_t)op[2] * 3119192727961219111UL) + ((uint64_t)op[3] * 7243046096440637249UL) + ((uint64_t)op[4] * 11861273293012357964UL) + ((uint64_t)op[5] * 10768416583599282701UL) + ((uint64_t)op[6] * 14270945964062237637UL) + ((((uint64_t)op[7] * 2669329516822033678UL) + ((uint64_t)op[8] * 12018502344290750430UL) + ((uint64_t)op[9] * 13568889618045874363UL) + ((uint64_t)op[10] * 12180621306793228983UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 12180621306793228983UL) + ((uint64_t)op[1] * 11363757983100489948UL) + ((uint64_t)op[2] * 233692023464269553UL) + ((uint64_t)op[3] * 3119192727961219111UL) + ((uint64_t)op[4] * 7243046096440637249UL) + ((uint64_t)op[5] * 11861273293012357964UL) + ((uint64_t)op[6] * 10768416583599282701UL) + ((uint64_t)op[7] * 14270945964062237637UL) + ((((uint64_t)op[8] * 2669329516822033678UL) + ((uint64_t)op[9] * 12018502344290750430UL) + ((uint64_t)op[10] * 13568889618045874363UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 13568889618045874363UL) + ((uint64_t)op[1] * 12180621306793228983UL) + ((uint64_t)op[2] * 11363757983100489948UL) + ((uint64_t)op[3] * 233692023464269553UL) + ((uint64_t)op[4] * 3119192727961219111UL) + ((uint64_t)op[5] * 7243046096440637249UL) + ((uint64_t)op[6] * 11861273293012357964UL) + ((uint64_t)op[7] * 10768416583599282701UL) + ((uint64_t)op[8] * 14270945964062237637UL) + ((((uint64_t)op[9] * 2669329516822033678UL) + ((uint64_t)op[10] * 12018502344290750430UL)) * 18446744073709551614);
	tmp_q[9] = ((uint64_t)op[0] * 12018502344290750430UL) + ((uint64_t)op[1] * 13568889618045874363UL) + ((uint64_t)op[2] * 12180621306793228983UL) + ((uint64_t)op[3] * 11363757983100489948UL) + ((uint64_t)op[4] * 233692023464269553UL) + ((uint64_t)op[5] * 3119192727961219111UL) + ((uint64_t)op[6] * 7243046096440637249UL) + ((uint64_t)op[7] * 11861273293012357964UL) + ((uint64_t)op[8] * 10768416583599282701UL) + ((uint64_t)op[9] * 14270945964062237637UL) + ((uint64_t)op[10] * 13108085040065484260UL);
	tmp_q[10] = ((uint64_t)op[0] * 2669329516822033678UL) + ((uint64_t)op[1] * 12018502344290750430UL) + ((uint64_t)op[2] * 13568889618045874363UL) + ((uint64_t)op[3] * 12180621306793228983UL) + ((uint64_t)op[4] * 11363757983100489948UL) + ((uint64_t)op[5] * 233692023464269553UL) + ((uint64_t)op[6] * 3119192727961219111UL) + ((uint64_t)op[7] * 7243046096440637249UL) + ((uint64_t)op[8] * 11861273293012357964UL) + ((uint64_t)op[9] * 10768416583599282701UL) + ((uint64_t)op[10] * 14270945964062237637UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 79077553721769L) - ((-((int128)tmp_q[1] * 73413668535250L) + ((int128)tmp_q[2] * 92693321186596L) + ((int128)tmp_q[3] * 12728146028306L) + ((int128)tmp_q[4] * 81184355768469L) + ((int128)tmp_q[5] * 37101720706767L) - ((int128)tmp_q[6] * 66358004134435L) - ((int128)tmp_q[7] * 60689900466982L) - ((int128)tmp_q[8] * 99980209118254L) - ((int128)tmp_q[9] * 61945083820369L) - ((int128)tmp_q[10] * 19742727880395L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 19742727880395L) - ((int128)tmp_q[1] * 79077553721769L) - ((-((int128)tmp_q[2] * 73413668535250L) + ((int128)tmp_q[3] * 92693321186596L) + ((int128)tmp_q[4] * 12728146028306L) + ((int128)tmp_q[5] * 81184355768469L) + ((int128)tmp_q[6] * 37101720706767L) - ((int128)tmp_q[7] * 66358004134435L) - ((int128)tmp_q[8] * 60689900466982L) - ((int128)tmp_q[9] * 99980209118254L) - ((int128)tmp_q[10] * 61945083820369L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 61945083820369L) - ((int128)tmp_q[1] * 19742727880395L) - ((int128)tmp_q[2] * 79077553721769L) - ((-((int128)tmp_q[3] * 73413668535250L) + ((int128)tmp_q[4] * 92693321186596L) + ((int128)tmp_q[5] * 12728146028306L) + ((int128)tmp_q[6] * 81184355768469L) + ((int128)tmp_q[7] * 37101720706767L) - ((int128)tmp_q[8] * 66358004134435L) - ((int128)tmp_q[9] * 60689900466982L) - ((int128)tmp_q[10] * 99980209118254L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 99980209118254L) - ((int128)tmp_q[1] * 61945083820369L) - ((int128)tmp_q[2] * 19742727880395L) - ((int128)tmp_q[3] * 79077553721769L) - ((-((int128)tmp_q[4] * 73413668535250L) + ((int128)tmp_q[5] * 92693321186596L) + ((int128)tmp_q[6] * 12728146028306L) + ((int128)tmp_q[7] * 81184355768469L) + ((int128)tmp_q[8] * 37101720706767L) - ((int128)tmp_q[9] * 66358004134435L) - ((int128)tmp_q[10] * 60689900466982L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 60689900466982L) - ((int128)tmp_q[1] * 99980209118254L) - ((int128)tmp_q[2] * 61945083820369L) - ((int128)tmp_q[3] * 19742727880395L) - ((int128)tmp_q[4] * 79077553721769L) - ((-((int128)tmp_q[5] * 73413668535250L) + ((int128)tmp_q[6] * 92693321186596L) + ((int128)tmp_q[7] * 12728146028306L) + ((int128)tmp_q[8] * 81184355768469L) + ((int128)tmp_q[9] * 37101720706767L) - ((int128)tmp_q[10] * 66358004134435L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 66358004134435L) - ((int128)tmp_q[1] * 60689900466982L) - ((int128)tmp_q[2] * 99980209118254L) - ((int128)tmp_q[3] * 61945083820369L) - ((int128)tmp_q[4] * 19742727880395L) - ((int128)tmp_q[5] * 79077553721769L) - ((-((int128)tmp_q[6] * 73413668535250L) + ((int128)tmp_q[7] * 92693321186596L) + ((int128)tmp_q[8] * 12728146028306L) + ((int128)tmp_q[9] * 81184355768469L) + ((int128)tmp_q[10] * 37101720706767L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 37101720706767L) - ((int128)tmp_q[1] * 66358004134435L) - ((int128)tmp_q[2] * 60689900466982L) - ((int128)tmp_q[3] * 99980209118254L) - ((int128)tmp_q[4] * 61945083820369L) - ((int128)tmp_q[5] * 19742727880395L) - ((int128)tmp_q[6] * 79077553721769L) - ((-((int128)tmp_q[7] * 73413668535250L) + ((int128)tmp_q[8] * 92693321186596L) + ((int128)tmp_q[9] * 12728146028306L) + ((int128)tmp_q[10] * 81184355768469L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 81184355768469L) + ((int128)tmp_q[1] * 37101720706767L) - ((int128)tmp_q[2] * 66358004134435L) - ((int128)tmp_q[3] * 60689900466982L) - ((int128)tmp_q[4] * 99980209118254L) - ((int128)tmp_q[5] * 61945083820369L) - ((int128)tmp_q[6] * 19742727880395L) - ((int128)tmp_q[7] * 79077553721769L) - ((-((int128)tmp_q[8] * 73413668535250L) + ((int128)tmp_q[9] * 92693321186596L) + ((int128)tmp_q[10] * 12728146028306L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 12728146028306L) + ((int128)tmp_q[1] * 81184355768469L) + ((int128)tmp_q[2] * 37101720706767L) - ((int128)tmp_q[3] * 66358004134435L) - ((int128)tmp_q[4] * 60689900466982L) - ((int128)tmp_q[5] * 99980209118254L) - ((int128)tmp_q[6] * 61945083820369L) - ((int128)tmp_q[7] * 19742727880395L) - ((int128)tmp_q[8] * 79077553721769L) - ((-((int128)tmp_q[9] * 73413668535250L) + ((int128)tmp_q[10] * 92693321186596L)) * 2);
	tmp_zero[9] = ((int128)tmp_q[0] * 92693321186596L) + ((int128)tmp_q[1] * 12728146028306L) + ((int128)tmp_q[2] * 81184355768469L) + ((int128)tmp_q[3] * 37101720706767L) - ((int128)tmp_q[4] * 66358004134435L) - ((int128)tmp_q[5] * 60689900466982L) - ((int128)tmp_q[6] * 99980209118254L) - ((int128)tmp_q[7] * 61945083820369L) - ((int128)tmp_q[8] * 19742727880395L) - ((int128)tmp_q[9] * 79077553721769L) + ((int128)tmp_q[10] * 146827337070500L);
	tmp_zero[10] = -((int128)tmp_q[0] * 73413668535250L) + ((int128)tmp_q[1] * 92693321186596L) + ((int128)tmp_q[2] * 12728146028306L) + ((int128)tmp_q[3] * 81184355768469L) + ((int128)tmp_q[4] * 37101720706767L) - ((int128)tmp_q[5] * 66358004134435L) - ((int128)tmp_q[6] * 60689900466982L) - ((int128)tmp_q[7] * 99980209118254L) - ((int128)tmp_q[8] * 61945083820369L) - ((int128)tmp_q[9] * 19742727880395L) - ((int128)tmp_q[10] * 79077553721769L);

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

