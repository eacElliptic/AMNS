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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4072912300693416333UL) + ((((uint64_t)op[1] * 887849535823447690UL) + ((uint64_t)op[2] * 11505017989523014452UL) + ((uint64_t)op[3] * 17086279776000340079UL) + ((uint64_t)op[4] * 12094334716622833041UL) + ((uint64_t)op[5] * 7515338909455392165UL) + ((uint64_t)op[6] * 13019615329185958376UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 13019615329185958376UL) + ((uint64_t)op[1] * 4072912300693416333UL) + ((((uint64_t)op[2] * 887849535823447690UL) + ((uint64_t)op[3] * 11505017989523014452UL) + ((uint64_t)op[4] * 17086279776000340079UL) + ((uint64_t)op[5] * 12094334716622833041UL) + ((uint64_t)op[6] * 7515338909455392165UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 7515338909455392165UL) + ((uint64_t)op[1] * 13019615329185958376UL) + ((uint64_t)op[2] * 4072912300693416333UL) + ((((uint64_t)op[3] * 887849535823447690UL) + ((uint64_t)op[4] * 11505017989523014452UL) + ((uint64_t)op[5] * 17086279776000340079UL) + ((uint64_t)op[6] * 12094334716622833041UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 12094334716622833041UL) + ((uint64_t)op[1] * 7515338909455392165UL) + ((uint64_t)op[2] * 13019615329185958376UL) + ((uint64_t)op[3] * 4072912300693416333UL) + ((((uint64_t)op[4] * 887849535823447690UL) + ((uint64_t)op[5] * 11505017989523014452UL) + ((uint64_t)op[6] * 17086279776000340079UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 17086279776000340079UL) + ((uint64_t)op[1] * 12094334716622833041UL) + ((uint64_t)op[2] * 7515338909455392165UL) + ((uint64_t)op[3] * 13019615329185958376UL) + ((uint64_t)op[4] * 4072912300693416333UL) + ((((uint64_t)op[5] * 887849535823447690UL) + ((uint64_t)op[6] * 11505017989523014452UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 11505017989523014452UL) + ((uint64_t)op[1] * 17086279776000340079UL) + ((uint64_t)op[2] * 12094334716622833041UL) + ((uint64_t)op[3] * 7515338909455392165UL) + ((uint64_t)op[4] * 13019615329185958376UL) + ((uint64_t)op[5] * 4072912300693416333UL) + ((uint64_t)op[6] * 16671045002062656236UL);
	tmp_q[6] = ((uint64_t)op[0] * 887849535823447690UL) + ((uint64_t)op[1] * 11505017989523014452UL) + ((uint64_t)op[2] * 17086279776000340079UL) + ((uint64_t)op[3] * 12094334716622833041UL) + ((uint64_t)op[4] * 7515338909455392165UL) + ((uint64_t)op[5] * 13019615329185958376UL) + ((uint64_t)op[6] * 4072912300693416333UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 7177305999237621L) - ((-((int128)tmp_q[1] * 2547009489666076L) - ((int128)tmp_q[2] * 17865456038846762L) + ((int128)tmp_q[3] * 7932623374667850L) - ((int128)tmp_q[4] * 135105747378269L) + ((int128)tmp_q[5] * 12827490240307099L) - ((int128)tmp_q[6] * 2883751709374212L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 2883751709374212L) + ((int128)tmp_q[1] * 7177305999237621L) - ((-((int128)tmp_q[2] * 2547009489666076L) - ((int128)tmp_q[3] * 17865456038846762L) + ((int128)tmp_q[4] * 7932623374667850L) - ((int128)tmp_q[5] * 135105747378269L) + ((int128)tmp_q[6] * 12827490240307099L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 12827490240307099L) - ((int128)tmp_q[1] * 2883751709374212L) + ((int128)tmp_q[2] * 7177305999237621L) - ((-((int128)tmp_q[3] * 2547009489666076L) - ((int128)tmp_q[4] * 17865456038846762L) + ((int128)tmp_q[5] * 7932623374667850L) - ((int128)tmp_q[6] * 135105747378269L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 135105747378269L) + ((int128)tmp_q[1] * 12827490240307099L) - ((int128)tmp_q[2] * 2883751709374212L) + ((int128)tmp_q[3] * 7177305999237621L) - ((-((int128)tmp_q[4] * 2547009489666076L) - ((int128)tmp_q[5] * 17865456038846762L) + ((int128)tmp_q[6] * 7932623374667850L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 7932623374667850L) - ((int128)tmp_q[1] * 135105747378269L) + ((int128)tmp_q[2] * 12827490240307099L) - ((int128)tmp_q[3] * 2883751709374212L) + ((int128)tmp_q[4] * 7177305999237621L) - ((-((int128)tmp_q[5] * 2547009489666076L) - ((int128)tmp_q[6] * 17865456038846762L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 17865456038846762L) + ((int128)tmp_q[1] * 7932623374667850L) - ((int128)tmp_q[2] * 135105747378269L) + ((int128)tmp_q[3] * 12827490240307099L) - ((int128)tmp_q[4] * 2883751709374212L) + ((int128)tmp_q[5] * 7177305999237621L) + ((int128)tmp_q[6] * 5094018979332152L);
	tmp_zero[6] = -((int128)tmp_q[0] * 2547009489666076L) - ((int128)tmp_q[1] * 17865456038846762L) + ((int128)tmp_q[2] * 7932623374667850L) - ((int128)tmp_q[3] * 135105747378269L) + ((int128)tmp_q[4] * 12827490240307099L) - ((int128)tmp_q[5] * 2883751709374212L) + ((int128)tmp_q[6] * 7177305999237621L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

