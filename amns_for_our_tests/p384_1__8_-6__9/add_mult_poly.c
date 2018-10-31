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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11617722941532875511UL) + ((((uint64_t)op[1] * 16811528191888812485UL) + ((uint64_t)op[2] * 14570706581193857442UL) + ((uint64_t)op[3] * 6157234916588657861UL) + ((uint64_t)op[4] * 15400984152547628956UL) + ((uint64_t)op[5] * 18429584488250575437UL) + ((uint64_t)op[6] * 17291877074909172719UL) + ((uint64_t)op[7] * 16878977632347238449UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16878977632347238449UL) + ((uint64_t)op[1] * 11617722941532875511UL) + ((((uint64_t)op[2] * 16811528191888812485UL) + ((uint64_t)op[3] * 14570706581193857442UL) + ((uint64_t)op[4] * 6157234916588657861UL) + ((uint64_t)op[5] * 15400984152547628956UL) + ((uint64_t)op[6] * 18429584488250575437UL) + ((uint64_t)op[7] * 17291877074909172719UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 17291877074909172719UL) + ((uint64_t)op[1] * 16878977632347238449UL) + ((uint64_t)op[2] * 11617722941532875511UL) + ((((uint64_t)op[3] * 16811528191888812485UL) + ((uint64_t)op[4] * 14570706581193857442UL) + ((uint64_t)op[5] * 6157234916588657861UL) + ((uint64_t)op[6] * 15400984152547628956UL) + ((uint64_t)op[7] * 18429584488250575437UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 18429584488250575437UL) + ((uint64_t)op[1] * 17291877074909172719UL) + ((uint64_t)op[2] * 16878977632347238449UL) + ((uint64_t)op[3] * 11617722941532875511UL) + ((((uint64_t)op[4] * 16811528191888812485UL) + ((uint64_t)op[5] * 14570706581193857442UL) + ((uint64_t)op[6] * 6157234916588657861UL) + ((uint64_t)op[7] * 15400984152547628956UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 15400984152547628956UL) + ((uint64_t)op[1] * 18429584488250575437UL) + ((uint64_t)op[2] * 17291877074909172719UL) + ((uint64_t)op[3] * 16878977632347238449UL) + ((uint64_t)op[4] * 11617722941532875511UL) + ((((uint64_t)op[5] * 16811528191888812485UL) + ((uint64_t)op[6] * 14570706581193857442UL) + ((uint64_t)op[7] * 6157234916588657861UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 6157234916588657861UL) + ((uint64_t)op[1] * 15400984152547628956UL) + ((uint64_t)op[2] * 18429584488250575437UL) + ((uint64_t)op[3] * 17291877074909172719UL) + ((uint64_t)op[4] * 16878977632347238449UL) + ((uint64_t)op[5] * 11617722941532875511UL) + ((((uint64_t)op[6] * 16811528191888812485UL) + ((uint64_t)op[7] * 14570706581193857442UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 14570706581193857442UL) + ((uint64_t)op[1] * 6157234916588657861UL) + ((uint64_t)op[2] * 15400984152547628956UL) + ((uint64_t)op[3] * 18429584488250575437UL) + ((uint64_t)op[4] * 17291877074909172719UL) + ((uint64_t)op[5] * 16878977632347238449UL) + ((uint64_t)op[6] * 11617722941532875511UL) + ((uint64_t)op[7] * 9811295290924434786UL);
	tmp_q[7] = ((uint64_t)op[0] * 16811528191888812485UL) + ((uint64_t)op[1] * 14570706581193857442UL) + ((uint64_t)op[2] * 6157234916588657861UL) + ((uint64_t)op[3] * 15400984152547628956UL) + ((uint64_t)op[4] * 18429584488250575437UL) + ((uint64_t)op[5] * 17291877074909172719UL) + ((uint64_t)op[6] * 16878977632347238449UL) + ((uint64_t)op[7] * 11617722941532875511UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 146854571124185L) - ((-((int128)tmp_q[1] * 340269465044L) - ((int128)tmp_q[2] * 2137298416292L) + ((int128)tmp_q[3] * 9453567338900L) + ((int128)tmp_q[4] * 63018474410805L) - ((int128)tmp_q[5] * 120573694899098L) + ((int128)tmp_q[6] * 149925550693366L) - ((int128)tmp_q[7] * 14423918401507L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 14423918401507L) - ((int128)tmp_q[1] * 146854571124185L) - ((-((int128)tmp_q[2] * 340269465044L) - ((int128)tmp_q[3] * 2137298416292L) + ((int128)tmp_q[4] * 9453567338900L) + ((int128)tmp_q[5] * 63018474410805L) - ((int128)tmp_q[6] * 120573694899098L) + ((int128)tmp_q[7] * 149925550693366L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 149925550693366L) - ((int128)tmp_q[1] * 14423918401507L) - ((int128)tmp_q[2] * 146854571124185L) - ((-((int128)tmp_q[3] * 340269465044L) - ((int128)tmp_q[4] * 2137298416292L) + ((int128)tmp_q[5] * 9453567338900L) + ((int128)tmp_q[6] * 63018474410805L) - ((int128)tmp_q[7] * 120573694899098L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 120573694899098L) + ((int128)tmp_q[1] * 149925550693366L) - ((int128)tmp_q[2] * 14423918401507L) - ((int128)tmp_q[3] * 146854571124185L) - ((-((int128)tmp_q[4] * 340269465044L) - ((int128)tmp_q[5] * 2137298416292L) + ((int128)tmp_q[6] * 9453567338900L) + ((int128)tmp_q[7] * 63018474410805L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 63018474410805L) - ((int128)tmp_q[1] * 120573694899098L) + ((int128)tmp_q[2] * 149925550693366L) - ((int128)tmp_q[3] * 14423918401507L) - ((int128)tmp_q[4] * 146854571124185L) - ((-((int128)tmp_q[5] * 340269465044L) - ((int128)tmp_q[6] * 2137298416292L) + ((int128)tmp_q[7] * 9453567338900L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 9453567338900L) + ((int128)tmp_q[1] * 63018474410805L) - ((int128)tmp_q[2] * 120573694899098L) + ((int128)tmp_q[3] * 149925550693366L) - ((int128)tmp_q[4] * 14423918401507L) - ((int128)tmp_q[5] * 146854571124185L) - ((-((int128)tmp_q[6] * 340269465044L) - ((int128)tmp_q[7] * 2137298416292L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 2137298416292L) + ((int128)tmp_q[1] * 9453567338900L) + ((int128)tmp_q[2] * 63018474410805L) - ((int128)tmp_q[3] * 120573694899098L) + ((int128)tmp_q[4] * 149925550693366L) - ((int128)tmp_q[5] * 14423918401507L) - ((int128)tmp_q[6] * 146854571124185L) + ((int128)tmp_q[7] * 2041616790264L);
	tmp_zero[7] = -((int128)tmp_q[0] * 340269465044L) - ((int128)tmp_q[1] * 2137298416292L) + ((int128)tmp_q[2] * 9453567338900L) + ((int128)tmp_q[3] * 63018474410805L) - ((int128)tmp_q[4] * 120573694899098L) + ((int128)tmp_q[5] * 149925550693366L) - ((int128)tmp_q[6] * 14423918401507L) - ((int128)tmp_q[7] * 146854571124185L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

