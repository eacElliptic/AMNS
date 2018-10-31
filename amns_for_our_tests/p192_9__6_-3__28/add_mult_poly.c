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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3273232495486360784UL) + ((((uint64_t)op[1] * 2452362575839749371UL) + ((uint64_t)op[2] * 10358213276685351769UL) + ((uint64_t)op[3] * 7804288293212756619UL) + ((uint64_t)op[4] * 13361991645593519353UL) + ((uint64_t)op[5] * 10308219354373867131UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 10308219354373867131UL) + ((uint64_t)op[1] * 3273232495486360784UL) + ((((uint64_t)op[2] * 2452362575839749371UL) + ((uint64_t)op[3] * 10358213276685351769UL) + ((uint64_t)op[4] * 7804288293212756619UL) + ((uint64_t)op[5] * 13361991645593519353UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 13361991645593519353UL) + ((uint64_t)op[1] * 10308219354373867131UL) + ((uint64_t)op[2] * 3273232495486360784UL) + ((((uint64_t)op[3] * 2452362575839749371UL) + ((uint64_t)op[4] * 10358213276685351769UL) + ((uint64_t)op[5] * 7804288293212756619UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7804288293212756619UL) + ((uint64_t)op[1] * 13361991645593519353UL) + ((uint64_t)op[2] * 10308219354373867131UL) + ((uint64_t)op[3] * 3273232495486360784UL) + ((((uint64_t)op[4] * 2452362575839749371UL) + ((uint64_t)op[5] * 10358213276685351769UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 10358213276685351769UL) + ((uint64_t)op[1] * 7804288293212756619UL) + ((uint64_t)op[2] * 13361991645593519353UL) + ((uint64_t)op[3] * 10308219354373867131UL) + ((uint64_t)op[4] * 3273232495486360784UL) + ((uint64_t)op[5] * 11089656346190303503UL);
	tmp_q[5] = ((uint64_t)op[0] * 2452362575839749371UL) + ((uint64_t)op[1] * 10358213276685351769UL) + ((uint64_t)op[2] * 7804288293212756619UL) + ((uint64_t)op[3] * 13361991645593519353UL) + ((uint64_t)op[4] * 10308219354373867131UL) + ((uint64_t)op[5] * 3273232495486360784UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 525594884L) - ((((int128)tmp_q[1] * 915541349L) - ((int128)tmp_q[2] * 1767659309L) + ((int128)tmp_q[3] * 167651901L) + ((int128)tmp_q[4] * 226995147L) + ((int128)tmp_q[5] * 2172749237L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 2172749237L) - ((int128)tmp_q[1] * 525594884L) - ((((int128)tmp_q[2] * 915541349L) - ((int128)tmp_q[3] * 1767659309L) + ((int128)tmp_q[4] * 167651901L) + ((int128)tmp_q[5] * 226995147L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 226995147L) + ((int128)tmp_q[1] * 2172749237L) - ((int128)tmp_q[2] * 525594884L) - ((((int128)tmp_q[3] * 915541349L) - ((int128)tmp_q[4] * 1767659309L) + ((int128)tmp_q[5] * 167651901L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 167651901L) + ((int128)tmp_q[1] * 226995147L) + ((int128)tmp_q[2] * 2172749237L) - ((int128)tmp_q[3] * 525594884L) - ((((int128)tmp_q[4] * 915541349L) - ((int128)tmp_q[5] * 1767659309L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 1767659309L) + ((int128)tmp_q[1] * 167651901L) + ((int128)tmp_q[2] * 226995147L) + ((int128)tmp_q[3] * 2172749237L) - ((int128)tmp_q[4] * 525594884L) - ((int128)tmp_q[5] * 2746624047L);
	tmp_zero[5] = ((int128)tmp_q[0] * 915541349L) - ((int128)tmp_q[1] * 1767659309L) + ((int128)tmp_q[2] * 167651901L) + ((int128)tmp_q[3] * 226995147L) + ((int128)tmp_q[4] * 2172749237L) - ((int128)tmp_q[5] * 525594884L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

