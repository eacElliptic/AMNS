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
	tmp_q[0] = ((uint64_t)op[0] * 14956118940274522428UL) + ((((uint64_t)op[1] * 10456490974239222630UL) + ((uint64_t)op[2] * 11934202902924023557UL) + ((uint64_t)op[3] * 4818345728764407180UL) + ((uint64_t)op[4] * 9502947831648082818UL) + ((uint64_t)op[5] * 14811373265844506896UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 14811373265844506896UL) + ((uint64_t)op[1] * 14956118940274522428UL) + ((((uint64_t)op[2] * 10456490974239222630UL) + ((uint64_t)op[3] * 11934202902924023557UL) + ((uint64_t)op[4] * 4818345728764407180UL) + ((uint64_t)op[5] * 9502947831648082818UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 9502947831648082818UL) + ((uint64_t)op[1] * 14811373265844506896UL) + ((uint64_t)op[2] * 14956118940274522428UL) + ((((uint64_t)op[3] * 10456490974239222630UL) + ((uint64_t)op[4] * 11934202902924023557UL) + ((uint64_t)op[5] * 4818345728764407180UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 4818345728764407180UL) + ((uint64_t)op[1] * 9502947831648082818UL) + ((uint64_t)op[2] * 14811373265844506896UL) + ((uint64_t)op[3] * 14956118940274522428UL) + ((((uint64_t)op[4] * 10456490974239222630UL) + ((uint64_t)op[5] * 11934202902924023557UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 11934202902924023557UL) + ((uint64_t)op[1] * 4818345728764407180UL) + ((uint64_t)op[2] * 9502947831648082818UL) + ((uint64_t)op[3] * 14811373265844506896UL) + ((uint64_t)op[4] * 14956118940274522428UL) + ((uint64_t)op[5] * 5524015224701435342UL);
	tmp_q[5] = ((uint64_t)op[0] * 10456490974239222630UL) + ((uint64_t)op[1] * 11934202902924023557UL) + ((uint64_t)op[2] * 4818345728764407180UL) + ((uint64_t)op[3] * 9502947831648082818UL) + ((uint64_t)op[4] * 14811373265844506896UL) + ((uint64_t)op[5] * 14956118940274522428UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1479087194L) - ((-((int128)tmp_q[1] * 1463514672L) + ((int128)tmp_q[2] * 836680196L) + ((int128)tmp_q[3] * 1068948206L) - ((int128)tmp_q[4] * 1594745809L) - ((int128)tmp_q[5] * 141502588L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 141502588L) + ((int128)tmp_q[1] * 1479087194L) - ((-((int128)tmp_q[2] * 1463514672L) + ((int128)tmp_q[3] * 836680196L) + ((int128)tmp_q[4] * 1068948206L) - ((int128)tmp_q[5] * 1594745809L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1594745809L) - ((int128)tmp_q[1] * 141502588L) + ((int128)tmp_q[2] * 1479087194L) - ((-((int128)tmp_q[3] * 1463514672L) + ((int128)tmp_q[4] * 836680196L) + ((int128)tmp_q[5] * 1068948206L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1068948206L) - ((int128)tmp_q[1] * 1594745809L) - ((int128)tmp_q[2] * 141502588L) + ((int128)tmp_q[3] * 1479087194L) - ((-((int128)tmp_q[4] * 1463514672L) + ((int128)tmp_q[5] * 836680196L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 836680196L) + ((int128)tmp_q[1] * 1068948206L) - ((int128)tmp_q[2] * 1594745809L) - ((int128)tmp_q[3] * 141502588L) + ((int128)tmp_q[4] * 1479087194L) + ((int128)tmp_q[5] * 4390544016L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1463514672L) + ((int128)tmp_q[1] * 836680196L) + ((int128)tmp_q[2] * 1068948206L) - ((int128)tmp_q[3] * 1594745809L) - ((int128)tmp_q[4] * 141502588L) + ((int128)tmp_q[5] * 1479087194L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

