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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12961445870018685433UL) + ((((uint64_t)op[1] * 4364958576476909966UL) + ((uint64_t)op[2] * 1095846286521641593UL) + ((uint64_t)op[3] * 3100851547238619757UL) + ((uint64_t)op[4] * 7783989413292010539UL) + ((uint64_t)op[5] * 17755879977489701250UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 17755879977489701250UL) + ((uint64_t)op[1] * 12961445870018685433UL) + ((((uint64_t)op[2] * 4364958576476909966UL) + ((uint64_t)op[3] * 1095846286521641593UL) + ((uint64_t)op[4] * 3100851547238619757UL) + ((uint64_t)op[5] * 7783989413292010539UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 7783989413292010539UL) + ((uint64_t)op[1] * 17755879977489701250UL) + ((uint64_t)op[2] * 12961445870018685433UL) + ((((uint64_t)op[3] * 4364958576476909966UL) + ((uint64_t)op[4] * 1095846286521641593UL) + ((uint64_t)op[5] * 3100851547238619757UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 3100851547238619757UL) + ((uint64_t)op[1] * 7783989413292010539UL) + ((uint64_t)op[2] * 17755879977489701250UL) + ((uint64_t)op[3] * 12961445870018685433UL) + ((((uint64_t)op[4] * 4364958576476909966UL) + ((uint64_t)op[5] * 1095846286521641593UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 1095846286521641593UL) + ((uint64_t)op[1] * 3100851547238619757UL) + ((uint64_t)op[2] * 7783989413292010539UL) + ((uint64_t)op[3] * 17755879977489701250UL) + ((uint64_t)op[4] * 12961445870018685433UL) + ((uint64_t)op[5] * 16472924538105728112UL);
	tmp_q[5] = ((uint64_t)op[0] * 4364958576476909966UL) + ((uint64_t)op[1] * 1095846286521641593UL) + ((uint64_t)op[2] * 3100851547238619757UL) + ((uint64_t)op[3] * 7783989413292010539UL) + ((uint64_t)op[4] * 17755879977489701250UL) + ((uint64_t)op[5] * 12961445870018685433UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 14825037849961L) + ((((int128)tmp_q[1] * 26692735063966L) - ((int128)tmp_q[2] * 37557165891296L) - ((int128)tmp_q[3] * 19249211813823L) + ((int128)tmp_q[4] * 48428658645575L) + ((int128)tmp_q[5] * 2951422140922L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 2951422140922L) - ((int128)tmp_q[1] * 14825037849961L) + ((((int128)tmp_q[2] * 26692735063966L) - ((int128)tmp_q[3] * 37557165891296L) - ((int128)tmp_q[4] * 19249211813823L) + ((int128)tmp_q[5] * 48428658645575L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 48428658645575L) + ((int128)tmp_q[1] * 2951422140922L) - ((int128)tmp_q[2] * 14825037849961L) + ((((int128)tmp_q[3] * 26692735063966L) - ((int128)tmp_q[4] * 37557165891296L) - ((int128)tmp_q[5] * 19249211813823L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 19249211813823L) + ((int128)tmp_q[1] * 48428658645575L) + ((int128)tmp_q[2] * 2951422140922L) - ((int128)tmp_q[3] * 14825037849961L) + ((((int128)tmp_q[4] * 26692735063966L) - ((int128)tmp_q[5] * 37557165891296L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 37557165891296L) - ((int128)tmp_q[1] * 19249211813823L) + ((int128)tmp_q[2] * 48428658645575L) + ((int128)tmp_q[3] * 2951422140922L) - ((int128)tmp_q[4] * 14825037849961L) + ((int128)tmp_q[5] * 213541880511728L);
	tmp_zero[5] = ((int128)tmp_q[0] * 26692735063966L) - ((int128)tmp_q[1] * 37557165891296L) - ((int128)tmp_q[2] * 19249211813823L) + ((int128)tmp_q[3] * 48428658645575L) + ((int128)tmp_q[4] * 2951422140922L) - ((int128)tmp_q[5] * 14825037849961L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

