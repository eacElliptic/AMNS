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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 207000411646587965UL) + ((((uint64_t)op[1] * 13657861021357832552UL) + ((uint64_t)op[2] * 5877309459636118592UL) + ((uint64_t)op[3] * 2555195175938225946UL) + ((uint64_t)op[4] * 3676128987081672850UL) + ((uint64_t)op[5] * 17695232017106501378UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 17695232017106501378UL) + ((uint64_t)op[1] * 207000411646587965UL) + ((((uint64_t)op[2] * 13657861021357832552UL) + ((uint64_t)op[3] * 5877309459636118592UL) + ((uint64_t)op[4] * 2555195175938225946UL) + ((uint64_t)op[5] * 3676128987081672850UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 3676128987081672850UL) + ((uint64_t)op[1] * 17695232017106501378UL) + ((uint64_t)op[2] * 207000411646587965UL) + ((((uint64_t)op[3] * 13657861021357832552UL) + ((uint64_t)op[4] * 5877309459636118592UL) + ((uint64_t)op[5] * 2555195175938225946UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 2555195175938225946UL) + ((uint64_t)op[1] * 3676128987081672850UL) + ((uint64_t)op[2] * 17695232017106501378UL) + ((uint64_t)op[3] * 207000411646587965UL) + ((((uint64_t)op[4] * 13657861021357832552UL) + ((uint64_t)op[5] * 5877309459636118592UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 5877309459636118592UL) + ((uint64_t)op[1] * 2555195175938225946UL) + ((uint64_t)op[2] * 3676128987081672850UL) + ((uint64_t)op[3] * 17695232017106501378UL) + ((uint64_t)op[4] * 207000411646587965UL) + ((uint64_t)op[5] * 4080094916654394424UL);
	tmp_q[5] = ((uint64_t)op[0] * 13657861021357832552UL) + ((uint64_t)op[1] * 5877309459636118592UL) + ((uint64_t)op[2] * 2555195175938225946UL) + ((uint64_t)op[3] * 3676128987081672850UL) + ((uint64_t)op[4] * 17695232017106501378UL) + ((uint64_t)op[5] * 207000411646587965UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 701625969511L) + ((((int128)tmp_q[1] * 5573116323728L) + ((int128)tmp_q[2] * 918091570444L) + ((int128)tmp_q[3] * 3473515248642L) - ((int128)tmp_q[4] * 896009224506L) - ((int128)tmp_q[5] * 1547222595998L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 1547222595998L) + ((int128)tmp_q[1] * 701625969511L) + ((((int128)tmp_q[2] * 5573116323728L) + ((int128)tmp_q[3] * 918091570444L) + ((int128)tmp_q[4] * 3473515248642L) - ((int128)tmp_q[5] * 896009224506L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 896009224506L) - ((int128)tmp_q[1] * 1547222595998L) + ((int128)tmp_q[2] * 701625969511L) + ((((int128)tmp_q[3] * 5573116323728L) + ((int128)tmp_q[4] * 918091570444L) + ((int128)tmp_q[5] * 3473515248642L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 3473515248642L) - ((int128)tmp_q[1] * 896009224506L) - ((int128)tmp_q[2] * 1547222595998L) + ((int128)tmp_q[3] * 701625969511L) + ((((int128)tmp_q[4] * 5573116323728L) + ((int128)tmp_q[5] * 918091570444L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 918091570444L) + ((int128)tmp_q[1] * 3473515248642L) - ((int128)tmp_q[2] * 896009224506L) - ((int128)tmp_q[3] * 1547222595998L) + ((int128)tmp_q[4] * 701625969511L) + ((int128)tmp_q[5] * 16719348971184L);
	tmp_zero[5] = ((int128)tmp_q[0] * 5573116323728L) + ((int128)tmp_q[1] * 918091570444L) + ((int128)tmp_q[2] * 3473515248642L) - ((int128)tmp_q[3] * 896009224506L) - ((int128)tmp_q[4] * 1547222595998L) + ((int128)tmp_q[5] * 701625969511L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

