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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11879720740706589755UL) + ((((uint64_t)op[1] * 14538597360121901225UL) + ((uint64_t)op[2] * 6819717601932065354UL) + ((uint64_t)op[3] * 18247351760293261487UL) + ((uint64_t)op[4] * 7333887826004046506UL) + ((uint64_t)op[5] * 18183785257552255362UL) + ((uint64_t)op[6] * 14236102892384403928UL) + ((uint64_t)op[7] * 7855547280233996288UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 7855547280233996288UL) + ((uint64_t)op[1] * 11879720740706589755UL) + ((((uint64_t)op[2] * 14538597360121901225UL) + ((uint64_t)op[3] * 6819717601932065354UL) + ((uint64_t)op[4] * 18247351760293261487UL) + ((uint64_t)op[5] * 7333887826004046506UL) + ((uint64_t)op[6] * 18183785257552255362UL) + ((uint64_t)op[7] * 14236102892384403928UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 14236102892384403928UL) + ((uint64_t)op[1] * 7855547280233996288UL) + ((uint64_t)op[2] * 11879720740706589755UL) + ((((uint64_t)op[3] * 14538597360121901225UL) + ((uint64_t)op[4] * 6819717601932065354UL) + ((uint64_t)op[5] * 18247351760293261487UL) + ((uint64_t)op[6] * 7333887826004046506UL) + ((uint64_t)op[7] * 18183785257552255362UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 18183785257552255362UL) + ((uint64_t)op[1] * 14236102892384403928UL) + ((uint64_t)op[2] * 7855547280233996288UL) + ((uint64_t)op[3] * 11879720740706589755UL) + ((((uint64_t)op[4] * 14538597360121901225UL) + ((uint64_t)op[5] * 6819717601932065354UL) + ((uint64_t)op[6] * 18247351760293261487UL) + ((uint64_t)op[7] * 7333887826004046506UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 7333887826004046506UL) + ((uint64_t)op[1] * 18183785257552255362UL) + ((uint64_t)op[2] * 14236102892384403928UL) + ((uint64_t)op[3] * 7855547280233996288UL) + ((uint64_t)op[4] * 11879720740706589755UL) + ((((uint64_t)op[5] * 14538597360121901225UL) + ((uint64_t)op[6] * 6819717601932065354UL) + ((uint64_t)op[7] * 18247351760293261487UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 18247351760293261487UL) + ((uint64_t)op[1] * 7333887826004046506UL) + ((uint64_t)op[2] * 18183785257552255362UL) + ((uint64_t)op[3] * 14236102892384403928UL) + ((uint64_t)op[4] * 7855547280233996288UL) + ((uint64_t)op[5] * 11879720740706589755UL) + ((((uint64_t)op[6] * 14538597360121901225UL) + ((uint64_t)op[7] * 6819717601932065354UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 6819717601932065354UL) + ((uint64_t)op[1] * 18247351760293261487UL) + ((uint64_t)op[2] * 7333887826004046506UL) + ((uint64_t)op[3] * 18183785257552255362UL) + ((uint64_t)op[4] * 14236102892384403928UL) + ((uint64_t)op[5] * 7855547280233996288UL) + ((uint64_t)op[6] * 11879720740706589755UL) + ((uint64_t)op[7] * 10630450646534250834UL);
	tmp_q[7] = ((uint64_t)op[0] * 14538597360121901225UL) + ((uint64_t)op[1] * 6819717601932065354UL) + ((uint64_t)op[2] * 18247351760293261487UL) + ((uint64_t)op[3] * 7333887826004046506UL) + ((uint64_t)op[4] * 18183785257552255362UL) + ((uint64_t)op[5] * 14236102892384403928UL) + ((uint64_t)op[6] * 7855547280233996288UL) + ((uint64_t)op[7] * 11879720740706589755UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3984467733557L) + ((-((int128)tmp_q[1] * 41485375477377L) + ((int128)tmp_q[2] * 51673748971212L) - ((int128)tmp_q[3] * 79721514460289L) + ((int128)tmp_q[4] * 90125397350338L) + ((int128)tmp_q[5] * 57986452960022L) - ((int128)tmp_q[6] * 82147710974414L) - ((int128)tmp_q[7] * 125998554405932L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 125998554405932L) + ((int128)tmp_q[1] * 3984467733557L) + ((-((int128)tmp_q[2] * 41485375477377L) + ((int128)tmp_q[3] * 51673748971212L) - ((int128)tmp_q[4] * 79721514460289L) + ((int128)tmp_q[5] * 90125397350338L) + ((int128)tmp_q[6] * 57986452960022L) - ((int128)tmp_q[7] * 82147710974414L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 82147710974414L) - ((int128)tmp_q[1] * 125998554405932L) + ((int128)tmp_q[2] * 3984467733557L) + ((-((int128)tmp_q[3] * 41485375477377L) + ((int128)tmp_q[4] * 51673748971212L) - ((int128)tmp_q[5] * 79721514460289L) + ((int128)tmp_q[6] * 90125397350338L) + ((int128)tmp_q[7] * 57986452960022L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 57986452960022L) - ((int128)tmp_q[1] * 82147710974414L) - ((int128)tmp_q[2] * 125998554405932L) + ((int128)tmp_q[3] * 3984467733557L) + ((-((int128)tmp_q[4] * 41485375477377L) + ((int128)tmp_q[5] * 51673748971212L) - ((int128)tmp_q[6] * 79721514460289L) + ((int128)tmp_q[7] * 90125397350338L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 90125397350338L) + ((int128)tmp_q[1] * 57986452960022L) - ((int128)tmp_q[2] * 82147710974414L) - ((int128)tmp_q[3] * 125998554405932L) + ((int128)tmp_q[4] * 3984467733557L) + ((-((int128)tmp_q[5] * 41485375477377L) + ((int128)tmp_q[6] * 51673748971212L) - ((int128)tmp_q[7] * 79721514460289L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 79721514460289L) + ((int128)tmp_q[1] * 90125397350338L) + ((int128)tmp_q[2] * 57986452960022L) - ((int128)tmp_q[3] * 82147710974414L) - ((int128)tmp_q[4] * 125998554405932L) + ((int128)tmp_q[5] * 3984467733557L) + ((-((int128)tmp_q[6] * 41485375477377L) + ((int128)tmp_q[7] * 51673748971212L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 51673748971212L) - ((int128)tmp_q[1] * 79721514460289L) + ((int128)tmp_q[2] * 90125397350338L) + ((int128)tmp_q[3] * 57986452960022L) - ((int128)tmp_q[4] * 82147710974414L) - ((int128)tmp_q[5] * 125998554405932L) + ((int128)tmp_q[6] * 3984467733557L) - ((int128)tmp_q[7] * 82970750954754L);
	tmp_zero[7] = -((int128)tmp_q[0] * 41485375477377L) + ((int128)tmp_q[1] * 51673748971212L) - ((int128)tmp_q[2] * 79721514460289L) + ((int128)tmp_q[3] * 90125397350338L) + ((int128)tmp_q[4] * 57986452960022L) - ((int128)tmp_q[5] * 82147710974414L) - ((int128)tmp_q[6] * 125998554405932L) + ((int128)tmp_q[7] * 3984467733557L);

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

