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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14625558725039029793UL) + ((((uint64_t)op[1] * 5171730250707446614UL) + ((uint64_t)op[2] * 9568679680864233134UL) + ((uint64_t)op[3] * 9273387351002546400UL) + ((uint64_t)op[4] * 11823454437519385156UL) + ((uint64_t)op[5] * 6196589352587992976UL) + ((uint64_t)op[6] * 5585154714805464070UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 5585154714805464070UL) + ((uint64_t)op[1] * 14625558725039029793UL) + ((((uint64_t)op[2] * 5171730250707446614UL) + ((uint64_t)op[3] * 9568679680864233134UL) + ((uint64_t)op[4] * 9273387351002546400UL) + ((uint64_t)op[5] * 11823454437519385156UL) + ((uint64_t)op[6] * 6196589352587992976UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 6196589352587992976UL) + ((uint64_t)op[1] * 5585154714805464070UL) + ((uint64_t)op[2] * 14625558725039029793UL) + ((((uint64_t)op[3] * 5171730250707446614UL) + ((uint64_t)op[4] * 9568679680864233134UL) + ((uint64_t)op[5] * 9273387351002546400UL) + ((uint64_t)op[6] * 11823454437519385156UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 11823454437519385156UL) + ((uint64_t)op[1] * 6196589352587992976UL) + ((uint64_t)op[2] * 5585154714805464070UL) + ((uint64_t)op[3] * 14625558725039029793UL) + ((((uint64_t)op[4] * 5171730250707446614UL) + ((uint64_t)op[5] * 9568679680864233134UL) + ((uint64_t)op[6] * 9273387351002546400UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 9273387351002546400UL) + ((uint64_t)op[1] * 11823454437519385156UL) + ((uint64_t)op[2] * 6196589352587992976UL) + ((uint64_t)op[3] * 5585154714805464070UL) + ((uint64_t)op[4] * 14625558725039029793UL) + ((((uint64_t)op[5] * 5171730250707446614UL) + ((uint64_t)op[6] * 9568679680864233134UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 9568679680864233134UL) + ((uint64_t)op[1] * 9273387351002546400UL) + ((uint64_t)op[2] * 11823454437519385156UL) + ((uint64_t)op[3] * 6196589352587992976UL) + ((uint64_t)op[4] * 5585154714805464070UL) + ((uint64_t)op[5] * 14625558725039029793UL) + ((uint64_t)op[6] * 691376392466976934UL);
	tmp_q[6] = ((uint64_t)op[0] * 5171730250707446614UL) + ((uint64_t)op[1] * 9568679680864233134UL) + ((uint64_t)op[2] * 9273387351002546400UL) + ((uint64_t)op[3] * 11823454437519385156UL) + ((uint64_t)op[4] * 6196589352587992976UL) + ((uint64_t)op[5] * 5585154714805464070UL) + ((uint64_t)op[6] * 14625558725039029793UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 42083925377L) - ((((int128)tmp_q[1] * 42855841862L) - ((int128)tmp_q[2] * 23345497238L) + ((int128)tmp_q[3] * 33471279736L) + ((int128)tmp_q[4] * 33742320336L) + ((int128)tmp_q[5] * 11822898148L) + ((int128)tmp_q[6] * 34535492038L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 34535492038L) - ((int128)tmp_q[1] * 42083925377L) - ((((int128)tmp_q[2] * 42855841862L) - ((int128)tmp_q[3] * 23345497238L) + ((int128)tmp_q[4] * 33471279736L) + ((int128)tmp_q[5] * 33742320336L) + ((int128)tmp_q[6] * 11822898148L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 11822898148L) + ((int128)tmp_q[1] * 34535492038L) - ((int128)tmp_q[2] * 42083925377L) - ((((int128)tmp_q[3] * 42855841862L) - ((int128)tmp_q[4] * 23345497238L) + ((int128)tmp_q[5] * 33471279736L) + ((int128)tmp_q[6] * 33742320336L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 33742320336L) + ((int128)tmp_q[1] * 11822898148L) + ((int128)tmp_q[2] * 34535492038L) - ((int128)tmp_q[3] * 42083925377L) - ((((int128)tmp_q[4] * 42855841862L) - ((int128)tmp_q[5] * 23345497238L) + ((int128)tmp_q[6] * 33471279736L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 33471279736L) + ((int128)tmp_q[1] * 33742320336L) + ((int128)tmp_q[2] * 11822898148L) + ((int128)tmp_q[3] * 34535492038L) - ((int128)tmp_q[4] * 42083925377L) - ((((int128)tmp_q[5] * 42855841862L) - ((int128)tmp_q[6] * 23345497238L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 23345497238L) + ((int128)tmp_q[1] * 33471279736L) + ((int128)tmp_q[2] * 33742320336L) + ((int128)tmp_q[3] * 11822898148L) + ((int128)tmp_q[4] * 34535492038L) - ((int128)tmp_q[5] * 42083925377L) - ((int128)tmp_q[6] * 299990893034L);
	tmp_zero[6] = ((int128)tmp_q[0] * 42855841862L) - ((int128)tmp_q[1] * 23345497238L) + ((int128)tmp_q[2] * 33471279736L) + ((int128)tmp_q[3] * 33742320336L) + ((int128)tmp_q[4] * 11822898148L) + ((int128)tmp_q[5] * 34535492038L) - ((int128)tmp_q[6] * 42083925377L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

