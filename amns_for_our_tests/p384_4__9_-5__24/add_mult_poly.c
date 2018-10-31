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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18172535352960714616UL) + ((((uint64_t)op[1] * 15463150451162056597UL) + ((uint64_t)op[2] * 6421925713528488922UL) + ((uint64_t)op[3] * 9133045761256993702UL) + ((uint64_t)op[4] * 407508576896164176UL) + ((uint64_t)op[5] * 10986848096820725223UL) + ((uint64_t)op[6] * 16710207295713611012UL) + ((uint64_t)op[7] * 1304644891464255643UL) + ((uint64_t)op[8] * 420076336158509092UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 420076336158509092UL) + ((uint64_t)op[1] * 18172535352960714616UL) + ((((uint64_t)op[2] * 15463150451162056597UL) + ((uint64_t)op[3] * 6421925713528488922UL) + ((uint64_t)op[4] * 9133045761256993702UL) + ((uint64_t)op[5] * 407508576896164176UL) + ((uint64_t)op[6] * 10986848096820725223UL) + ((uint64_t)op[7] * 16710207295713611012UL) + ((uint64_t)op[8] * 1304644891464255643UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1304644891464255643UL) + ((uint64_t)op[1] * 420076336158509092UL) + ((uint64_t)op[2] * 18172535352960714616UL) + ((((uint64_t)op[3] * 15463150451162056597UL) + ((uint64_t)op[4] * 6421925713528488922UL) + ((uint64_t)op[5] * 9133045761256993702UL) + ((uint64_t)op[6] * 407508576896164176UL) + ((uint64_t)op[7] * 10986848096820725223UL) + ((uint64_t)op[8] * 16710207295713611012UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16710207295713611012UL) + ((uint64_t)op[1] * 1304644891464255643UL) + ((uint64_t)op[2] * 420076336158509092UL) + ((uint64_t)op[3] * 18172535352960714616UL) + ((((uint64_t)op[4] * 15463150451162056597UL) + ((uint64_t)op[5] * 6421925713528488922UL) + ((uint64_t)op[6] * 9133045761256993702UL) + ((uint64_t)op[7] * 407508576896164176UL) + ((uint64_t)op[8] * 10986848096820725223UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 10986848096820725223UL) + ((uint64_t)op[1] * 16710207295713611012UL) + ((uint64_t)op[2] * 1304644891464255643UL) + ((uint64_t)op[3] * 420076336158509092UL) + ((uint64_t)op[4] * 18172535352960714616UL) + ((((uint64_t)op[5] * 15463150451162056597UL) + ((uint64_t)op[6] * 6421925713528488922UL) + ((uint64_t)op[7] * 9133045761256993702UL) + ((uint64_t)op[8] * 407508576896164176UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 407508576896164176UL) + ((uint64_t)op[1] * 10986848096820725223UL) + ((uint64_t)op[2] * 16710207295713611012UL) + ((uint64_t)op[3] * 1304644891464255643UL) + ((uint64_t)op[4] * 420076336158509092UL) + ((uint64_t)op[5] * 18172535352960714616UL) + ((((uint64_t)op[6] * 15463150451162056597UL) + ((uint64_t)op[7] * 6421925713528488922UL) + ((uint64_t)op[8] * 9133045761256993702UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 9133045761256993702UL) + ((uint64_t)op[1] * 407508576896164176UL) + ((uint64_t)op[2] * 10986848096820725223UL) + ((uint64_t)op[3] * 16710207295713611012UL) + ((uint64_t)op[4] * 1304644891464255643UL) + ((uint64_t)op[5] * 420076336158509092UL) + ((uint64_t)op[6] * 18172535352960714616UL) + ((((uint64_t)op[7] * 15463150451162056597UL) + ((uint64_t)op[8] * 6421925713528488922UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 6421925713528488922UL) + ((uint64_t)op[1] * 9133045761256993702UL) + ((uint64_t)op[2] * 407508576896164176UL) + ((uint64_t)op[3] * 10986848096820725223UL) + ((uint64_t)op[4] * 16710207295713611012UL) + ((uint64_t)op[5] * 1304644891464255643UL) + ((uint64_t)op[6] * 420076336158509092UL) + ((uint64_t)op[7] * 18172535352960714616UL) + ((uint64_t)op[8] * 14917968112737475095UL);
	tmp_q[8] = ((uint64_t)op[0] * 15463150451162056597UL) + ((uint64_t)op[1] * 6421925713528488922UL) + ((uint64_t)op[2] * 9133045761256993702UL) + ((uint64_t)op[3] * 407508576896164176UL) + ((uint64_t)op[4] * 10986848096820725223UL) + ((uint64_t)op[5] * 16710207295713611012UL) + ((uint64_t)op[6] * 1304644891464255643UL) + ((uint64_t)op[7] * 420076336158509092UL) + ((uint64_t)op[8] * 18172535352960714616UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3088720013608L) - ((-((int128)tmp_q[1] * 2779920500689L) - ((int128)tmp_q[2] * 3055878985581L) - ((int128)tmp_q[3] * 1133422641739L) + ((int128)tmp_q[4] * 2707348594728L) - ((int128)tmp_q[5] * 521972795257L) - ((int128)tmp_q[6] * 874099304139L) - ((int128)tmp_q[7] * 4413304507920L) - ((int128)tmp_q[8] * 491989106626L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 491989106626L) + ((int128)tmp_q[1] * 3088720013608L) - ((-((int128)tmp_q[2] * 2779920500689L) - ((int128)tmp_q[3] * 3055878985581L) - ((int128)tmp_q[4] * 1133422641739L) + ((int128)tmp_q[5] * 2707348594728L) - ((int128)tmp_q[6] * 521972795257L) - ((int128)tmp_q[7] * 874099304139L) - ((int128)tmp_q[8] * 4413304507920L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 4413304507920L) - ((int128)tmp_q[1] * 491989106626L) + ((int128)tmp_q[2] * 3088720013608L) - ((-((int128)tmp_q[3] * 2779920500689L) - ((int128)tmp_q[4] * 3055878985581L) - ((int128)tmp_q[5] * 1133422641739L) + ((int128)tmp_q[6] * 2707348594728L) - ((int128)tmp_q[7] * 521972795257L) - ((int128)tmp_q[8] * 874099304139L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 874099304139L) - ((int128)tmp_q[1] * 4413304507920L) - ((int128)tmp_q[2] * 491989106626L) + ((int128)tmp_q[3] * 3088720013608L) - ((-((int128)tmp_q[4] * 2779920500689L) - ((int128)tmp_q[5] * 3055878985581L) - ((int128)tmp_q[6] * 1133422641739L) + ((int128)tmp_q[7] * 2707348594728L) - ((int128)tmp_q[8] * 521972795257L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 521972795257L) - ((int128)tmp_q[1] * 874099304139L) - ((int128)tmp_q[2] * 4413304507920L) - ((int128)tmp_q[3] * 491989106626L) + ((int128)tmp_q[4] * 3088720013608L) - ((-((int128)tmp_q[5] * 2779920500689L) - ((int128)tmp_q[6] * 3055878985581L) - ((int128)tmp_q[7] * 1133422641739L) + ((int128)tmp_q[8] * 2707348594728L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 2707348594728L) - ((int128)tmp_q[1] * 521972795257L) - ((int128)tmp_q[2] * 874099304139L) - ((int128)tmp_q[3] * 4413304507920L) - ((int128)tmp_q[4] * 491989106626L) + ((int128)tmp_q[5] * 3088720013608L) - ((-((int128)tmp_q[6] * 2779920500689L) - ((int128)tmp_q[7] * 3055878985581L) - ((int128)tmp_q[8] * 1133422641739L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 1133422641739L) + ((int128)tmp_q[1] * 2707348594728L) - ((int128)tmp_q[2] * 521972795257L) - ((int128)tmp_q[3] * 874099304139L) - ((int128)tmp_q[4] * 4413304507920L) - ((int128)tmp_q[5] * 491989106626L) + ((int128)tmp_q[6] * 3088720013608L) - ((-((int128)tmp_q[7] * 2779920500689L) - ((int128)tmp_q[8] * 3055878985581L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 3055878985581L) - ((int128)tmp_q[1] * 1133422641739L) + ((int128)tmp_q[2] * 2707348594728L) - ((int128)tmp_q[3] * 521972795257L) - ((int128)tmp_q[4] * 874099304139L) - ((int128)tmp_q[5] * 4413304507920L) - ((int128)tmp_q[6] * 491989106626L) + ((int128)tmp_q[7] * 3088720013608L) + ((int128)tmp_q[8] * 13899602503445L);
	tmp_zero[8] = -((int128)tmp_q[0] * 2779920500689L) - ((int128)tmp_q[1] * 3055878985581L) - ((int128)tmp_q[2] * 1133422641739L) + ((int128)tmp_q[3] * 2707348594728L) - ((int128)tmp_q[4] * 521972795257L) - ((int128)tmp_q[5] * 874099304139L) - ((int128)tmp_q[6] * 4413304507920L) - ((int128)tmp_q[7] * 491989106626L) + ((int128)tmp_q[8] * 3088720013608L);

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
}

