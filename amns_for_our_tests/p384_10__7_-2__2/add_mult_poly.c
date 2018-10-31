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
	tmp_q[0] = ((uint64_t)op[0] * 14874681314762754793UL) + ((((uint64_t)op[1] * 18160695831890842531UL) + ((uint64_t)op[2] * 13676226128641456744UL) + ((uint64_t)op[3] * 9831571049943062898UL) + ((uint64_t)op[4] * 12747208419559577652UL) + ((uint64_t)op[5] * 4316597382310217839UL) + ((uint64_t)op[6] * 10217220473946984383UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 10217220473946984383UL) + ((uint64_t)op[1] * 14874681314762754793UL) + ((((uint64_t)op[2] * 18160695831890842531UL) + ((uint64_t)op[3] * 13676226128641456744UL) + ((uint64_t)op[4] * 9831571049943062898UL) + ((uint64_t)op[5] * 12747208419559577652UL) + ((uint64_t)op[6] * 4316597382310217839UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4316597382310217839UL) + ((uint64_t)op[1] * 10217220473946984383UL) + ((uint64_t)op[2] * 14874681314762754793UL) + ((((uint64_t)op[3] * 18160695831890842531UL) + ((uint64_t)op[4] * 13676226128641456744UL) + ((uint64_t)op[5] * 9831571049943062898UL) + ((uint64_t)op[6] * 12747208419559577652UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 12747208419559577652UL) + ((uint64_t)op[1] * 4316597382310217839UL) + ((uint64_t)op[2] * 10217220473946984383UL) + ((uint64_t)op[3] * 14874681314762754793UL) + ((((uint64_t)op[4] * 18160695831890842531UL) + ((uint64_t)op[5] * 13676226128641456744UL) + ((uint64_t)op[6] * 9831571049943062898UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 9831571049943062898UL) + ((uint64_t)op[1] * 12747208419559577652UL) + ((uint64_t)op[2] * 4316597382310217839UL) + ((uint64_t)op[3] * 10217220473946984383UL) + ((uint64_t)op[4] * 14874681314762754793UL) + ((((uint64_t)op[5] * 18160695831890842531UL) + ((uint64_t)op[6] * 13676226128641456744UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 13676226128641456744UL) + ((uint64_t)op[1] * 9831571049943062898UL) + ((uint64_t)op[2] * 12747208419559577652UL) + ((uint64_t)op[3] * 4316597382310217839UL) + ((uint64_t)op[4] * 10217220473946984383UL) + ((uint64_t)op[5] * 14874681314762754793UL) + ((uint64_t)op[6] * 572096483637418170UL);
	tmp_q[6] = ((uint64_t)op[0] * 18160695831890842531UL) + ((uint64_t)op[1] * 13676226128641456744UL) + ((uint64_t)op[2] * 9831571049943062898UL) + ((uint64_t)op[3] * 12747208419559577652UL) + ((uint64_t)op[4] * 4316597382310217839UL) + ((uint64_t)op[5] * 10217220473946984383UL) + ((uint64_t)op[6] * 14874681314762754793UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 18337291199286301L) - ((((int128)tmp_q[1] * 15420138742655248L) + ((int128)tmp_q[2] * 889438059707314L) + ((int128)tmp_q[3] * 1242529167527209L) + ((int128)tmp_q[4] * 3642825764878579L) + ((int128)tmp_q[5] * 11292250417308784L) + ((int128)tmp_q[6] * 8767258154805881L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 8767258154805881L) + ((int128)tmp_q[1] * 18337291199286301L) - ((((int128)tmp_q[2] * 15420138742655248L) + ((int128)tmp_q[3] * 889438059707314L) + ((int128)tmp_q[4] * 1242529167527209L) + ((int128)tmp_q[5] * 3642825764878579L) + ((int128)tmp_q[6] * 11292250417308784L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 11292250417308784L) + ((int128)tmp_q[1] * 8767258154805881L) + ((int128)tmp_q[2] * 18337291199286301L) - ((((int128)tmp_q[3] * 15420138742655248L) + ((int128)tmp_q[4] * 889438059707314L) + ((int128)tmp_q[5] * 1242529167527209L) + ((int128)tmp_q[6] * 3642825764878579L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 3642825764878579L) + ((int128)tmp_q[1] * 11292250417308784L) + ((int128)tmp_q[2] * 8767258154805881L) + ((int128)tmp_q[3] * 18337291199286301L) - ((((int128)tmp_q[4] * 15420138742655248L) + ((int128)tmp_q[5] * 889438059707314L) + ((int128)tmp_q[6] * 1242529167527209L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 1242529167527209L) + ((int128)tmp_q[1] * 3642825764878579L) + ((int128)tmp_q[2] * 11292250417308784L) + ((int128)tmp_q[3] * 8767258154805881L) + ((int128)tmp_q[4] * 18337291199286301L) - ((((int128)tmp_q[5] * 15420138742655248L) + ((int128)tmp_q[6] * 889438059707314L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 889438059707314L) + ((int128)tmp_q[1] * 1242529167527209L) + ((int128)tmp_q[2] * 3642825764878579L) + ((int128)tmp_q[3] * 11292250417308784L) + ((int128)tmp_q[4] * 8767258154805881L) + ((int128)tmp_q[5] * 18337291199286301L) - ((int128)tmp_q[6] * 30840277485310496L);
	tmp_zero[6] = ((int128)tmp_q[0] * 15420138742655248L) + ((int128)tmp_q[1] * 889438059707314L) + ((int128)tmp_q[2] * 1242529167527209L) + ((int128)tmp_q[3] * 3642825764878579L) + ((int128)tmp_q[4] * 11292250417308784L) + ((int128)tmp_q[5] * 8767258154805881L) + ((int128)tmp_q[6] * 18337291199286301L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

