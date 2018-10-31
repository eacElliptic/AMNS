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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1057039848886429626UL) + ((((uint64_t)op[1] * 10503580497126942092UL) + ((uint64_t)op[2] * 3285914013282360846UL) + ((uint64_t)op[3] * 9303322612983703062UL) + ((uint64_t)op[4] * 17505502466521394923UL) + ((uint64_t)op[5] * 3230433559227871469UL) + ((uint64_t)op[6] * 2910468138164965403UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 2910468138164965403UL) + ((uint64_t)op[1] * 1057039848886429626UL) + ((((uint64_t)op[2] * 10503580497126942092UL) + ((uint64_t)op[3] * 3285914013282360846UL) + ((uint64_t)op[4] * 9303322612983703062UL) + ((uint64_t)op[5] * 17505502466521394923UL) + ((uint64_t)op[6] * 3230433559227871469UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 3230433559227871469UL) + ((uint64_t)op[1] * 2910468138164965403UL) + ((uint64_t)op[2] * 1057039848886429626UL) + ((((uint64_t)op[3] * 10503580497126942092UL) + ((uint64_t)op[4] * 3285914013282360846UL) + ((uint64_t)op[5] * 9303322612983703062UL) + ((uint64_t)op[6] * 17505502466521394923UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 17505502466521394923UL) + ((uint64_t)op[1] * 3230433559227871469UL) + ((uint64_t)op[2] * 2910468138164965403UL) + ((uint64_t)op[3] * 1057039848886429626UL) + ((((uint64_t)op[4] * 10503580497126942092UL) + ((uint64_t)op[5] * 3285914013282360846UL) + ((uint64_t)op[6] * 9303322612983703062UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 9303322612983703062UL) + ((uint64_t)op[1] * 17505502466521394923UL) + ((uint64_t)op[2] * 3230433559227871469UL) + ((uint64_t)op[3] * 2910468138164965403UL) + ((uint64_t)op[4] * 1057039848886429626UL) + ((((uint64_t)op[5] * 10503580497126942092UL) + ((uint64_t)op[6] * 3285914013282360846UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 3285914013282360846UL) + ((uint64_t)op[1] * 9303322612983703062UL) + ((uint64_t)op[2] * 17505502466521394923UL) + ((uint64_t)op[3] * 3230433559227871469UL) + ((uint64_t)op[4] * 2910468138164965403UL) + ((uint64_t)op[5] * 1057039848886429626UL) + ((uint64_t)op[6] * 13063997417671274660UL);
	tmp_q[6] = ((uint64_t)op[0] * 10503580497126942092UL) + ((uint64_t)op[1] * 3285914013282360846UL) + ((uint64_t)op[2] * 9303322612983703062UL) + ((uint64_t)op[3] * 17505502466521394923UL) + ((uint64_t)op[4] * 3230433559227871469UL) + ((uint64_t)op[5] * 2910468138164965403UL) + ((uint64_t)op[6] * 1057039848886429626UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 23775499538L) + ((((int128)tmp_q[1] * 19615527813L) + ((int128)tmp_q[2] * 25971920435L) - ((int128)tmp_q[3] * 23772878141L) + ((int128)tmp_q[4] * 28199379802L) - ((int128)tmp_q[5] * 31551733341L) + ((int128)tmp_q[6] * 13920978737L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 13920978737L) - ((int128)tmp_q[1] * 23775499538L) + ((((int128)tmp_q[2] * 19615527813L) + ((int128)tmp_q[3] * 25971920435L) - ((int128)tmp_q[4] * 23772878141L) + ((int128)tmp_q[5] * 28199379802L) - ((int128)tmp_q[6] * 31551733341L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 31551733341L) + ((int128)tmp_q[1] * 13920978737L) - ((int128)tmp_q[2] * 23775499538L) + ((((int128)tmp_q[3] * 19615527813L) + ((int128)tmp_q[4] * 25971920435L) - ((int128)tmp_q[5] * 23772878141L) + ((int128)tmp_q[6] * 28199379802L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 28199379802L) - ((int128)tmp_q[1] * 31551733341L) + ((int128)tmp_q[2] * 13920978737L) - ((int128)tmp_q[3] * 23775499538L) + ((((int128)tmp_q[4] * 19615527813L) + ((int128)tmp_q[5] * 25971920435L) - ((int128)tmp_q[6] * 23772878141L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 23772878141L) + ((int128)tmp_q[1] * 28199379802L) - ((int128)tmp_q[2] * 31551733341L) + ((int128)tmp_q[3] * 13920978737L) - ((int128)tmp_q[4] * 23775499538L) + ((((int128)tmp_q[5] * 19615527813L) + ((int128)tmp_q[6] * 25971920435L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 25971920435L) - ((int128)tmp_q[1] * 23772878141L) + ((int128)tmp_q[2] * 28199379802L) - ((int128)tmp_q[3] * 31551733341L) + ((int128)tmp_q[4] * 13920978737L) - ((int128)tmp_q[5] * 23775499538L) + ((int128)tmp_q[6] * 58846583439L);
	tmp_zero[6] = ((int128)tmp_q[0] * 19615527813L) + ((int128)tmp_q[1] * 25971920435L) - ((int128)tmp_q[2] * 23772878141L) + ((int128)tmp_q[3] * 28199379802L) - ((int128)tmp_q[4] * 31551733341L) + ((int128)tmp_q[5] * 13920978737L) - ((int128)tmp_q[6] * 23775499538L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

