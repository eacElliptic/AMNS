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
	tmp_q[0] = ((uint64_t)op[0] * 16844708163453615069UL) + ((((uint64_t)op[1] * 6893457011716114035UL) + ((uint64_t)op[2] * 5340817801823996773UL) + ((uint64_t)op[3] * 14952745971586699518UL) + ((uint64_t)op[4] * 8580709140303993979UL) + ((uint64_t)op[5] * 16714637191724309594UL) + ((uint64_t)op[6] * 9492214483232924109UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 9492214483232924109UL) + ((uint64_t)op[1] * 16844708163453615069UL) + ((((uint64_t)op[2] * 6893457011716114035UL) + ((uint64_t)op[3] * 5340817801823996773UL) + ((uint64_t)op[4] * 14952745971586699518UL) + ((uint64_t)op[5] * 8580709140303993979UL) + ((uint64_t)op[6] * 16714637191724309594UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 16714637191724309594UL) + ((uint64_t)op[1] * 9492214483232924109UL) + ((uint64_t)op[2] * 16844708163453615069UL) + ((((uint64_t)op[3] * 6893457011716114035UL) + ((uint64_t)op[4] * 5340817801823996773UL) + ((uint64_t)op[5] * 14952745971586699518UL) + ((uint64_t)op[6] * 8580709140303993979UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 8580709140303993979UL) + ((uint64_t)op[1] * 16714637191724309594UL) + ((uint64_t)op[2] * 9492214483232924109UL) + ((uint64_t)op[3] * 16844708163453615069UL) + ((((uint64_t)op[4] * 6893457011716114035UL) + ((uint64_t)op[5] * 5340817801823996773UL) + ((uint64_t)op[6] * 14952745971586699518UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 14952745971586699518UL) + ((uint64_t)op[1] * 8580709140303993979UL) + ((uint64_t)op[2] * 16714637191724309594UL) + ((uint64_t)op[3] * 9492214483232924109UL) + ((uint64_t)op[4] * 16844708163453615069UL) + ((((uint64_t)op[5] * 6893457011716114035UL) + ((uint64_t)op[6] * 5340817801823996773UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 5340817801823996773UL) + ((uint64_t)op[1] * 14952745971586699518UL) + ((uint64_t)op[2] * 8580709140303993979UL) + ((uint64_t)op[3] * 16714637191724309594UL) + ((uint64_t)op[4] * 9492214483232924109UL) + ((uint64_t)op[5] * 16844708163453615069UL) + ((uint64_t)op[6] * 2233626961438790489UL);
	tmp_q[6] = ((uint64_t)op[0] * 6893457011716114035UL) + ((uint64_t)op[1] * 5340817801823996773UL) + ((uint64_t)op[2] * 14952745971586699518UL) + ((uint64_t)op[3] * 8580709140303993979UL) + ((uint64_t)op[4] * 16714637191724309594UL) + ((uint64_t)op[5] * 9492214483232924109UL) + ((uint64_t)op[6] * 16844708163453615069UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 48502609915L) + ((-((int128)tmp_q[1] * 23190727080L) - ((int128)tmp_q[2] * 8291780850L) - ((int128)tmp_q[3] * 11986555999L) - ((int128)tmp_q[4] * 73879487938L) - ((int128)tmp_q[5] * 26245521622L) - ((int128)tmp_q[6] * 36738071353L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 36738071353L) - ((int128)tmp_q[1] * 48502609915L) + ((-((int128)tmp_q[2] * 23190727080L) - ((int128)tmp_q[3] * 8291780850L) - ((int128)tmp_q[4] * 11986555999L) - ((int128)tmp_q[5] * 73879487938L) - ((int128)tmp_q[6] * 26245521622L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 26245521622L) - ((int128)tmp_q[1] * 36738071353L) - ((int128)tmp_q[2] * 48502609915L) + ((-((int128)tmp_q[3] * 23190727080L) - ((int128)tmp_q[4] * 8291780850L) - ((int128)tmp_q[5] * 11986555999L) - ((int128)tmp_q[6] * 73879487938L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 73879487938L) - ((int128)tmp_q[1] * 26245521622L) - ((int128)tmp_q[2] * 36738071353L) - ((int128)tmp_q[3] * 48502609915L) + ((-((int128)tmp_q[4] * 23190727080L) - ((int128)tmp_q[5] * 8291780850L) - ((int128)tmp_q[6] * 11986555999L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 11986555999L) - ((int128)tmp_q[1] * 73879487938L) - ((int128)tmp_q[2] * 26245521622L) - ((int128)tmp_q[3] * 36738071353L) - ((int128)tmp_q[4] * 48502609915L) + ((-((int128)tmp_q[5] * 23190727080L) - ((int128)tmp_q[6] * 8291780850L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 8291780850L) - ((int128)tmp_q[1] * 11986555999L) - ((int128)tmp_q[2] * 73879487938L) - ((int128)tmp_q[3] * 26245521622L) - ((int128)tmp_q[4] * 36738071353L) - ((int128)tmp_q[5] * 48502609915L) - ((int128)tmp_q[6] * 69572181240L);
	tmp_zero[6] = -((int128)tmp_q[0] * 23190727080L) - ((int128)tmp_q[1] * 8291780850L) - ((int128)tmp_q[2] * 11986555999L) - ((int128)tmp_q[3] * 73879487938L) - ((int128)tmp_q[4] * 26245521622L) - ((int128)tmp_q[5] * 36738071353L) - ((int128)tmp_q[6] * 48502609915L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

