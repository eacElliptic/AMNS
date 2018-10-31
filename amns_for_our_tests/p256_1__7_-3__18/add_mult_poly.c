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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8016070588041942479UL) + ((((uint64_t)op[1] * 18157510457872760563UL) + ((uint64_t)op[2] * 17198345207836073112UL) + ((uint64_t)op[3] * 8712575983539471666UL) + ((uint64_t)op[4] * 15634368577905084684UL) + ((uint64_t)op[5] * 2808138930377981944UL) + ((uint64_t)op[6] * 16385699860357762861UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16385699860357762861UL) + ((uint64_t)op[1] * 8016070588041942479UL) + ((((uint64_t)op[2] * 18157510457872760563UL) + ((uint64_t)op[3] * 17198345207836073112UL) + ((uint64_t)op[4] * 8712575983539471666UL) + ((uint64_t)op[5] * 15634368577905084684UL) + ((uint64_t)op[6] * 2808138930377981944UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 2808138930377981944UL) + ((uint64_t)op[1] * 16385699860357762861UL) + ((uint64_t)op[2] * 8016070588041942479UL) + ((((uint64_t)op[3] * 18157510457872760563UL) + ((uint64_t)op[4] * 17198345207836073112UL) + ((uint64_t)op[5] * 8712575983539471666UL) + ((uint64_t)op[6] * 15634368577905084684UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 15634368577905084684UL) + ((uint64_t)op[1] * 2808138930377981944UL) + ((uint64_t)op[2] * 16385699860357762861UL) + ((uint64_t)op[3] * 8016070588041942479UL) + ((((uint64_t)op[4] * 18157510457872760563UL) + ((uint64_t)op[5] * 17198345207836073112UL) + ((uint64_t)op[6] * 8712575983539471666UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 8712575983539471666UL) + ((uint64_t)op[1] * 15634368577905084684UL) + ((uint64_t)op[2] * 2808138930377981944UL) + ((uint64_t)op[3] * 16385699860357762861UL) + ((uint64_t)op[4] * 8016070588041942479UL) + ((((uint64_t)op[5] * 18157510457872760563UL) + ((uint64_t)op[6] * 17198345207836073112UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 17198345207836073112UL) + ((uint64_t)op[1] * 8712575983539471666UL) + ((uint64_t)op[2] * 15634368577905084684UL) + ((uint64_t)op[3] * 2808138930377981944UL) + ((uint64_t)op[4] * 16385699860357762861UL) + ((uint64_t)op[5] * 8016070588041942479UL) + ((uint64_t)op[6] * 867700847510373159UL);
	tmp_q[6] = ((uint64_t)op[0] * 18157510457872760563UL) + ((uint64_t)op[1] * 17198345207836073112UL) + ((uint64_t)op[2] * 8712575983539471666UL) + ((uint64_t)op[3] * 15634368577905084684UL) + ((uint64_t)op[4] * 2808138930377981944UL) + ((uint64_t)op[5] * 16385699860357762861UL) + ((uint64_t)op[6] * 8016070588041942479UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 39670660835L) - ((((int128)tmp_q[1] * 26102760231L) - ((int128)tmp_q[2] * 59152951840L) - ((int128)tmp_q[3] * 15989836815L) + ((int128)tmp_q[4] * 50173650639L) + ((int128)tmp_q[5] * 8241069510L) + ((int128)tmp_q[6] * 38380259597L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 38380259597L) - ((int128)tmp_q[1] * 39670660835L) - ((((int128)tmp_q[2] * 26102760231L) - ((int128)tmp_q[3] * 59152951840L) - ((int128)tmp_q[4] * 15989836815L) + ((int128)tmp_q[5] * 50173650639L) + ((int128)tmp_q[6] * 8241069510L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 8241069510L) + ((int128)tmp_q[1] * 38380259597L) - ((int128)tmp_q[2] * 39670660835L) - ((((int128)tmp_q[3] * 26102760231L) - ((int128)tmp_q[4] * 59152951840L) - ((int128)tmp_q[5] * 15989836815L) + ((int128)tmp_q[6] * 50173650639L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 50173650639L) + ((int128)tmp_q[1] * 8241069510L) + ((int128)tmp_q[2] * 38380259597L) - ((int128)tmp_q[3] * 39670660835L) - ((((int128)tmp_q[4] * 26102760231L) - ((int128)tmp_q[5] * 59152951840L) - ((int128)tmp_q[6] * 15989836815L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 15989836815L) + ((int128)tmp_q[1] * 50173650639L) + ((int128)tmp_q[2] * 8241069510L) + ((int128)tmp_q[3] * 38380259597L) - ((int128)tmp_q[4] * 39670660835L) - ((((int128)tmp_q[5] * 26102760231L) - ((int128)tmp_q[6] * 59152951840L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 59152951840L) - ((int128)tmp_q[1] * 15989836815L) + ((int128)tmp_q[2] * 50173650639L) + ((int128)tmp_q[3] * 8241069510L) + ((int128)tmp_q[4] * 38380259597L) - ((int128)tmp_q[5] * 39670660835L) - ((int128)tmp_q[6] * 78308280693L);
	tmp_zero[6] = ((int128)tmp_q[0] * 26102760231L) - ((int128)tmp_q[1] * 59152951840L) - ((int128)tmp_q[2] * 15989836815L) + ((int128)tmp_q[3] * 50173650639L) + ((int128)tmp_q[4] * 8241069510L) + ((int128)tmp_q[5] * 38380259597L) - ((int128)tmp_q[6] * 39670660835L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

