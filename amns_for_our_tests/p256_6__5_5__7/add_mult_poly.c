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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3674476817345213575UL) + ((((uint64_t)op[1] * 5607106479988302703UL) + ((uint64_t)op[2] * 11926026464258504054UL) + ((uint64_t)op[3] * 14177522570414834979UL) + ((uint64_t)op[4] * 5744739837092369156UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 5744739837092369156UL) + ((uint64_t)op[1] * 3674476817345213575UL) + ((((uint64_t)op[2] * 5607106479988302703UL) + ((uint64_t)op[3] * 11926026464258504054UL) + ((uint64_t)op[4] * 14177522570414834979UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 14177522570414834979UL) + ((uint64_t)op[1] * 5744739837092369156UL) + ((uint64_t)op[2] * 3674476817345213575UL) + ((((uint64_t)op[3] * 5607106479988302703UL) + ((uint64_t)op[4] * 11926026464258504054UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 11926026464258504054UL) + ((uint64_t)op[1] * 14177522570414834979UL) + ((uint64_t)op[2] * 5744739837092369156UL) + ((uint64_t)op[3] * 3674476817345213575UL) + ((uint64_t)op[4] * 9588788326231961899UL);
	tmp_q[4] = ((uint64_t)op[0] * 5607106479988302703UL) + ((uint64_t)op[1] * 11926026464258504054UL) + ((uint64_t)op[2] * 14177522570414834979UL) + ((uint64_t)op[3] * 5744739837092369156UL) + ((uint64_t)op[4] * 3674476817345213575UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 983500133380506L) + ((-((int128)tmp_q[1] * 649925106678301L) + ((int128)tmp_q[2] * 335422847989793L) - ((int128)tmp_q[3] * 2189289658398017L) + ((int128)tmp_q[4] * 119657016898600L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 119657016898600L) - ((int128)tmp_q[1] * 983500133380506L) + ((-((int128)tmp_q[2] * 649925106678301L) + ((int128)tmp_q[3] * 335422847989793L) - ((int128)tmp_q[4] * 2189289658398017L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 2189289658398017L) + ((int128)tmp_q[1] * 119657016898600L) - ((int128)tmp_q[2] * 983500133380506L) + ((-((int128)tmp_q[3] * 649925106678301L) + ((int128)tmp_q[4] * 335422847989793L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 335422847989793L) - ((int128)tmp_q[1] * 2189289658398017L) + ((int128)tmp_q[2] * 119657016898600L) - ((int128)tmp_q[3] * 983500133380506L) - ((int128)tmp_q[4] * 3249625533391505L);
	tmp_zero[4] = -((int128)tmp_q[0] * 649925106678301L) + ((int128)tmp_q[1] * 335422847989793L) - ((int128)tmp_q[2] * 2189289658398017L) + ((int128)tmp_q[3] * 119657016898600L) - ((int128)tmp_q[4] * 983500133380506L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

