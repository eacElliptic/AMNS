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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5245257215339080807UL) + ((((uint64_t)op[1] * 13461415224211826485UL) + ((uint64_t)op[2] * 10468906763370886009UL) + ((uint64_t)op[3] * 16866287556698432170UL) + ((uint64_t)op[4] * 13743670162701355294UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 13743670162701355294UL) + ((uint64_t)op[1] * 5245257215339080807UL) + ((((uint64_t)op[2] * 13461415224211826485UL) + ((uint64_t)op[3] * 10468906763370886009UL) + ((uint64_t)op[4] * 16866287556698432170UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 16866287556698432170UL) + ((uint64_t)op[1] * 13743670162701355294UL) + ((uint64_t)op[2] * 5245257215339080807UL) + ((((uint64_t)op[3] * 13461415224211826485UL) + ((uint64_t)op[4] * 10468906763370886009UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 10468906763370886009UL) + ((uint64_t)op[1] * 16866287556698432170UL) + ((uint64_t)op[2] * 13743670162701355294UL) + ((uint64_t)op[3] * 5245257215339080807UL) + ((uint64_t)op[4] * 16952172749428202708UL);
	tmp_q[4] = ((uint64_t)op[0] * 13461415224211826485UL) + ((uint64_t)op[1] * 10468906763370886009UL) + ((uint64_t)op[2] * 16866287556698432170UL) + ((uint64_t)op[3] * 13743670162701355294UL) + ((uint64_t)op[4] * 5245257215339080807UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 8435347520905L) + ((-((int128)tmp_q[1] * 4969932413815L) - ((int128)tmp_q[2] * 12625207432683L) - ((int128)tmp_q[3] * 6540521952722L) - ((int128)tmp_q[4] * 792290124782L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 792290124782L) + ((int128)tmp_q[1] * 8435347520905L) + ((-((int128)tmp_q[2] * 4969932413815L) - ((int128)tmp_q[3] * 12625207432683L) - ((int128)tmp_q[4] * 6540521952722L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 6540521952722L) - ((int128)tmp_q[1] * 792290124782L) + ((int128)tmp_q[2] * 8435347520905L) + ((-((int128)tmp_q[3] * 4969932413815L) - ((int128)tmp_q[4] * 12625207432683L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 12625207432683L) - ((int128)tmp_q[1] * 6540521952722L) - ((int128)tmp_q[2] * 792290124782L) + ((int128)tmp_q[3] * 8435347520905L) - ((int128)tmp_q[4] * 19879729655260L);
	tmp_zero[4] = -((int128)tmp_q[0] * 4969932413815L) - ((int128)tmp_q[1] * 12625207432683L) - ((int128)tmp_q[2] * 6540521952722L) - ((int128)tmp_q[3] * 792290124782L) + ((int128)tmp_q[4] * 8435347520905L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

