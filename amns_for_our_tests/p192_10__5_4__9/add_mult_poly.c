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
	tmp_q[0] = ((uint64_t)op[0] * 3429347481371803953UL) + ((((uint64_t)op[1] * 9941048980067155984UL) + ((uint64_t)op[2] * 320716702727001402UL) + ((uint64_t)op[3] * 7698391259374352424UL) + ((uint64_t)op[4] * 4683281208283613430UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 4683281208283613430UL) + ((uint64_t)op[1] * 3429347481371803953UL) + ((((uint64_t)op[2] * 9941048980067155984UL) + ((uint64_t)op[3] * 320716702727001402UL) + ((uint64_t)op[4] * 7698391259374352424UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 7698391259374352424UL) + ((uint64_t)op[1] * 4683281208283613430UL) + ((uint64_t)op[2] * 3429347481371803953UL) + ((((uint64_t)op[3] * 9941048980067155984UL) + ((uint64_t)op[4] * 320716702727001402UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 320716702727001402UL) + ((uint64_t)op[1] * 7698391259374352424UL) + ((uint64_t)op[2] * 4683281208283613430UL) + ((uint64_t)op[3] * 3429347481371803953UL) + ((uint64_t)op[4] * 2870707772849520704UL);
	tmp_q[4] = ((uint64_t)op[0] * 9941048980067155984UL) + ((uint64_t)op[1] * 320716702727001402UL) + ((uint64_t)op[2] * 7698391259374352424UL) + ((uint64_t)op[3] * 4683281208283613430UL) + ((uint64_t)op[4] * 3429347481371803953UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 48130585585L) + ((-((int128)tmp_q[1] * 146850294712L) - ((int128)tmp_q[2] * 48967426894L) + ((int128)tmp_q[3] * 158801773220L) - ((int128)tmp_q[4] * 117637458650L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 117637458650L) - ((int128)tmp_q[1] * 48130585585L) + ((-((int128)tmp_q[2] * 146850294712L) - ((int128)tmp_q[3] * 48967426894L) + ((int128)tmp_q[4] * 158801773220L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 158801773220L) - ((int128)tmp_q[1] * 117637458650L) - ((int128)tmp_q[2] * 48130585585L) + ((-((int128)tmp_q[3] * 146850294712L) - ((int128)tmp_q[4] * 48967426894L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 48967426894L) + ((int128)tmp_q[1] * 158801773220L) - ((int128)tmp_q[2] * 117637458650L) - ((int128)tmp_q[3] * 48130585585L) - ((int128)tmp_q[4] * 587401178848L);
	tmp_zero[4] = -((int128)tmp_q[0] * 146850294712L) - ((int128)tmp_q[1] * 48967426894L) + ((int128)tmp_q[2] * 158801773220L) - ((int128)tmp_q[3] * 117637458650L) - ((int128)tmp_q[4] * 48130585585L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

