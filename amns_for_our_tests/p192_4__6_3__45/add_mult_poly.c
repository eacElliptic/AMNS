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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2977767978344375888UL) + ((((uint64_t)op[1] * 8486492614343142867UL) + ((uint64_t)op[2] * 5277847773540846226UL) + ((uint64_t)op[3] * 9972061412922412937UL) + ((uint64_t)op[4] * 13893625031020077345UL) + ((uint64_t)op[5] * 4683231677061486880UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 4683231677061486880UL) + ((uint64_t)op[1] * 2977767978344375888UL) + ((((uint64_t)op[2] * 8486492614343142867UL) + ((uint64_t)op[3] * 5277847773540846226UL) + ((uint64_t)op[4] * 9972061412922412937UL) + ((uint64_t)op[5] * 13893625031020077345UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 13893625031020077345UL) + ((uint64_t)op[1] * 4683231677061486880UL) + ((uint64_t)op[2] * 2977767978344375888UL) + ((((uint64_t)op[3] * 8486492614343142867UL) + ((uint64_t)op[4] * 5277847773540846226UL) + ((uint64_t)op[5] * 9972061412922412937UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 9972061412922412937UL) + ((uint64_t)op[1] * 13893625031020077345UL) + ((uint64_t)op[2] * 4683231677061486880UL) + ((uint64_t)op[3] * 2977767978344375888UL) + ((((uint64_t)op[4] * 8486492614343142867UL) + ((uint64_t)op[5] * 5277847773540846226UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 5277847773540846226UL) + ((uint64_t)op[1] * 9972061412922412937UL) + ((uint64_t)op[2] * 13893625031020077345UL) + ((uint64_t)op[3] * 4683231677061486880UL) + ((uint64_t)op[4] * 2977767978344375888UL) + ((uint64_t)op[5] * 7012733769319876985UL);
	tmp_q[5] = ((uint64_t)op[0] * 8486492614343142867UL) + ((uint64_t)op[1] * 5277847773540846226UL) + ((uint64_t)op[2] * 9972061412922412937UL) + ((uint64_t)op[3] * 13893625031020077345UL) + ((uint64_t)op[4] * 4683231677061486880UL) + ((uint64_t)op[5] * 2977767978344375888UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 231694910L) + ((((int128)tmp_q[1] * 1358597917L) - ((int128)tmp_q[2] * 449829740L) - ((int128)tmp_q[3] * 1862513729L) + ((int128)tmp_q[4] * 2110984123L) + ((int128)tmp_q[5] * 215689228L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 215689228L) - ((int128)tmp_q[1] * 231694910L) + ((((int128)tmp_q[2] * 1358597917L) - ((int128)tmp_q[3] * 449829740L) - ((int128)tmp_q[4] * 1862513729L) + ((int128)tmp_q[5] * 2110984123L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2110984123L) + ((int128)tmp_q[1] * 215689228L) - ((int128)tmp_q[2] * 231694910L) + ((((int128)tmp_q[3] * 1358597917L) - ((int128)tmp_q[4] * 449829740L) - ((int128)tmp_q[5] * 1862513729L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 1862513729L) + ((int128)tmp_q[1] * 2110984123L) + ((int128)tmp_q[2] * 215689228L) - ((int128)tmp_q[3] * 231694910L) + ((((int128)tmp_q[4] * 1358597917L) - ((int128)tmp_q[5] * 449829740L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 449829740L) - ((int128)tmp_q[1] * 1862513729L) + ((int128)tmp_q[2] * 2110984123L) + ((int128)tmp_q[3] * 215689228L) - ((int128)tmp_q[4] * 231694910L) + ((int128)tmp_q[5] * 4075793751L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1358597917L) - ((int128)tmp_q[1] * 449829740L) - ((int128)tmp_q[2] * 1862513729L) + ((int128)tmp_q[3] * 2110984123L) + ((int128)tmp_q[4] * 215689228L) - ((int128)tmp_q[5] * 231694910L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

