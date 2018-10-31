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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7695298333509096396UL) + ((((uint64_t)op[1] * 1902270364079965663UL) + ((uint64_t)op[2] * 18118976333733900773UL) + ((uint64_t)op[3] * 52106244904864951UL) + ((uint64_t)op[4] * 13726472851460177869UL) + ((uint64_t)op[5] * 14127763428685501817UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 14127763428685501817UL) + ((uint64_t)op[1] * 7695298333509096396UL) + ((((uint64_t)op[2] * 1902270364079965663UL) + ((uint64_t)op[3] * 18118976333733900773UL) + ((uint64_t)op[4] * 52106244904864951UL) + ((uint64_t)op[5] * 13726472851460177869UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 13726472851460177869UL) + ((uint64_t)op[1] * 14127763428685501817UL) + ((uint64_t)op[2] * 7695298333509096396UL) + ((((uint64_t)op[3] * 1902270364079965663UL) + ((uint64_t)op[4] * 18118976333733900773UL) + ((uint64_t)op[5] * 52106244904864951UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 52106244904864951UL) + ((uint64_t)op[1] * 13726472851460177869UL) + ((uint64_t)op[2] * 14127763428685501817UL) + ((uint64_t)op[3] * 7695298333509096396UL) + ((((uint64_t)op[4] * 1902270364079965663UL) + ((uint64_t)op[5] * 18118976333733900773UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 18118976333733900773UL) + ((uint64_t)op[1] * 52106244904864951UL) + ((uint64_t)op[2] * 13726472851460177869UL) + ((uint64_t)op[3] * 14127763428685501817UL) + ((uint64_t)op[4] * 7695298333509096396UL) + ((uint64_t)op[5] * 13315892548559759641UL);
	tmp_q[5] = ((uint64_t)op[0] * 1902270364079965663UL) + ((uint64_t)op[1] * 18118976333733900773UL) + ((uint64_t)op[2] * 52106244904864951UL) + ((uint64_t)op[3] * 13726472851460177869UL) + ((uint64_t)op[4] * 14127763428685501817UL) + ((uint64_t)op[5] * 7695298333509096396UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2535152880L) + ((-((int128)tmp_q[1] * 1217572737L) + ((int128)tmp_q[2] * 2185586839L) - ((int128)tmp_q[3] * 2134517401L) + ((int128)tmp_q[4] * 1240199851L) + ((int128)tmp_q[5] * 2177203297L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 2177203297L) - ((int128)tmp_q[1] * 2535152880L) + ((-((int128)tmp_q[2] * 1217572737L) + ((int128)tmp_q[3] * 2185586839L) - ((int128)tmp_q[4] * 2134517401L) + ((int128)tmp_q[5] * 1240199851L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 1240199851L) + ((int128)tmp_q[1] * 2177203297L) - ((int128)tmp_q[2] * 2535152880L) + ((-((int128)tmp_q[3] * 1217572737L) + ((int128)tmp_q[4] * 2185586839L) - ((int128)tmp_q[5] * 2134517401L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 2134517401L) + ((int128)tmp_q[1] * 1240199851L) + ((int128)tmp_q[2] * 2177203297L) - ((int128)tmp_q[3] * 2535152880L) + ((-((int128)tmp_q[4] * 1217572737L) + ((int128)tmp_q[5] * 2185586839L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 2185586839L) - ((int128)tmp_q[1] * 2134517401L) + ((int128)tmp_q[2] * 1240199851L) + ((int128)tmp_q[3] * 2177203297L) - ((int128)tmp_q[4] * 2535152880L) - ((int128)tmp_q[5] * 8523009159L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1217572737L) + ((int128)tmp_q[1] * 2185586839L) - ((int128)tmp_q[2] * 2134517401L) + ((int128)tmp_q[3] * 1240199851L) + ((int128)tmp_q[4] * 2177203297L) - ((int128)tmp_q[5] * 2535152880L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

