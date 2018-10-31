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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 164764989800463981UL) + ((((uint64_t)op[1] * 8014050892508533112UL) + ((uint64_t)op[2] * 2207732639735558995UL) + ((uint64_t)op[3] * 13311084321199601676UL) + ((uint64_t)op[4] * 17926323916765344497UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 17926323916765344497UL) + ((uint64_t)op[1] * 164764989800463981UL) + ((((uint64_t)op[2] * 8014050892508533112UL) + ((uint64_t)op[3] * 2207732639735558995UL) + ((uint64_t)op[4] * 13311084321199601676UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 13311084321199601676UL) + ((uint64_t)op[1] * 17926323916765344497UL) + ((uint64_t)op[2] * 164764989800463981UL) + ((((uint64_t)op[3] * 8014050892508533112UL) + ((uint64_t)op[4] * 2207732639735558995UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 2207732639735558995UL) + ((uint64_t)op[1] * 13311084321199601676UL) + ((uint64_t)op[2] * 17926323916765344497UL) + ((uint64_t)op[3] * 164764989800463981UL) + ((uint64_t)op[4] * 15269977758585989288UL);
	tmp_q[4] = ((uint64_t)op[0] * 8014050892508533112UL) + ((uint64_t)op[1] * 2207732639735558995UL) + ((uint64_t)op[2] * 13311084321199601676UL) + ((uint64_t)op[3] * 17926323916765344497UL) + ((uint64_t)op[4] * 164764989800463981UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 569599961708792L) - ((((int128)tmp_q[1] * 488705589303398L) + ((int128)tmp_q[2] * 920410023665511L) - ((int128)tmp_q[3] * 1445698921211727L) - ((int128)tmp_q[4] * 2130938985425835L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 2130938985425835L) - ((int128)tmp_q[1] * 569599961708792L) - ((((int128)tmp_q[2] * 488705589303398L) + ((int128)tmp_q[3] * 920410023665511L) - ((int128)tmp_q[4] * 1445698921211727L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 1445698921211727L) - ((int128)tmp_q[1] * 2130938985425835L) - ((int128)tmp_q[2] * 569599961708792L) - ((((int128)tmp_q[3] * 488705589303398L) + ((int128)tmp_q[4] * 920410023665511L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 920410023665511L) - ((int128)tmp_q[1] * 1445698921211727L) - ((int128)tmp_q[2] * 2130938985425835L) - ((int128)tmp_q[3] * 569599961708792L) - ((int128)tmp_q[4] * 2443527946516990L);
	tmp_zero[4] = ((int128)tmp_q[0] * 488705589303398L) + ((int128)tmp_q[1] * 920410023665511L) - ((int128)tmp_q[2] * 1445698921211727L) - ((int128)tmp_q[3] * 2130938985425835L) - ((int128)tmp_q[4] * 569599961708792L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

