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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5277251372589387974UL) + ((((uint64_t)op[1] * 14950516434219846670UL) + ((uint64_t)op[2] * 12872402411383603099UL) + ((uint64_t)op[3] * 5241582139333608828UL) + ((uint64_t)op[4] * 335406096112389790UL) + ((uint64_t)op[5] * 16436183372595813519UL) + ((uint64_t)op[6] * 739611233332211369UL) + ((uint64_t)op[7] * 10436486552720950840UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 10436486552720950840UL) + ((uint64_t)op[1] * 5277251372589387974UL) + ((((uint64_t)op[2] * 14950516434219846670UL) + ((uint64_t)op[3] * 12872402411383603099UL) + ((uint64_t)op[4] * 5241582139333608828UL) + ((uint64_t)op[5] * 335406096112389790UL) + ((uint64_t)op[6] * 16436183372595813519UL) + ((uint64_t)op[7] * 739611233332211369UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 739611233332211369UL) + ((uint64_t)op[1] * 10436486552720950840UL) + ((uint64_t)op[2] * 5277251372589387974UL) + ((((uint64_t)op[3] * 14950516434219846670UL) + ((uint64_t)op[4] * 12872402411383603099UL) + ((uint64_t)op[5] * 5241582139333608828UL) + ((uint64_t)op[6] * 335406096112389790UL) + ((uint64_t)op[7] * 16436183372595813519UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 16436183372595813519UL) + ((uint64_t)op[1] * 739611233332211369UL) + ((uint64_t)op[2] * 10436486552720950840UL) + ((uint64_t)op[3] * 5277251372589387974UL) + ((((uint64_t)op[4] * 14950516434219846670UL) + ((uint64_t)op[5] * 12872402411383603099UL) + ((uint64_t)op[6] * 5241582139333608828UL) + ((uint64_t)op[7] * 335406096112389790UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 335406096112389790UL) + ((uint64_t)op[1] * 16436183372595813519UL) + ((uint64_t)op[2] * 739611233332211369UL) + ((uint64_t)op[3] * 10436486552720950840UL) + ((uint64_t)op[4] * 5277251372589387974UL) + ((((uint64_t)op[5] * 14950516434219846670UL) + ((uint64_t)op[6] * 12872402411383603099UL) + ((uint64_t)op[7] * 5241582139333608828UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 5241582139333608828UL) + ((uint64_t)op[1] * 335406096112389790UL) + ((uint64_t)op[2] * 16436183372595813519UL) + ((uint64_t)op[3] * 739611233332211369UL) + ((uint64_t)op[4] * 10436486552720950840UL) + ((uint64_t)op[5] * 5277251372589387974UL) + ((((uint64_t)op[6] * 14950516434219846670UL) + ((uint64_t)op[7] * 12872402411383603099UL)) * 18446744073709551609);
	tmp_q[6] = ((uint64_t)op[0] * 12872402411383603099UL) + ((uint64_t)op[1] * 5241582139333608828UL) + ((uint64_t)op[2] * 335406096112389790UL) + ((uint64_t)op[3] * 16436183372595813519UL) + ((uint64_t)op[4] * 739611233332211369UL) + ((uint64_t)op[5] * 10436486552720950840UL) + ((uint64_t)op[6] * 5277251372589387974UL) + ((uint64_t)op[7] * 6026849402718383006UL);
	tmp_q[7] = ((uint64_t)op[0] * 14950516434219846670UL) + ((uint64_t)op[1] * 12872402411383603099UL) + ((uint64_t)op[2] * 5241582139333608828UL) + ((uint64_t)op[3] * 335406096112389790UL) + ((uint64_t)op[4] * 16436183372595813519UL) + ((uint64_t)op[5] * 739611233332211369UL) + ((uint64_t)op[6] * 10436486552720950840UL) + ((uint64_t)op[7] * 5277251372589387974UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 55556900733627L) - ((((int128)tmp_q[1] * 170507189643102L) + ((int128)tmp_q[2] * 62519349263394L) - ((int128)tmp_q[3] * 109138255718643L) + ((int128)tmp_q[4] * 16037098347889L) - ((int128)tmp_q[5] * 76840802931618L) + ((int128)tmp_q[6] * 53828688025458L) + ((int128)tmp_q[7] * 125872103696950L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 125872103696950L) + ((int128)tmp_q[1] * 55556900733627L) - ((((int128)tmp_q[2] * 170507189643102L) + ((int128)tmp_q[3] * 62519349263394L) - ((int128)tmp_q[4] * 109138255718643L) + ((int128)tmp_q[5] * 16037098347889L) - ((int128)tmp_q[6] * 76840802931618L) + ((int128)tmp_q[7] * 53828688025458L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 53828688025458L) + ((int128)tmp_q[1] * 125872103696950L) + ((int128)tmp_q[2] * 55556900733627L) - ((((int128)tmp_q[3] * 170507189643102L) + ((int128)tmp_q[4] * 62519349263394L) - ((int128)tmp_q[5] * 109138255718643L) + ((int128)tmp_q[6] * 16037098347889L) - ((int128)tmp_q[7] * 76840802931618L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 76840802931618L) + ((int128)tmp_q[1] * 53828688025458L) + ((int128)tmp_q[2] * 125872103696950L) + ((int128)tmp_q[3] * 55556900733627L) - ((((int128)tmp_q[4] * 170507189643102L) + ((int128)tmp_q[5] * 62519349263394L) - ((int128)tmp_q[6] * 109138255718643L) + ((int128)tmp_q[7] * 16037098347889L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 16037098347889L) - ((int128)tmp_q[1] * 76840802931618L) + ((int128)tmp_q[2] * 53828688025458L) + ((int128)tmp_q[3] * 125872103696950L) + ((int128)tmp_q[4] * 55556900733627L) - ((((int128)tmp_q[5] * 170507189643102L) + ((int128)tmp_q[6] * 62519349263394L) - ((int128)tmp_q[7] * 109138255718643L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 109138255718643L) + ((int128)tmp_q[1] * 16037098347889L) - ((int128)tmp_q[2] * 76840802931618L) + ((int128)tmp_q[3] * 53828688025458L) + ((int128)tmp_q[4] * 125872103696950L) + ((int128)tmp_q[5] * 55556900733627L) - ((((int128)tmp_q[6] * 170507189643102L) + ((int128)tmp_q[7] * 62519349263394L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 62519349263394L) - ((int128)tmp_q[1] * 109138255718643L) + ((int128)tmp_q[2] * 16037098347889L) - ((int128)tmp_q[3] * 76840802931618L) + ((int128)tmp_q[4] * 53828688025458L) + ((int128)tmp_q[5] * 125872103696950L) + ((int128)tmp_q[6] * 55556900733627L) - ((int128)tmp_q[7] * 1193550327501714L);
	tmp_zero[7] = ((int128)tmp_q[0] * 170507189643102L) + ((int128)tmp_q[1] * 62519349263394L) - ((int128)tmp_q[2] * 109138255718643L) + ((int128)tmp_q[3] * 16037098347889L) - ((int128)tmp_q[4] * 76840802931618L) + ((int128)tmp_q[5] * 53828688025458L) + ((int128)tmp_q[6] * 125872103696950L) + ((int128)tmp_q[7] * 55556900733627L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

