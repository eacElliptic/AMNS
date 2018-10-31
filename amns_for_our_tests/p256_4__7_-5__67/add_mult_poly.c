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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3020856079824767406UL) + ((((uint64_t)op[1] * 6820716610515491968UL) + ((uint64_t)op[2] * 15333838507835274135UL) + ((uint64_t)op[3] * 15119700646151992782UL) + ((uint64_t)op[4] * 9079245494297303832UL) + ((uint64_t)op[5] * 536682177756858561UL) + ((uint64_t)op[6] * 5131166852086022705UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 5131166852086022705UL) + ((uint64_t)op[1] * 3020856079824767406UL) + ((((uint64_t)op[2] * 6820716610515491968UL) + ((uint64_t)op[3] * 15333838507835274135UL) + ((uint64_t)op[4] * 15119700646151992782UL) + ((uint64_t)op[5] * 9079245494297303832UL) + ((uint64_t)op[6] * 536682177756858561UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 536682177756858561UL) + ((uint64_t)op[1] * 5131166852086022705UL) + ((uint64_t)op[2] * 3020856079824767406UL) + ((((uint64_t)op[3] * 6820716610515491968UL) + ((uint64_t)op[4] * 15333838507835274135UL) + ((uint64_t)op[5] * 15119700646151992782UL) + ((uint64_t)op[6] * 9079245494297303832UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 9079245494297303832UL) + ((uint64_t)op[1] * 536682177756858561UL) + ((uint64_t)op[2] * 5131166852086022705UL) + ((uint64_t)op[3] * 3020856079824767406UL) + ((((uint64_t)op[4] * 6820716610515491968UL) + ((uint64_t)op[5] * 15333838507835274135UL) + ((uint64_t)op[6] * 15119700646151992782UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 15119700646151992782UL) + ((uint64_t)op[1] * 9079245494297303832UL) + ((uint64_t)op[2] * 536682177756858561UL) + ((uint64_t)op[3] * 5131166852086022705UL) + ((uint64_t)op[4] * 3020856079824767406UL) + ((((uint64_t)op[5] * 6820716610515491968UL) + ((uint64_t)op[6] * 15333838507835274135UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 15333838507835274135UL) + ((uint64_t)op[1] * 15119700646151992782UL) + ((uint64_t)op[2] * 9079245494297303832UL) + ((uint64_t)op[3] * 536682177756858561UL) + ((uint64_t)op[4] * 5131166852086022705UL) + ((uint64_t)op[5] * 3020856079824767406UL) + ((uint64_t)op[6] * 2789905094841643392UL);
	tmp_q[6] = ((uint64_t)op[0] * 6820716610515491968UL) + ((uint64_t)op[1] * 15333838507835274135UL) + ((uint64_t)op[2] * 15119700646151992782UL) + ((uint64_t)op[3] * 9079245494297303832UL) + ((uint64_t)op[4] * 536682177756858561UL) + ((uint64_t)op[5] * 5131166852086022705UL) + ((uint64_t)op[6] * 3020856079824767406UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 77078201659L) - ((-((int128)tmp_q[1] * 41215576737L) - ((int128)tmp_q[2] * 31446446983L) - ((int128)tmp_q[3] * 4152906129L) + ((int128)tmp_q[4] * 16194329624L) + ((int128)tmp_q[5] * 24457339485L) - ((int128)tmp_q[6] * 38150576896L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 38150576896L) - ((int128)tmp_q[1] * 77078201659L) - ((-((int128)tmp_q[2] * 41215576737L) - ((int128)tmp_q[3] * 31446446983L) - ((int128)tmp_q[4] * 4152906129L) + ((int128)tmp_q[5] * 16194329624L) + ((int128)tmp_q[6] * 24457339485L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 24457339485L) - ((int128)tmp_q[1] * 38150576896L) - ((int128)tmp_q[2] * 77078201659L) - ((-((int128)tmp_q[3] * 41215576737L) - ((int128)tmp_q[4] * 31446446983L) - ((int128)tmp_q[5] * 4152906129L) + ((int128)tmp_q[6] * 16194329624L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 16194329624L) + ((int128)tmp_q[1] * 24457339485L) - ((int128)tmp_q[2] * 38150576896L) - ((int128)tmp_q[3] * 77078201659L) - ((-((int128)tmp_q[4] * 41215576737L) - ((int128)tmp_q[5] * 31446446983L) - ((int128)tmp_q[6] * 4152906129L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 4152906129L) + ((int128)tmp_q[1] * 16194329624L) + ((int128)tmp_q[2] * 24457339485L) - ((int128)tmp_q[3] * 38150576896L) - ((int128)tmp_q[4] * 77078201659L) - ((-((int128)tmp_q[5] * 41215576737L) - ((int128)tmp_q[6] * 31446446983L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 31446446983L) - ((int128)tmp_q[1] * 4152906129L) + ((int128)tmp_q[2] * 16194329624L) + ((int128)tmp_q[3] * 24457339485L) - ((int128)tmp_q[4] * 38150576896L) - ((int128)tmp_q[5] * 77078201659L) + ((int128)tmp_q[6] * 206077883685L);
	tmp_zero[6] = -((int128)tmp_q[0] * 41215576737L) - ((int128)tmp_q[1] * 31446446983L) - ((int128)tmp_q[2] * 4152906129L) + ((int128)tmp_q[3] * 16194329624L) + ((int128)tmp_q[4] * 24457339485L) - ((int128)tmp_q[5] * 38150576896L) - ((int128)tmp_q[6] * 77078201659L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

