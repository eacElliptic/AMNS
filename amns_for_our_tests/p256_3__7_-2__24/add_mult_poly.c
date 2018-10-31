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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1295668563987385287UL) + ((((uint64_t)op[1] * 2001709999650063250UL) + ((uint64_t)op[2] * 17247567604054878147UL) + ((uint64_t)op[3] * 18335961478096673915UL) + ((uint64_t)op[4] * 4467496868157037645UL) + ((uint64_t)op[5] * 16092105766705864488UL) + ((uint64_t)op[6] * 2278572629573953882UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 2278572629573953882UL) + ((uint64_t)op[1] * 1295668563987385287UL) + ((((uint64_t)op[2] * 2001709999650063250UL) + ((uint64_t)op[3] * 17247567604054878147UL) + ((uint64_t)op[4] * 18335961478096673915UL) + ((uint64_t)op[5] * 4467496868157037645UL) + ((uint64_t)op[6] * 16092105766705864488UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 16092105766705864488UL) + ((uint64_t)op[1] * 2278572629573953882UL) + ((uint64_t)op[2] * 1295668563987385287UL) + ((((uint64_t)op[3] * 2001709999650063250UL) + ((uint64_t)op[4] * 17247567604054878147UL) + ((uint64_t)op[5] * 18335961478096673915UL) + ((uint64_t)op[6] * 4467496868157037645UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 4467496868157037645UL) + ((uint64_t)op[1] * 16092105766705864488UL) + ((uint64_t)op[2] * 2278572629573953882UL) + ((uint64_t)op[3] * 1295668563987385287UL) + ((((uint64_t)op[4] * 2001709999650063250UL) + ((uint64_t)op[5] * 17247567604054878147UL) + ((uint64_t)op[6] * 18335961478096673915UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 18335961478096673915UL) + ((uint64_t)op[1] * 4467496868157037645UL) + ((uint64_t)op[2] * 16092105766705864488UL) + ((uint64_t)op[3] * 2278572629573953882UL) + ((uint64_t)op[4] * 1295668563987385287UL) + ((((uint64_t)op[5] * 2001709999650063250UL) + ((uint64_t)op[6] * 17247567604054878147UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 17247567604054878147UL) + ((uint64_t)op[1] * 18335961478096673915UL) + ((uint64_t)op[2] * 4467496868157037645UL) + ((uint64_t)op[3] * 16092105766705864488UL) + ((uint64_t)op[4] * 2278572629573953882UL) + ((uint64_t)op[5] * 1295668563987385287UL) + ((uint64_t)op[6] * 14443324074409425116UL);
	tmp_q[6] = ((uint64_t)op[0] * 2001709999650063250UL) + ((uint64_t)op[1] * 17247567604054878147UL) + ((uint64_t)op[2] * 18335961478096673915UL) + ((uint64_t)op[3] * 4467496868157037645UL) + ((uint64_t)op[4] * 16092105766705864488UL) + ((uint64_t)op[5] * 2278572629573953882UL) + ((uint64_t)op[6] * 1295668563987385287UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 32046951139L) - ((((int128)tmp_q[1] * 35303289771L) - ((int128)tmp_q[2] * 11013815481L) - ((int128)tmp_q[3] * 33269879689L) + ((int128)tmp_q[4] * 54246295789L) + ((int128)tmp_q[5] * 33396968058L) + ((int128)tmp_q[6] * 36894473308L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 36894473308L) - ((int128)tmp_q[1] * 32046951139L) - ((((int128)tmp_q[2] * 35303289771L) - ((int128)tmp_q[3] * 11013815481L) - ((int128)tmp_q[4] * 33269879689L) + ((int128)tmp_q[5] * 54246295789L) + ((int128)tmp_q[6] * 33396968058L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 33396968058L) + ((int128)tmp_q[1] * 36894473308L) - ((int128)tmp_q[2] * 32046951139L) - ((((int128)tmp_q[3] * 35303289771L) - ((int128)tmp_q[4] * 11013815481L) - ((int128)tmp_q[5] * 33269879689L) + ((int128)tmp_q[6] * 54246295789L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 54246295789L) + ((int128)tmp_q[1] * 33396968058L) + ((int128)tmp_q[2] * 36894473308L) - ((int128)tmp_q[3] * 32046951139L) - ((((int128)tmp_q[4] * 35303289771L) - ((int128)tmp_q[5] * 11013815481L) - ((int128)tmp_q[6] * 33269879689L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 33269879689L) + ((int128)tmp_q[1] * 54246295789L) + ((int128)tmp_q[2] * 33396968058L) + ((int128)tmp_q[3] * 36894473308L) - ((int128)tmp_q[4] * 32046951139L) - ((((int128)tmp_q[5] * 35303289771L) - ((int128)tmp_q[6] * 11013815481L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 11013815481L) - ((int128)tmp_q[1] * 33269879689L) + ((int128)tmp_q[2] * 54246295789L) + ((int128)tmp_q[3] * 33396968058L) + ((int128)tmp_q[4] * 36894473308L) - ((int128)tmp_q[5] * 32046951139L) - ((int128)tmp_q[6] * 70606579542L);
	tmp_zero[6] = ((int128)tmp_q[0] * 35303289771L) - ((int128)tmp_q[1] * 11013815481L) - ((int128)tmp_q[2] * 33269879689L) + ((int128)tmp_q[3] * 54246295789L) + ((int128)tmp_q[4] * 33396968058L) + ((int128)tmp_q[5] * 36894473308L) - ((int128)tmp_q[6] * 32046951139L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

