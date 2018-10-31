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
	tmp_q[0] = ((uint64_t)op[0] * 6080105471648560523UL) + ((((uint64_t)op[1] * 8925767035915050407UL) + ((uint64_t)op[2] * 8868931626286951923UL) + ((uint64_t)op[3] * 3474739280752361114UL) + ((uint64_t)op[4] * 3143752576029775864UL) + ((uint64_t)op[5] * 9565729858729384182UL) + ((uint64_t)op[6] * 2890823269783851972UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 2890823269783851972UL) + ((uint64_t)op[1] * 6080105471648560523UL) + ((((uint64_t)op[2] * 8925767035915050407UL) + ((uint64_t)op[3] * 8868931626286951923UL) + ((uint64_t)op[4] * 3474739280752361114UL) + ((uint64_t)op[5] * 3143752576029775864UL) + ((uint64_t)op[6] * 9565729858729384182UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 9565729858729384182UL) + ((uint64_t)op[1] * 2890823269783851972UL) + ((uint64_t)op[2] * 6080105471648560523UL) + ((((uint64_t)op[3] * 8925767035915050407UL) + ((uint64_t)op[4] * 8868931626286951923UL) + ((uint64_t)op[5] * 3474739280752361114UL) + ((uint64_t)op[6] * 3143752576029775864UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 3143752576029775864UL) + ((uint64_t)op[1] * 9565729858729384182UL) + ((uint64_t)op[2] * 2890823269783851972UL) + ((uint64_t)op[3] * 6080105471648560523UL) + ((((uint64_t)op[4] * 8925767035915050407UL) + ((uint64_t)op[5] * 8868931626286951923UL) + ((uint64_t)op[6] * 3474739280752361114UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 3474739280752361114UL) + ((uint64_t)op[1] * 3143752576029775864UL) + ((uint64_t)op[2] * 9565729858729384182UL) + ((uint64_t)op[3] * 2890823269783851972UL) + ((uint64_t)op[4] * 6080105471648560523UL) + ((((uint64_t)op[5] * 8925767035915050407UL) + ((uint64_t)op[6] * 8868931626286951923UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 8868931626286951923UL) + ((uint64_t)op[1] * 3474739280752361114UL) + ((uint64_t)op[2] * 3143752576029775864UL) + ((uint64_t)op[3] * 9565729858729384182UL) + ((uint64_t)op[4] * 2890823269783851972UL) + ((uint64_t)op[5] * 6080105471648560523UL) + ((uint64_t)op[6] * 10711397041553402813UL);
	tmp_q[6] = ((uint64_t)op[0] * 8925767035915050407UL) + ((uint64_t)op[1] * 8868931626286951923UL) + ((uint64_t)op[2] * 3474739280752361114UL) + ((uint64_t)op[3] * 3143752576029775864UL) + ((uint64_t)op[4] * 9565729858729384182UL) + ((uint64_t)op[5] * 2890823269783851972UL) + ((uint64_t)op[6] * 6080105471648560523UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 34866378257L) - ((-((int128)tmp_q[1] * 38298593424L) + ((int128)tmp_q[2] * 2508234167L) + ((int128)tmp_q[3] * 47862971783L) + ((int128)tmp_q[4] * 9924054204L) - ((int128)tmp_q[5] * 80211867321L) - ((int128)tmp_q[6] * 38238276935L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 38238276935L) - ((int128)tmp_q[1] * 34866378257L) - ((-((int128)tmp_q[2] * 38298593424L) + ((int128)tmp_q[3] * 2508234167L) + ((int128)tmp_q[4] * 47862971783L) + ((int128)tmp_q[5] * 9924054204L) - ((int128)tmp_q[6] * 80211867321L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 80211867321L) - ((int128)tmp_q[1] * 38238276935L) - ((int128)tmp_q[2] * 34866378257L) - ((-((int128)tmp_q[3] * 38298593424L) + ((int128)tmp_q[4] * 2508234167L) + ((int128)tmp_q[5] * 47862971783L) + ((int128)tmp_q[6] * 9924054204L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 9924054204L) - ((int128)tmp_q[1] * 80211867321L) - ((int128)tmp_q[2] * 38238276935L) - ((int128)tmp_q[3] * 34866378257L) - ((-((int128)tmp_q[4] * 38298593424L) + ((int128)tmp_q[5] * 2508234167L) + ((int128)tmp_q[6] * 47862971783L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 47862971783L) + ((int128)tmp_q[1] * 9924054204L) - ((int128)tmp_q[2] * 80211867321L) - ((int128)tmp_q[3] * 38238276935L) - ((int128)tmp_q[4] * 34866378257L) - ((-((int128)tmp_q[5] * 38298593424L) + ((int128)tmp_q[6] * 2508234167L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 2508234167L) + ((int128)tmp_q[1] * 47862971783L) + ((int128)tmp_q[2] * 9924054204L) - ((int128)tmp_q[3] * 80211867321L) - ((int128)tmp_q[4] * 38238276935L) - ((int128)tmp_q[5] * 34866378257L) + ((int128)tmp_q[6] * 191492967120L);
	tmp_zero[6] = -((int128)tmp_q[0] * 38298593424L) + ((int128)tmp_q[1] * 2508234167L) + ((int128)tmp_q[2] * 47862971783L) + ((int128)tmp_q[3] * 9924054204L) - ((int128)tmp_q[4] * 80211867321L) - ((int128)tmp_q[5] * 38238276935L) - ((int128)tmp_q[6] * 34866378257L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

