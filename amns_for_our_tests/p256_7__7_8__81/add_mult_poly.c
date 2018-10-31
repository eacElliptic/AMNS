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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 76867222858261677UL) + ((((uint64_t)op[1] * 16953905332306389088UL) + ((uint64_t)op[2] * 3256929337873934150UL) + ((uint64_t)op[3] * 11567135212537927370UL) + ((uint64_t)op[4] * 2451960023941413437UL) + ((uint64_t)op[5] * 15256621588651275518UL) + ((uint64_t)op[6] * 7993145877170752947UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 7993145877170752947UL) + ((uint64_t)op[1] * 76867222858261677UL) + ((((uint64_t)op[2] * 16953905332306389088UL) + ((uint64_t)op[3] * 3256929337873934150UL) + ((uint64_t)op[4] * 11567135212537927370UL) + ((uint64_t)op[5] * 2451960023941413437UL) + ((uint64_t)op[6] * 15256621588651275518UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 15256621588651275518UL) + ((uint64_t)op[1] * 7993145877170752947UL) + ((uint64_t)op[2] * 76867222858261677UL) + ((((uint64_t)op[3] * 16953905332306389088UL) + ((uint64_t)op[4] * 3256929337873934150UL) + ((uint64_t)op[5] * 11567135212537927370UL) + ((uint64_t)op[6] * 2451960023941413437UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 2451960023941413437UL) + ((uint64_t)op[1] * 15256621588651275518UL) + ((uint64_t)op[2] * 7993145877170752947UL) + ((uint64_t)op[3] * 76867222858261677UL) + ((((uint64_t)op[4] * 16953905332306389088UL) + ((uint64_t)op[5] * 3256929337873934150UL) + ((uint64_t)op[6] * 11567135212537927370UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 11567135212537927370UL) + ((uint64_t)op[1] * 2451960023941413437UL) + ((uint64_t)op[2] * 15256621588651275518UL) + ((uint64_t)op[3] * 7993145877170752947UL) + ((uint64_t)op[4] * 76867222858261677UL) + ((((uint64_t)op[5] * 16953905332306389088UL) + ((uint64_t)op[6] * 3256929337873934150UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 3256929337873934150UL) + ((uint64_t)op[1] * 11567135212537927370UL) + ((uint64_t)op[2] * 2451960023941413437UL) + ((uint64_t)op[3] * 15256621588651275518UL) + ((uint64_t)op[4] * 7993145877170752947UL) + ((uint64_t)op[5] * 76867222858261677UL) + ((uint64_t)op[6] * 6504034142484251392UL);
	tmp_q[6] = ((uint64_t)op[0] * 16953905332306389088UL) + ((uint64_t)op[1] * 3256929337873934150UL) + ((uint64_t)op[2] * 11567135212537927370UL) + ((uint64_t)op[3] * 2451960023941413437UL) + ((uint64_t)op[4] * 15256621588651275518UL) + ((uint64_t)op[5] * 7993145877170752947UL) + ((uint64_t)op[6] * 76867222858261677UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 17742008179L) + ((-((int128)tmp_q[1] * 44763679282L) + ((int128)tmp_q[2] * 53763738324L) + ((int128)tmp_q[3] * 18779700093L) + ((int128)tmp_q[4] * 8801072700L) + ((int128)tmp_q[5] * 39394791361L) + ((int128)tmp_q[6] * 10463708931L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 10463708931L) + ((int128)tmp_q[1] * 17742008179L) + ((-((int128)tmp_q[2] * 44763679282L) + ((int128)tmp_q[3] * 53763738324L) + ((int128)tmp_q[4] * 18779700093L) + ((int128)tmp_q[5] * 8801072700L) + ((int128)tmp_q[6] * 39394791361L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 39394791361L) + ((int128)tmp_q[1] * 10463708931L) + ((int128)tmp_q[2] * 17742008179L) + ((-((int128)tmp_q[3] * 44763679282L) + ((int128)tmp_q[4] * 53763738324L) + ((int128)tmp_q[5] * 18779700093L) + ((int128)tmp_q[6] * 8801072700L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 8801072700L) + ((int128)tmp_q[1] * 39394791361L) + ((int128)tmp_q[2] * 10463708931L) + ((int128)tmp_q[3] * 17742008179L) + ((-((int128)tmp_q[4] * 44763679282L) + ((int128)tmp_q[5] * 53763738324L) + ((int128)tmp_q[6] * 18779700093L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 18779700093L) + ((int128)tmp_q[1] * 8801072700L) + ((int128)tmp_q[2] * 39394791361L) + ((int128)tmp_q[3] * 10463708931L) + ((int128)tmp_q[4] * 17742008179L) + ((-((int128)tmp_q[5] * 44763679282L) + ((int128)tmp_q[6] * 53763738324L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 53763738324L) + ((int128)tmp_q[1] * 18779700093L) + ((int128)tmp_q[2] * 8801072700L) + ((int128)tmp_q[3] * 39394791361L) + ((int128)tmp_q[4] * 10463708931L) + ((int128)tmp_q[5] * 17742008179L) - ((int128)tmp_q[6] * 358109434256L);
	tmp_zero[6] = -((int128)tmp_q[0] * 44763679282L) + ((int128)tmp_q[1] * 53763738324L) + ((int128)tmp_q[2] * 18779700093L) + ((int128)tmp_q[3] * 8801072700L) + ((int128)tmp_q[4] * 39394791361L) + ((int128)tmp_q[5] * 10463708931L) + ((int128)tmp_q[6] * 17742008179L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

