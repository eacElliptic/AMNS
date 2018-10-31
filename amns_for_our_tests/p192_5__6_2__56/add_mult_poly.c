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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6604146163712459249UL) + ((((uint64_t)op[1] * 11838216228545338313UL) + ((uint64_t)op[2] * 1759339048569902927UL) + ((uint64_t)op[3] * 14019204304259730047UL) + ((uint64_t)op[4] * 4124543574042569656UL) + ((uint64_t)op[5] * 6737052871188609364UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 6737052871188609364UL) + ((uint64_t)op[1] * 6604146163712459249UL) + ((((uint64_t)op[2] * 11838216228545338313UL) + ((uint64_t)op[3] * 1759339048569902927UL) + ((uint64_t)op[4] * 14019204304259730047UL) + ((uint64_t)op[5] * 4124543574042569656UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 4124543574042569656UL) + ((uint64_t)op[1] * 6737052871188609364UL) + ((uint64_t)op[2] * 6604146163712459249UL) + ((((uint64_t)op[3] * 11838216228545338313UL) + ((uint64_t)op[4] * 1759339048569902927UL) + ((uint64_t)op[5] * 14019204304259730047UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 14019204304259730047UL) + ((uint64_t)op[1] * 4124543574042569656UL) + ((uint64_t)op[2] * 6737052871188609364UL) + ((uint64_t)op[3] * 6604146163712459249UL) + ((((uint64_t)op[4] * 11838216228545338313UL) + ((uint64_t)op[5] * 1759339048569902927UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 1759339048569902927UL) + ((uint64_t)op[1] * 14019204304259730047UL) + ((uint64_t)op[2] * 4124543574042569656UL) + ((uint64_t)op[3] * 6737052871188609364UL) + ((uint64_t)op[4] * 6604146163712459249UL) + ((uint64_t)op[5] * 5229688383381125010UL);
	tmp_q[5] = ((uint64_t)op[0] * 11838216228545338313UL) + ((uint64_t)op[1] * 1759339048569902927UL) + ((uint64_t)op[2] * 14019204304259730047UL) + ((uint64_t)op[3] * 4124543574042569656UL) + ((uint64_t)op[4] * 6737052871188609364UL) + ((uint64_t)op[5] * 6604146163712459249UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1396397019L) + ((-((int128)tmp_q[1] * 2095815811L) + ((int128)tmp_q[2] * 213351623L) + ((int128)tmp_q[3] * 1865344761L) + ((int128)tmp_q[4] * 144122158L) - ((int128)tmp_q[5] * 1023713584L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1023713584L) - ((int128)tmp_q[1] * 1396397019L) + ((-((int128)tmp_q[2] * 2095815811L) + ((int128)tmp_q[3] * 213351623L) + ((int128)tmp_q[4] * 1865344761L) + ((int128)tmp_q[5] * 144122158L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 144122158L) - ((int128)tmp_q[1] * 1023713584L) - ((int128)tmp_q[2] * 1396397019L) + ((-((int128)tmp_q[3] * 2095815811L) + ((int128)tmp_q[4] * 213351623L) + ((int128)tmp_q[5] * 1865344761L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 1865344761L) + ((int128)tmp_q[1] * 144122158L) - ((int128)tmp_q[2] * 1023713584L) - ((int128)tmp_q[3] * 1396397019L) + ((-((int128)tmp_q[4] * 2095815811L) + ((int128)tmp_q[5] * 213351623L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 213351623L) + ((int128)tmp_q[1] * 1865344761L) + ((int128)tmp_q[2] * 144122158L) - ((int128)tmp_q[3] * 1023713584L) - ((int128)tmp_q[4] * 1396397019L) - ((int128)tmp_q[5] * 4191631622L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2095815811L) + ((int128)tmp_q[1] * 213351623L) + ((int128)tmp_q[2] * 1865344761L) + ((int128)tmp_q[3] * 144122158L) - ((int128)tmp_q[4] * 1023713584L) - ((int128)tmp_q[5] * 1396397019L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

