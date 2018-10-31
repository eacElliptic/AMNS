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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13441460072836382991UL) + ((((uint64_t)op[1] * 4332765571443392199UL) + ((uint64_t)op[2] * 2670222829598437723UL) + ((uint64_t)op[3] * 15448466844864124447UL) + ((uint64_t)op[4] * 4449484157402785662UL) + ((uint64_t)op[5] * 16534772701008525406UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 16534772701008525406UL) + ((uint64_t)op[1] * 13441460072836382991UL) + ((((uint64_t)op[2] * 4332765571443392199UL) + ((uint64_t)op[3] * 2670222829598437723UL) + ((uint64_t)op[4] * 15448466844864124447UL) + ((uint64_t)op[5] * 4449484157402785662UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4449484157402785662UL) + ((uint64_t)op[1] * 16534772701008525406UL) + ((uint64_t)op[2] * 13441460072836382991UL) + ((((uint64_t)op[3] * 4332765571443392199UL) + ((uint64_t)op[4] * 2670222829598437723UL) + ((uint64_t)op[5] * 15448466844864124447UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 15448466844864124447UL) + ((uint64_t)op[1] * 4449484157402785662UL) + ((uint64_t)op[2] * 16534772701008525406UL) + ((uint64_t)op[3] * 13441460072836382991UL) + ((((uint64_t)op[4] * 4332765571443392199UL) + ((uint64_t)op[5] * 2670222829598437723UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 2670222829598437723UL) + ((uint64_t)op[1] * 15448466844864124447UL) + ((uint64_t)op[2] * 4449484157402785662UL) + ((uint64_t)op[3] * 16534772701008525406UL) + ((uint64_t)op[4] * 13441460072836382991UL) + ((uint64_t)op[5] * 9781212930822767218UL);
	tmp_q[5] = ((uint64_t)op[0] * 4332765571443392199UL) + ((uint64_t)op[1] * 2670222829598437723UL) + ((uint64_t)op[2] * 15448466844864124447UL) + ((uint64_t)op[3] * 4449484157402785662UL) + ((uint64_t)op[4] * 16534772701008525406UL) + ((uint64_t)op[5] * 13441460072836382991UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 94040772609L) - ((-((int128)tmp_q[1] * 35497602433L) - ((int128)tmp_q[2] * 121946650237L) - ((int128)tmp_q[3] * 2656576183L) + ((int128)tmp_q[4] * 61470068756L) - ((int128)tmp_q[5] * 69617768890L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 69617768890L) - ((int128)tmp_q[1] * 94040772609L) - ((-((int128)tmp_q[2] * 35497602433L) - ((int128)tmp_q[3] * 121946650237L) - ((int128)tmp_q[4] * 2656576183L) + ((int128)tmp_q[5] * 61470068756L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 61470068756L) - ((int128)tmp_q[1] * 69617768890L) - ((int128)tmp_q[2] * 94040772609L) - ((-((int128)tmp_q[3] * 35497602433L) - ((int128)tmp_q[4] * 121946650237L) - ((int128)tmp_q[5] * 2656576183L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2656576183L) + ((int128)tmp_q[1] * 61470068756L) - ((int128)tmp_q[2] * 69617768890L) - ((int128)tmp_q[3] * 94040772609L) - ((-((int128)tmp_q[4] * 35497602433L) - ((int128)tmp_q[5] * 121946650237L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 121946650237L) - ((int128)tmp_q[1] * 2656576183L) + ((int128)tmp_q[2] * 61470068756L) - ((int128)tmp_q[3] * 69617768890L) - ((int128)tmp_q[4] * 94040772609L) + ((int128)tmp_q[5] * 70995204866L);
	tmp_zero[5] = -((int128)tmp_q[0] * 35497602433L) - ((int128)tmp_q[1] * 121946650237L) - ((int128)tmp_q[2] * 2656576183L) + ((int128)tmp_q[3] * 61470068756L) - ((int128)tmp_q[4] * 69617768890L) - ((int128)tmp_q[5] * 94040772609L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

