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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2188859952055661593UL) + ((((uint64_t)op[1] * 8422923970540486684UL) + ((uint64_t)op[2] * 15343533475014245180UL) + ((uint64_t)op[3] * 17870504161271606307UL) + ((uint64_t)op[4] * 2121024367687498883UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 2121024367687498883UL) + ((uint64_t)op[1] * 2188859952055661593UL) + ((((uint64_t)op[2] * 8422923970540486684UL) + ((uint64_t)op[3] * 15343533475014245180UL) + ((uint64_t)op[4] * 17870504161271606307UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 17870504161271606307UL) + ((uint64_t)op[1] * 2121024367687498883UL) + ((uint64_t)op[2] * 2188859952055661593UL) + ((((uint64_t)op[3] * 8422923970540486684UL) + ((uint64_t)op[4] * 15343533475014245180UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 15343533475014245180UL) + ((uint64_t)op[1] * 17870504161271606307UL) + ((uint64_t)op[2] * 2121024367687498883UL) + ((uint64_t)op[3] * 2188859952055661593UL) + ((uint64_t)op[4] * 1600896132628578248UL);
	tmp_q[4] = ((uint64_t)op[0] * 8422923970540486684UL) + ((uint64_t)op[1] * 15343533475014245180UL) + ((uint64_t)op[2] * 17870504161271606307UL) + ((uint64_t)op[3] * 2121024367687498883UL) + ((uint64_t)op[4] * 2188859952055661593UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 485682794454013L) - ((((int128)tmp_q[1] * 1309035339947865L) + ((int128)tmp_q[2] * 883839629407877L) + ((int128)tmp_q[3] * 756266737590560L) - ((int128)tmp_q[4] * 1218172427169299L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1218172427169299L) - ((int128)tmp_q[1] * 485682794454013L) - ((((int128)tmp_q[2] * 1309035339947865L) + ((int128)tmp_q[3] * 883839629407877L) + ((int128)tmp_q[4] * 756266737590560L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 756266737590560L) - ((int128)tmp_q[1] * 1218172427169299L) - ((int128)tmp_q[2] * 485682794454013L) - ((((int128)tmp_q[3] * 1309035339947865L) + ((int128)tmp_q[4] * 883839629407877L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 883839629407877L) + ((int128)tmp_q[1] * 756266737590560L) - ((int128)tmp_q[2] * 1218172427169299L) - ((int128)tmp_q[3] * 485682794454013L) - ((int128)tmp_q[4] * 2618070679895730L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1309035339947865L) + ((int128)tmp_q[1] * 883839629407877L) + ((int128)tmp_q[2] * 756266737590560L) - ((int128)tmp_q[3] * 1218172427169299L) - ((int128)tmp_q[4] * 485682794454013L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

