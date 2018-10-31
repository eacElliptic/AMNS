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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - ((int128)pa[6] * pb[6]);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - ((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - ((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - ((int128)pa[6] * pa[6]);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14430079939325993909UL) + ((((uint64_t)op[1] * 16734810617739159705UL) + ((uint64_t)op[2] * 8754919553291770802UL) + ((uint64_t)op[3] * 1513090461220109411UL) + ((uint64_t)op[4] * 13076486581935683939UL) + ((uint64_t)op[5] * 14652896196379461299UL) + ((uint64_t)op[6] * 903812314491200762UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 903812314491200762UL) + ((uint64_t)op[1] * 14430079939325993909UL) + ((((uint64_t)op[2] * 16734810617739159705UL) + ((uint64_t)op[3] * 8754919553291770802UL) + ((uint64_t)op[4] * 1513090461220109411UL) + ((uint64_t)op[5] * 13076486581935683939UL) + ((uint64_t)op[6] * 14652896196379461299UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 14652896196379461299UL) + ((uint64_t)op[1] * 903812314491200762UL) + ((uint64_t)op[2] * 14430079939325993909UL) + ((((uint64_t)op[3] * 16734810617739159705UL) + ((uint64_t)op[4] * 8754919553291770802UL) + ((uint64_t)op[5] * 1513090461220109411UL) + ((uint64_t)op[6] * 13076486581935683939UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 13076486581935683939UL) + ((uint64_t)op[1] * 14652896196379461299UL) + ((uint64_t)op[2] * 903812314491200762UL) + ((uint64_t)op[3] * 14430079939325993909UL) + ((((uint64_t)op[4] * 16734810617739159705UL) + ((uint64_t)op[5] * 8754919553291770802UL) + ((uint64_t)op[6] * 1513090461220109411UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 1513090461220109411UL) + ((uint64_t)op[1] * 13076486581935683939UL) + ((uint64_t)op[2] * 14652896196379461299UL) + ((uint64_t)op[3] * 903812314491200762UL) + ((uint64_t)op[4] * 14430079939325993909UL) + ((((uint64_t)op[5] * 16734810617739159705UL) + ((uint64_t)op[6] * 8754919553291770802UL)) * 18446744073709551615);
	tmp_q[5] = ((uint64_t)op[0] * 8754919553291770802UL) + ((uint64_t)op[1] * 1513090461220109411UL) + ((uint64_t)op[2] * 13076486581935683939UL) + ((uint64_t)op[3] * 14652896196379461299UL) + ((uint64_t)op[4] * 903812314491200762UL) + ((uint64_t)op[5] * 14430079939325993909UL) + ((uint64_t)op[6] * 1711933455970391911UL);
	tmp_q[6] = ((uint64_t)op[0] * 16734810617739159705UL) + ((uint64_t)op[1] * 8754919553291770802UL) + ((uint64_t)op[2] * 1513090461220109411UL) + ((uint64_t)op[3] * 13076486581935683939UL) + ((uint64_t)op[4] * 14652896196379461299UL) + ((uint64_t)op[5] * 903812314491200762UL) + ((uint64_t)op[6] * 14430079939325993909UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2101002494198L) - (-((int128)tmp_q[1] * 258617172668L) - ((int128)tmp_q[2] * 4372318292665L) + ((int128)tmp_q[3] * 32033487321L) + ((int128)tmp_q[4] * 729427631595L) - ((int128)tmp_q[5] * 4666366794184L) + ((int128)tmp_q[6] * 850942675734L));
	tmp_zero[1] = ((int128)tmp_q[0] * 850942675734L) + ((int128)tmp_q[1] * 2101002494198L) - (-((int128)tmp_q[2] * 258617172668L) - ((int128)tmp_q[3] * 4372318292665L) + ((int128)tmp_q[4] * 32033487321L) + ((int128)tmp_q[5] * 729427631595L) - ((int128)tmp_q[6] * 4666366794184L));
	tmp_zero[2] = -((int128)tmp_q[0] * 4666366794184L) + ((int128)tmp_q[1] * 850942675734L) + ((int128)tmp_q[2] * 2101002494198L) - (-((int128)tmp_q[3] * 258617172668L) - ((int128)tmp_q[4] * 4372318292665L) + ((int128)tmp_q[5] * 32033487321L) + ((int128)tmp_q[6] * 729427631595L));
	tmp_zero[3] = ((int128)tmp_q[0] * 729427631595L) - ((int128)tmp_q[1] * 4666366794184L) + ((int128)tmp_q[2] * 850942675734L) + ((int128)tmp_q[3] * 2101002494198L) - (-((int128)tmp_q[4] * 258617172668L) - ((int128)tmp_q[5] * 4372318292665L) + ((int128)tmp_q[6] * 32033487321L));
	tmp_zero[4] = ((int128)tmp_q[0] * 32033487321L) + ((int128)tmp_q[1] * 729427631595L) - ((int128)tmp_q[2] * 4666366794184L) + ((int128)tmp_q[3] * 850942675734L) + ((int128)tmp_q[4] * 2101002494198L) - (-((int128)tmp_q[5] * 258617172668L) - ((int128)tmp_q[6] * 4372318292665L));
	tmp_zero[5] = -((int128)tmp_q[0] * 4372318292665L) + ((int128)tmp_q[1] * 32033487321L) + ((int128)tmp_q[2] * 729427631595L) - ((int128)tmp_q[3] * 4666366794184L) + ((int128)tmp_q[4] * 850942675734L) + ((int128)tmp_q[5] * 2101002494198L) + ((int128)tmp_q[6] * 258617172668L);
	tmp_zero[6] = -((int128)tmp_q[0] * 258617172668L) - ((int128)tmp_q[1] * 4372318292665L) + ((int128)tmp_q[2] * 32033487321L) + ((int128)tmp_q[3] * 729427631595L) - ((int128)tmp_q[4] * 4666366794184L) + ((int128)tmp_q[5] * 850942675734L) + ((int128)tmp_q[6] * 2101002494198L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

