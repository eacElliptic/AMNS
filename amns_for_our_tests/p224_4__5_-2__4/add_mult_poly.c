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
	tmp_q[0] = ((uint64_t)op[0] * 6433345580943596991UL) + ((((uint64_t)op[1] * 4789593761812032968UL) + ((uint64_t)op[2] * 11289283647169930277UL) + ((uint64_t)op[3] * 2134428048376384230UL) + ((uint64_t)op[4] * 5139377456636727075UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 5139377456636727075UL) + ((uint64_t)op[1] * 6433345580943596991UL) + ((((uint64_t)op[2] * 4789593761812032968UL) + ((uint64_t)op[3] * 11289283647169930277UL) + ((uint64_t)op[4] * 2134428048376384230UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2134428048376384230UL) + ((uint64_t)op[1] * 5139377456636727075UL) + ((uint64_t)op[2] * 6433345580943596991UL) + ((((uint64_t)op[3] * 4789593761812032968UL) + ((uint64_t)op[4] * 11289283647169930277UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 11289283647169930277UL) + ((uint64_t)op[1] * 2134428048376384230UL) + ((uint64_t)op[2] * 5139377456636727075UL) + ((uint64_t)op[3] * 6433345580943596991UL) + ((uint64_t)op[4] * 8867556550085485680UL);
	tmp_q[4] = ((uint64_t)op[0] * 4789593761812032968UL) + ((uint64_t)op[1] * 11289283647169930277UL) + ((uint64_t)op[2] * 2134428048376384230UL) + ((uint64_t)op[3] * 5139377456636727075UL) + ((uint64_t)op[4] * 6433345580943596991UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 19620935084147L) - ((-((int128)tmp_q[1] * 6081793782337L) + ((int128)tmp_q[2] * 3024238154282L) - ((int128)tmp_q[3] * 16267764223083L) - ((int128)tmp_q[4] * 15766073304225L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 15766073304225L) - ((int128)tmp_q[1] * 19620935084147L) - ((-((int128)tmp_q[2] * 6081793782337L) + ((int128)tmp_q[3] * 3024238154282L) - ((int128)tmp_q[4] * 16267764223083L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 16267764223083L) - ((int128)tmp_q[1] * 15766073304225L) - ((int128)tmp_q[2] * 19620935084147L) - ((-((int128)tmp_q[3] * 6081793782337L) + ((int128)tmp_q[4] * 3024238154282L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 3024238154282L) - ((int128)tmp_q[1] * 16267764223083L) - ((int128)tmp_q[2] * 15766073304225L) - ((int128)tmp_q[3] * 19620935084147L) + ((int128)tmp_q[4] * 12163587564674L);
	tmp_zero[4] = -((int128)tmp_q[0] * 6081793782337L) + ((int128)tmp_q[1] * 3024238154282L) - ((int128)tmp_q[2] * 16267764223083L) - ((int128)tmp_q[3] * 15766073304225L) - ((int128)tmp_q[4] * 19620935084147L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

