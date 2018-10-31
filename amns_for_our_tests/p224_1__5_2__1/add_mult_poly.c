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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16553040569207393595UL) + ((((uint64_t)op[1] * 5795943044003542573UL) + ((uint64_t)op[2] * 5348521024396878410UL) + ((uint64_t)op[3] * 558170535149685429UL) + ((uint64_t)op[4] * 16496634865050854322UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 16496634865050854322UL) + ((uint64_t)op[1] * 16553040569207393595UL) + ((((uint64_t)op[2] * 5795943044003542573UL) + ((uint64_t)op[3] * 5348521024396878410UL) + ((uint64_t)op[4] * 558170535149685429UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 558170535149685429UL) + ((uint64_t)op[1] * 16496634865050854322UL) + ((uint64_t)op[2] * 16553040569207393595UL) + ((((uint64_t)op[3] * 5795943044003542573UL) + ((uint64_t)op[4] * 5348521024396878410UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 5348521024396878410UL) + ((uint64_t)op[1] * 558170535149685429UL) + ((uint64_t)op[2] * 16496634865050854322UL) + ((uint64_t)op[3] * 16553040569207393595UL) + ((uint64_t)op[4] * 11591886088007085146UL);
	tmp_q[4] = ((uint64_t)op[0] * 5795943044003542573UL) + ((uint64_t)op[1] * 5348521024396878410UL) + ((uint64_t)op[2] * 558170535149685429UL) + ((uint64_t)op[3] * 16496634865050854322UL) + ((uint64_t)op[4] * 16553040569207393595UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4408837444281L) + ((((int128)tmp_q[1] * 14286763516258L) - ((int128)tmp_q[2] * 7009225068944L) + ((int128)tmp_q[3] * 393624550705L) + ((int128)tmp_q[4] * 8179333441136L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 8179333441136L) + ((int128)tmp_q[1] * 4408837444281L) + ((((int128)tmp_q[2] * 14286763516258L) - ((int128)tmp_q[3] * 7009225068944L) + ((int128)tmp_q[4] * 393624550705L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 393624550705L) + ((int128)tmp_q[1] * 8179333441136L) + ((int128)tmp_q[2] * 4408837444281L) + ((((int128)tmp_q[3] * 14286763516258L) - ((int128)tmp_q[4] * 7009225068944L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 7009225068944L) + ((int128)tmp_q[1] * 393624550705L) + ((int128)tmp_q[2] * 8179333441136L) + ((int128)tmp_q[3] * 4408837444281L) + ((int128)tmp_q[4] * 28573527032516L);
	tmp_zero[4] = ((int128)tmp_q[0] * 14286763516258L) - ((int128)tmp_q[1] * 7009225068944L) + ((int128)tmp_q[2] * 393624550705L) + ((int128)tmp_q[3] * 8179333441136L) + ((int128)tmp_q[4] * 4408837444281L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

