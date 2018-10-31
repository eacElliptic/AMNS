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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8421140006141382512UL) + ((((uint64_t)op[1] * 1768846766534990688UL) + ((uint64_t)op[2] * 6381366926573354045UL) + ((uint64_t)op[3] * 15049723021932632948UL) + ((uint64_t)op[4] * 16375509950581339924UL) + ((uint64_t)op[5] * 10353042788260867254UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 10353042788260867254UL) + ((uint64_t)op[1] * 8421140006141382512UL) + ((((uint64_t)op[2] * 1768846766534990688UL) + ((uint64_t)op[3] * 6381366926573354045UL) + ((uint64_t)op[4] * 15049723021932632948UL) + ((uint64_t)op[5] * 16375509950581339924UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 16375509950581339924UL) + ((uint64_t)op[1] * 10353042788260867254UL) + ((uint64_t)op[2] * 8421140006141382512UL) + ((((uint64_t)op[3] * 1768846766534990688UL) + ((uint64_t)op[4] * 6381366926573354045UL) + ((uint64_t)op[5] * 15049723021932632948UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 15049723021932632948UL) + ((uint64_t)op[1] * 16375509950581339924UL) + ((uint64_t)op[2] * 10353042788260867254UL) + ((uint64_t)op[3] * 8421140006141382512UL) + ((((uint64_t)op[4] * 1768846766534990688UL) + ((uint64_t)op[5] * 6381366926573354045UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 6381366926573354045UL) + ((uint64_t)op[1] * 15049723021932632948UL) + ((uint64_t)op[2] * 16375509950581339924UL) + ((uint64_t)op[3] * 10353042788260867254UL) + ((uint64_t)op[4] * 8421140006141382512UL) + ((uint64_t)op[5] * 12381927365744934816UL);
	tmp_q[5] = ((uint64_t)op[0] * 1768846766534990688UL) + ((uint64_t)op[1] * 6381366926573354045UL) + ((uint64_t)op[2] * 15049723021932632948UL) + ((uint64_t)op[3] * 16375509950581339924UL) + ((uint64_t)op[4] * 10353042788260867254UL) + ((uint64_t)op[5] * 8421140006141382512UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 79050090476L) + ((((int128)tmp_q[1] * 55007612174L) - ((int128)tmp_q[2] * 98060691984L) + ((int128)tmp_q[3] * 92974941520L) + ((int128)tmp_q[4] * 77510671033L) + ((int128)tmp_q[5] * 17707167516L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 17707167516L) + ((int128)tmp_q[1] * 79050090476L) + ((((int128)tmp_q[2] * 55007612174L) - ((int128)tmp_q[3] * 98060691984L) + ((int128)tmp_q[4] * 92974941520L) + ((int128)tmp_q[5] * 77510671033L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 77510671033L) + ((int128)tmp_q[1] * 17707167516L) + ((int128)tmp_q[2] * 79050090476L) + ((((int128)tmp_q[3] * 55007612174L) - ((int128)tmp_q[4] * 98060691984L) + ((int128)tmp_q[5] * 92974941520L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 92974941520L) + ((int128)tmp_q[1] * 77510671033L) + ((int128)tmp_q[2] * 17707167516L) + ((int128)tmp_q[3] * 79050090476L) + ((((int128)tmp_q[4] * 55007612174L) - ((int128)tmp_q[5] * 98060691984L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 98060691984L) + ((int128)tmp_q[1] * 92974941520L) + ((int128)tmp_q[2] * 77510671033L) + ((int128)tmp_q[3] * 17707167516L) + ((int128)tmp_q[4] * 79050090476L) + ((int128)tmp_q[5] * 385053285218L);
	tmp_zero[5] = ((int128)tmp_q[0] * 55007612174L) - ((int128)tmp_q[1] * 98060691984L) + ((int128)tmp_q[2] * 92974941520L) + ((int128)tmp_q[3] * 77510671033L) + ((int128)tmp_q[4] * 17707167516L) + ((int128)tmp_q[5] * 79050090476L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

