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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7156180466544030577UL) + ((((uint64_t)op[1] * 16346285864076192303UL) + ((uint64_t)op[2] * 12275610049858110435UL) + ((uint64_t)op[3] * 11043553856313118224UL) + ((uint64_t)op[4] * 9755704761549381436UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 9755704761549381436UL) + ((uint64_t)op[1] * 7156180466544030577UL) + ((((uint64_t)op[2] * 16346285864076192303UL) + ((uint64_t)op[3] * 12275610049858110435UL) + ((uint64_t)op[4] * 11043553856313118224UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 11043553856313118224UL) + ((uint64_t)op[1] * 9755704761549381436UL) + ((uint64_t)op[2] * 7156180466544030577UL) + ((((uint64_t)op[3] * 16346285864076192303UL) + ((uint64_t)op[4] * 12275610049858110435UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 12275610049858110435UL) + ((uint64_t)op[1] * 11043553856313118224UL) + ((uint64_t)op[2] * 9755704761549381436UL) + ((uint64_t)op[3] * 7156180466544030577UL) + ((uint64_t)op[4] * 3743536606276036425UL);
	tmp_q[4] = ((uint64_t)op[0] * 16346285864076192303UL) + ((uint64_t)op[1] * 12275610049858110435UL) + ((uint64_t)op[2] * 11043553856313118224UL) + ((uint64_t)op[3] * 9755704761549381436UL) + ((uint64_t)op[4] * 7156180466544030577UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5714678129428L) + ((((int128)tmp_q[1] * 7155984959175L) - ((int128)tmp_q[2] * 9932928067435L) - ((int128)tmp_q[3] * 15354881156612L) - ((int128)tmp_q[4] * 12820057864005L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 12820057864005L) + ((int128)tmp_q[1] * 5714678129428L) + ((((int128)tmp_q[2] * 7155984959175L) - ((int128)tmp_q[3] * 9932928067435L) - ((int128)tmp_q[4] * 15354881156612L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 15354881156612L) - ((int128)tmp_q[1] * 12820057864005L) + ((int128)tmp_q[2] * 5714678129428L) + ((((int128)tmp_q[3] * 7155984959175L) - ((int128)tmp_q[4] * 9932928067435L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 9932928067435L) - ((int128)tmp_q[1] * 15354881156612L) - ((int128)tmp_q[2] * 12820057864005L) + ((int128)tmp_q[3] * 5714678129428L) + ((int128)tmp_q[4] * 50091894714225L);
	tmp_zero[4] = ((int128)tmp_q[0] * 7155984959175L) - ((int128)tmp_q[1] * 9932928067435L) - ((int128)tmp_q[2] * 15354881156612L) - ((int128)tmp_q[3] * 12820057864005L) + ((int128)tmp_q[4] * 5714678129428L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

