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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3273232495486360784UL) + ((((uint64_t)op[1] * 12991929538832666625UL) + ((uint64_t)op[2] * 8113527758179942166UL) + ((uint64_t)op[3] * 10642455780496794997UL) + ((uint64_t)op[4] * 8087204387153167883UL) + ((uint64_t)op[5] * 8063533835868457528UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 8063533835868457528UL) + ((uint64_t)op[1] * 3273232495486360784UL) + ((((uint64_t)op[2] * 12991929538832666625UL) + ((uint64_t)op[3] * 8113527758179942166UL) + ((uint64_t)op[4] * 10642455780496794997UL) + ((uint64_t)op[5] * 8087204387153167883UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 8087204387153167883UL) + ((uint64_t)op[1] * 8063533835868457528UL) + ((uint64_t)op[2] * 3273232495486360784UL) + ((((uint64_t)op[3] * 12991929538832666625UL) + ((uint64_t)op[4] * 8113527758179942166UL) + ((uint64_t)op[5] * 10642455780496794997UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 10642455780496794997UL) + ((uint64_t)op[1] * 8087204387153167883UL) + ((uint64_t)op[2] * 8063533835868457528UL) + ((uint64_t)op[3] * 3273232495486360784UL) + ((((uint64_t)op[4] * 12991929538832666625UL) + ((uint64_t)op[5] * 8113527758179942166UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 8113527758179942166UL) + ((uint64_t)op[1] * 10642455780496794997UL) + ((uint64_t)op[2] * 8087204387153167883UL) + ((uint64_t)op[3] * 8063533835868457528UL) + ((uint64_t)op[4] * 3273232495486360784UL) + ((uint64_t)op[5] * 16364443604630654973UL);
	tmp_q[5] = ((uint64_t)op[0] * 12991929538832666625UL) + ((uint64_t)op[1] * 8113527758179942166UL) + ((uint64_t)op[2] * 10642455780496794997UL) + ((uint64_t)op[3] * 8087204387153167883UL) + ((uint64_t)op[4] * 8063533835868457528UL) + ((uint64_t)op[5] * 3273232495486360784UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 525594884L) - ((((int128)tmp_q[1] * 344273101L) - ((int128)tmp_q[2] * 202544964L) - ((int128)tmp_q[3] * 167651901L) - ((int128)tmp_q[4] * 1486809597L) + ((int128)tmp_q[5] * 3737863582L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 3737863582L) - ((int128)tmp_q[1] * 525594884L) - ((((int128)tmp_q[2] * 344273101L) - ((int128)tmp_q[3] * 202544964L) - ((int128)tmp_q[4] * 167651901L) - ((int128)tmp_q[5] * 1486809597L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1486809597L) + ((int128)tmp_q[1] * 3737863582L) - ((int128)tmp_q[2] * 525594884L) - ((((int128)tmp_q[3] * 344273101L) - ((int128)tmp_q[4] * 202544964L) - ((int128)tmp_q[5] * 167651901L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 167651901L) - ((int128)tmp_q[1] * 1486809597L) + ((int128)tmp_q[2] * 3737863582L) - ((int128)tmp_q[3] * 525594884L) - ((((int128)tmp_q[4] * 344273101L) - ((int128)tmp_q[5] * 202544964L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 202544964L) - ((int128)tmp_q[1] * 167651901L) - ((int128)tmp_q[2] * 1486809597L) + ((int128)tmp_q[3] * 3737863582L) - ((int128)tmp_q[4] * 525594884L) - ((int128)tmp_q[5] * 1032819303L);
	tmp_zero[5] = ((int128)tmp_q[0] * 344273101L) - ((int128)tmp_q[1] * 202544964L) - ((int128)tmp_q[2] * 167651901L) - ((int128)tmp_q[3] * 1486809597L) + ((int128)tmp_q[4] * 3737863582L) - ((int128)tmp_q[5] * 525594884L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

