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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[3] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 362000387875473835UL) + ((((uint64_t)op[1] * 11429064691452592034UL) + ((uint64_t)op[2] * 17657565747270687557UL) + ((uint64_t)op[3] * 5391404411155036114UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 5391404411155036114UL) + ((uint64_t)op[1] * 362000387875473835UL) + ((((uint64_t)op[2] * 11429064691452592034UL) + ((uint64_t)op[3] * 17657565747270687557UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 17657565747270687557UL) + ((uint64_t)op[1] * 5391404411155036114UL) + ((uint64_t)op[2] * 362000387875473835UL) + ((uint64_t)op[3] * 17645541236782529808UL);
	tmp_q[3] = ((uint64_t)op[0] * 11429064691452592034UL) + ((uint64_t)op[1] * 17657565747270687557UL) + ((uint64_t)op[2] * 5391404411155036114UL) + ((uint64_t)op[3] * 362000387875473835UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 287092079088581L) + ((((int128)tmp_q[1] * 88845410694702L) - ((int128)tmp_q[2] * 123588477907063L) + ((int128)tmp_q[3] * 117544820143090L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 117544820143090L) + ((int128)tmp_q[1] * 287092079088581L) + ((((int128)tmp_q[2] * 88845410694702L) - ((int128)tmp_q[3] * 123588477907063L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 123588477907063L) + ((int128)tmp_q[1] * 117544820143090L) + ((int128)tmp_q[2] * 287092079088581L) + ((int128)tmp_q[3] * 710763285557616L);
	tmp_zero[3] = ((int128)tmp_q[0] * 88845410694702L) - ((int128)tmp_q[1] * 123588477907063L) + ((int128)tmp_q[2] * 117544820143090L) + ((int128)tmp_q[3] * 287092079088581L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

