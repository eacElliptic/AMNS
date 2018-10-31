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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15643430833819115965UL) + ((((uint64_t)op[1] * 13080102749201511787UL) + ((uint64_t)op[2] * 17598401686789740085UL) + ((uint64_t)op[3] * 4240440009166231706UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 4240440009166231706UL) + ((uint64_t)op[1] * 15643430833819115965UL) + ((((uint64_t)op[2] * 13080102749201511787UL) + ((uint64_t)op[3] * 17598401686789740085UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 17598401686789740085UL) + ((uint64_t)op[1] * 4240440009166231706UL) + ((uint64_t)op[2] * 15643430833819115965UL) + ((uint64_t)op[3] * 8386462548830647529UL);
	tmp_q[3] = ((uint64_t)op[0] * 13080102749201511787UL) + ((uint64_t)op[1] * 17598401686789740085UL) + ((uint64_t)op[2] * 4240440009166231706UL) + ((uint64_t)op[3] * 15643430833819115965UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 84929463253223L) - ((-((int128)tmp_q[1] * 77730427106030L) + ((int128)tmp_q[2] * 124824828094243L) - ((int128)tmp_q[3] * 15282062749183L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 15282062749183L) - ((int128)tmp_q[1] * 84929463253223L) - ((-((int128)tmp_q[2] * 77730427106030L) + ((int128)tmp_q[3] * 124824828094243L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 124824828094243L) - ((int128)tmp_q[1] * 15282062749183L) - ((int128)tmp_q[2] * 84929463253223L) + ((int128)tmp_q[3] * 388652135530150L);
	tmp_zero[3] = -((int128)tmp_q[0] * 77730427106030L) + ((int128)tmp_q[1] * 124824828094243L) - ((int128)tmp_q[2] * 15282062749183L) - ((int128)tmp_q[3] * 84929463253223L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

