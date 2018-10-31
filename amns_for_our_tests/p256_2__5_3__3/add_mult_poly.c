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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17166486105417144355UL) + ((((uint64_t)op[1] * 13186406228919152681UL) + ((uint64_t)op[2] * 14674995656935540949UL) + ((uint64_t)op[3] * 10661276947901259384UL) + ((uint64_t)op[4] * 12960043309879026576UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 12960043309879026576UL) + ((uint64_t)op[1] * 17166486105417144355UL) + ((((uint64_t)op[2] * 13186406228919152681UL) + ((uint64_t)op[3] * 14674995656935540949UL) + ((uint64_t)op[4] * 10661276947901259384UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 10661276947901259384UL) + ((uint64_t)op[1] * 12960043309879026576UL) + ((uint64_t)op[2] * 17166486105417144355UL) + ((((uint64_t)op[3] * 13186406228919152681UL) + ((uint64_t)op[4] * 14674995656935540949UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 14674995656935540949UL) + ((uint64_t)op[1] * 10661276947901259384UL) + ((uint64_t)op[2] * 12960043309879026576UL) + ((uint64_t)op[3] * 17166486105417144355UL) + ((uint64_t)op[4] * 2665730539338354811UL);
	tmp_q[4] = ((uint64_t)op[0] * 13186406228919152681UL) + ((uint64_t)op[1] * 14674995656935540949UL) + ((uint64_t)op[2] * 10661276947901259384UL) + ((uint64_t)op[3] * 12960043309879026576UL) + ((uint64_t)op[4] * 17166486105417144355UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1700159244628576L) + ((-((int128)tmp_q[1] * 105766970035815L) + ((int128)tmp_q[2] * 1290709848777687L) - ((int128)tmp_q[3] * 1476457361151084L) + ((int128)tmp_q[4] * 870581539434577L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 870581539434577L) + ((int128)tmp_q[1] * 1700159244628576L) + ((-((int128)tmp_q[2] * 105766970035815L) + ((int128)tmp_q[3] * 1290709848777687L) - ((int128)tmp_q[4] * 1476457361151084L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1476457361151084L) + ((int128)tmp_q[1] * 870581539434577L) + ((int128)tmp_q[2] * 1700159244628576L) + ((-((int128)tmp_q[3] * 105766970035815L) + ((int128)tmp_q[4] * 1290709848777687L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1290709848777687L) - ((int128)tmp_q[1] * 1476457361151084L) + ((int128)tmp_q[2] * 870581539434577L) + ((int128)tmp_q[3] * 1700159244628576L) - ((int128)tmp_q[4] * 317300910107445L);
	tmp_zero[4] = -((int128)tmp_q[0] * 105766970035815L) + ((int128)tmp_q[1] * 1290709848777687L) - ((int128)tmp_q[2] * 1476457361151084L) + ((int128)tmp_q[3] * 870581539434577L) + ((int128)tmp_q[4] * 1700159244628576L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

