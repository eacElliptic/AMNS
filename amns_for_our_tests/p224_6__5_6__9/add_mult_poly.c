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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8366232400863008899UL) + ((((uint64_t)op[1] * 3205606965409686778UL) + ((uint64_t)op[2] * 14533902511416087781UL) + ((uint64_t)op[3] * 17665193166791952864UL) + ((uint64_t)op[4] * 16009947092498503812UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 16009947092498503812UL) + ((uint64_t)op[1] * 8366232400863008899UL) + ((((uint64_t)op[2] * 3205606965409686778UL) + ((uint64_t)op[3] * 14533902511416087781UL) + ((uint64_t)op[4] * 17665193166791952864UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 17665193166791952864UL) + ((uint64_t)op[1] * 16009947092498503812UL) + ((uint64_t)op[2] * 8366232400863008899UL) + ((((uint64_t)op[3] * 3205606965409686778UL) + ((uint64_t)op[4] * 14533902511416087781UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 14533902511416087781UL) + ((uint64_t)op[1] * 17665193166791952864UL) + ((uint64_t)op[2] * 16009947092498503812UL) + ((uint64_t)op[3] * 8366232400863008899UL) + ((uint64_t)op[4] * 786897718748569052UL);
	tmp_q[4] = ((uint64_t)op[0] * 3205606965409686778UL) + ((uint64_t)op[1] * 14533902511416087781UL) + ((uint64_t)op[2] * 17665193166791952864UL) + ((uint64_t)op[3] * 16009947092498503812UL) + ((uint64_t)op[4] * 8366232400863008899UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2433543799355L) + ((((int128)tmp_q[1] * 7269471199056L) - ((int128)tmp_q[2] * 7872980688579L) + ((int128)tmp_q[3] * 8371119371780L) - ((int128)tmp_q[4] * 8404128093134L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 8404128093134L) - ((int128)tmp_q[1] * 2433543799355L) + ((((int128)tmp_q[2] * 7269471199056L) - ((int128)tmp_q[3] * 7872980688579L) + ((int128)tmp_q[4] * 8371119371780L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 8371119371780L) - ((int128)tmp_q[1] * 8404128093134L) - ((int128)tmp_q[2] * 2433543799355L) + ((((int128)tmp_q[3] * 7269471199056L) - ((int128)tmp_q[4] * 7872980688579L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 7872980688579L) + ((int128)tmp_q[1] * 8371119371780L) - ((int128)tmp_q[2] * 8404128093134L) - ((int128)tmp_q[3] * 2433543799355L) + ((int128)tmp_q[4] * 43616827194336L);
	tmp_zero[4] = ((int128)tmp_q[0] * 7269471199056L) - ((int128)tmp_q[1] * 7872980688579L) + ((int128)tmp_q[2] * 8371119371780L) - ((int128)tmp_q[3] * 8404128093134L) - ((int128)tmp_q[4] * 2433543799355L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

