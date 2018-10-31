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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[3] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14404731248266644585UL) + ((((uint64_t)op[1] * 3286953480751773897UL) + ((uint64_t)op[2] * 5408557101331171271UL) + ((uint64_t)op[3] * 6945810801382881972UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 6945810801382881972UL) + ((uint64_t)op[1] * 14404731248266644585UL) + ((((uint64_t)op[2] * 3286953480751773897UL) + ((uint64_t)op[3] * 5408557101331171271UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 5408557101331171271UL) + ((uint64_t)op[1] * 6945810801382881972UL) + ((uint64_t)op[2] * 14404731248266644585UL) + ((uint64_t)op[3] * 16434767403758869485UL);
	tmp_q[3] = ((uint64_t)op[0] * 3286953480751773897UL) + ((uint64_t)op[1] * 5408557101331171271UL) + ((uint64_t)op[2] * 6945810801382881972UL) + ((uint64_t)op[3] * 14404731248266644585UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 128659502648371L) + ((((int128)tmp_q[1] * 95477151974086L) + ((int128)tmp_q[2] * 95597758307649L) + ((int128)tmp_q[3] * 99992787138717L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 99992787138717L) + ((int128)tmp_q[1] * 128659502648371L) + ((((int128)tmp_q[2] * 95477151974086L) + ((int128)tmp_q[3] * 95597758307649L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 95597758307649L) + ((int128)tmp_q[1] * 99992787138717L) + ((int128)tmp_q[2] * 128659502648371L) + ((int128)tmp_q[3] * 477385759870430L);
	tmp_zero[3] = ((int128)tmp_q[0] * 95477151974086L) + ((int128)tmp_q[1] * 95597758307649L) + ((int128)tmp_q[2] * 99992787138717L) + ((int128)tmp_q[3] * 128659502648371L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

