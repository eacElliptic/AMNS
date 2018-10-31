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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16392443008074936907UL) + ((((uint64_t)op[1] * 3381611317052853511UL) + ((uint64_t)op[2] * 6678322058755481378UL) + ((uint64_t)op[3] * 16339283945330421219UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 16339283945330421219UL) + ((uint64_t)op[1] * 16392443008074936907UL) + ((((uint64_t)op[2] * 3381611317052853511UL) + ((uint64_t)op[3] * 6678322058755481378UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 6678322058755481378UL) + ((uint64_t)op[1] * 16339283945330421219UL) + ((uint64_t)op[2] * 16392443008074936907UL) + ((uint64_t)op[3] * 9840597610996275144UL);
	tmp_q[3] = ((uint64_t)op[0] * 3381611317052853511UL) + ((uint64_t)op[1] * 6678322058755481378UL) + ((uint64_t)op[2] * 16339283945330421219UL) + ((uint64_t)op[3] * 16392443008074936907UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 32589747138581L) - ((-((int128)tmp_q[1] * 74487986815954L) + ((int128)tmp_q[2] * 113630220479103L) - ((int128)tmp_q[3] * 124448669919877L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 124448669919877L) + ((int128)tmp_q[1] * 32589747138581L) - ((-((int128)tmp_q[2] * 74487986815954L) + ((int128)tmp_q[3] * 113630220479103L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 113630220479103L) - ((int128)tmp_q[1] * 124448669919877L) + ((int128)tmp_q[2] * 32589747138581L) + ((int128)tmp_q[3] * 595903894527632L);
	tmp_zero[3] = -((int128)tmp_q[0] * 74487986815954L) + ((int128)tmp_q[1] * 113630220479103L) - ((int128)tmp_q[2] * 124448669919877L) + ((int128)tmp_q[3] * 32589747138581L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

