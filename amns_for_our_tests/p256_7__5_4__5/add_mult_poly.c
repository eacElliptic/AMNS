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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10211738422238389713UL) + ((((uint64_t)op[1] * 540214429955852388UL) + ((uint64_t)op[2] * 4952234588169196913UL) + ((uint64_t)op[3] * 8602191454897955424UL) + ((uint64_t)op[4] * 3250352560508507174UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 3250352560508507174UL) + ((uint64_t)op[1] * 10211738422238389713UL) + ((((uint64_t)op[2] * 540214429955852388UL) + ((uint64_t)op[3] * 4952234588169196913UL) + ((uint64_t)op[4] * 8602191454897955424UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 8602191454897955424UL) + ((uint64_t)op[1] * 3250352560508507174UL) + ((uint64_t)op[2] * 10211738422238389713UL) + ((((uint64_t)op[3] * 540214429955852388UL) + ((uint64_t)op[4] * 4952234588169196913UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 4952234588169196913UL) + ((uint64_t)op[1] * 8602191454897955424UL) + ((uint64_t)op[2] * 3250352560508507174UL) + ((uint64_t)op[3] * 10211738422238389713UL) + ((uint64_t)op[4] * 2160857719823409552UL);
	tmp_q[4] = ((uint64_t)op[0] * 540214429955852388UL) + ((uint64_t)op[1] * 4952234588169196913UL) + ((uint64_t)op[2] * 8602191454897955424UL) + ((uint64_t)op[3] * 3250352560508507174UL) + ((uint64_t)op[4] * 10211738422238389713UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 249190801853887L) + ((((int128)tmp_q[1] * 141427559599180L) - ((int128)tmp_q[2] * 778492527340695L) - ((int128)tmp_q[3] * 875342413436076L) - ((int128)tmp_q[4] * 979268721475230L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 979268721475230L) + ((int128)tmp_q[1] * 249190801853887L) + ((((int128)tmp_q[2] * 141427559599180L) - ((int128)tmp_q[3] * 778492527340695L) - ((int128)tmp_q[4] * 875342413436076L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 875342413436076L) - ((int128)tmp_q[1] * 979268721475230L) + ((int128)tmp_q[2] * 249190801853887L) + ((((int128)tmp_q[3] * 141427559599180L) - ((int128)tmp_q[4] * 778492527340695L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 778492527340695L) - ((int128)tmp_q[1] * 875342413436076L) - ((int128)tmp_q[2] * 979268721475230L) + ((int128)tmp_q[3] * 249190801853887L) + ((int128)tmp_q[4] * 565710238396720L);
	tmp_zero[4] = ((int128)tmp_q[0] * 141427559599180L) - ((int128)tmp_q[1] * 778492527340695L) - ((int128)tmp_q[2] * 875342413436076L) - ((int128)tmp_q[3] * 979268721475230L) + ((int128)tmp_q[4] * 249190801853887L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

