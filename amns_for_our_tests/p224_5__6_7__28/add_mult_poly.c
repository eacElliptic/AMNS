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
	tmp_q[0] = ((uint64_t)op[0] * 15574730655776365842UL) + ((((uint64_t)op[1] * 16241749441165806690UL) + ((uint64_t)op[2] * 18061815049820671842UL) + ((uint64_t)op[3] * 66309605393520128UL) + ((uint64_t)op[4] * 584117201948519610UL) + ((uint64_t)op[5] * 16703335326215492735UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 16703335326215492735UL) + ((uint64_t)op[1] * 15574730655776365842UL) + ((((uint64_t)op[2] * 16241749441165806690UL) + ((uint64_t)op[3] * 18061815049820671842UL) + ((uint64_t)op[4] * 66309605393520128UL) + ((uint64_t)op[5] * 584117201948519610UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 584117201948519610UL) + ((uint64_t)op[1] * 16703335326215492735UL) + ((uint64_t)op[2] * 15574730655776365842UL) + ((((uint64_t)op[3] * 16241749441165806690UL) + ((uint64_t)op[4] * 18061815049820671842UL) + ((uint64_t)op[5] * 66309605393520128UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 66309605393520128UL) + ((uint64_t)op[1] * 584117201948519610UL) + ((uint64_t)op[2] * 16703335326215492735UL) + ((uint64_t)op[3] * 15574730655776365842UL) + ((((uint64_t)op[4] * 16241749441165806690UL) + ((uint64_t)op[5] * 18061815049820671842UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 18061815049820671842UL) + ((uint64_t)op[1] * 66309605393520128UL) + ((uint64_t)op[2] * 584117201948519610UL) + ((uint64_t)op[3] * 16703335326215492735UL) + ((uint64_t)op[4] * 15574730655776365842UL) + ((uint64_t)op[5] * 3011781645903337134UL);
	tmp_q[5] = ((uint64_t)op[0] * 16241749441165806690UL) + ((uint64_t)op[1] * 18061815049820671842UL) + ((uint64_t)op[2] * 66309605393520128UL) + ((uint64_t)op[3] * 584117201948519610UL) + ((uint64_t)op[4] * 16703335326215492735UL) + ((uint64_t)op[5] * 15574730655776365842UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 66835511342L) + ((((int128)tmp_q[1] * 53588843875L) - ((int128)tmp_q[2] * 141108321498L) + ((int128)tmp_q[3] * 87905381758L) - ((int128)tmp_q[4] * 39263408326L) + ((int128)tmp_q[5] * 23486795296L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 23486795296L) - ((int128)tmp_q[1] * 66835511342L) + ((((int128)tmp_q[2] * 53588843875L) - ((int128)tmp_q[3] * 141108321498L) + ((int128)tmp_q[4] * 87905381758L) - ((int128)tmp_q[5] * 39263408326L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 39263408326L) + ((int128)tmp_q[1] * 23486795296L) - ((int128)tmp_q[2] * 66835511342L) + ((((int128)tmp_q[3] * 53588843875L) - ((int128)tmp_q[4] * 141108321498L) + ((int128)tmp_q[5] * 87905381758L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 87905381758L) - ((int128)tmp_q[1] * 39263408326L) + ((int128)tmp_q[2] * 23486795296L) - ((int128)tmp_q[3] * 66835511342L) + ((((int128)tmp_q[4] * 53588843875L) - ((int128)tmp_q[5] * 141108321498L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 141108321498L) + ((int128)tmp_q[1] * 87905381758L) - ((int128)tmp_q[2] * 39263408326L) + ((int128)tmp_q[3] * 23486795296L) - ((int128)tmp_q[4] * 66835511342L) + ((int128)tmp_q[5] * 375121907125L);
	tmp_zero[5] = ((int128)tmp_q[0] * 53588843875L) - ((int128)tmp_q[1] * 141108321498L) + ((int128)tmp_q[2] * 87905381758L) - ((int128)tmp_q[3] * 39263408326L) + ((int128)tmp_q[4] * 23486795296L) - ((int128)tmp_q[5] * 66835511342L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

