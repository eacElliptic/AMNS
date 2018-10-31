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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 306084390924305781UL) + ((((uint64_t)op[1] * 16103599056347512265UL) + ((uint64_t)op[2] * 1254349270587588844UL) + ((uint64_t)op[3] * 8338158024365027328UL) + ((uint64_t)op[4] * 16191510846242955345UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 16191510846242955345UL) + ((uint64_t)op[1] * 306084390924305781UL) + ((((uint64_t)op[2] * 16103599056347512265UL) + ((uint64_t)op[3] * 1254349270587588844UL) + ((uint64_t)op[4] * 8338158024365027328UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 8338158024365027328UL) + ((uint64_t)op[1] * 16191510846242955345UL) + ((uint64_t)op[2] * 306084390924305781UL) + ((((uint64_t)op[3] * 16103599056347512265UL) + ((uint64_t)op[4] * 1254349270587588844UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 1254349270587588844UL) + ((uint64_t)op[1] * 8338158024365027328UL) + ((uint64_t)op[2] * 16191510846242955345UL) + ((uint64_t)op[3] * 306084390924305781UL) + ((uint64_t)op[4] * 11715725086810196755UL);
	tmp_q[4] = ((uint64_t)op[0] * 16103599056347512265UL) + ((uint64_t)op[1] * 1254349270587588844UL) + ((uint64_t)op[2] * 8338158024365027328UL) + ((uint64_t)op[3] * 16191510846242955345UL) + ((uint64_t)op[4] * 306084390924305781UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 21605338186537L) - ((-((int128)tmp_q[1] * 13915209661698L) + ((int128)tmp_q[2] * 10421528068441L) - ((int128)tmp_q[3] * 5956633765987L) - ((int128)tmp_q[4] * 11097336865620L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 11097336865620L) + ((int128)tmp_q[1] * 21605338186537L) - ((-((int128)tmp_q[2] * 13915209661698L) + ((int128)tmp_q[3] * 10421528068441L) - ((int128)tmp_q[4] * 5956633765987L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 5956633765987L) - ((int128)tmp_q[1] * 11097336865620L) + ((int128)tmp_q[2] * 21605338186537L) - ((-((int128)tmp_q[3] * 13915209661698L) + ((int128)tmp_q[4] * 10421528068441L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 10421528068441L) - ((int128)tmp_q[1] * 5956633765987L) - ((int128)tmp_q[2] * 11097336865620L) + ((int128)tmp_q[3] * 21605338186537L) + ((int128)tmp_q[4] * 69576048308490L);
	tmp_zero[4] = -((int128)tmp_q[0] * 13915209661698L) + ((int128)tmp_q[1] * 10421528068441L) - ((int128)tmp_q[2] * 5956633765987L) - ((int128)tmp_q[3] * 11097336865620L) + ((int128)tmp_q[4] * 21605338186537L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

