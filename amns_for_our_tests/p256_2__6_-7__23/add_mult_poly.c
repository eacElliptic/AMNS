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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3393341686523354415UL) + ((((uint64_t)op[1] * 3069161338010645369UL) + ((uint64_t)op[2] * 16957060210284376520UL) + ((uint64_t)op[3] * 9212445800735415868UL) + ((uint64_t)op[4] * 3858668959064655797UL) + ((uint64_t)op[5] * 11426562032955140868UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 11426562032955140868UL) + ((uint64_t)op[1] * 3393341686523354415UL) + ((((uint64_t)op[2] * 3069161338010645369UL) + ((uint64_t)op[3] * 16957060210284376520UL) + ((uint64_t)op[4] * 9212445800735415868UL) + ((uint64_t)op[5] * 3858668959064655797UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 3858668959064655797UL) + ((uint64_t)op[1] * 11426562032955140868UL) + ((uint64_t)op[2] * 3393341686523354415UL) + ((((uint64_t)op[3] * 3069161338010645369UL) + ((uint64_t)op[4] * 16957060210284376520UL) + ((uint64_t)op[5] * 9212445800735415868UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 9212445800735415868UL) + ((uint64_t)op[1] * 3858668959064655797UL) + ((uint64_t)op[2] * 11426562032955140868UL) + ((uint64_t)op[3] * 3393341686523354415UL) + ((((uint64_t)op[4] * 3069161338010645369UL) + ((uint64_t)op[5] * 16957060210284376520UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 16957060210284376520UL) + ((uint64_t)op[1] * 9212445800735415868UL) + ((uint64_t)op[2] * 3858668959064655797UL) + ((uint64_t)op[3] * 11426562032955140868UL) + ((uint64_t)op[4] * 3393341686523354415UL) + ((uint64_t)op[5] * 15409358781344585649UL);
	tmp_q[5] = ((uint64_t)op[0] * 3069161338010645369UL) + ((uint64_t)op[1] * 16957060210284376520UL) + ((uint64_t)op[2] * 9212445800735415868UL) + ((uint64_t)op[3] * 3858668959064655797UL) + ((uint64_t)op[4] * 11426562032955140868UL) + ((uint64_t)op[5] * 3393341686523354415UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1367763585443L) - ((((int128)tmp_q[1] * 2985528233393L) + ((int128)tmp_q[2] * 1253131150962L) - ((int128)tmp_q[3] * 1363806765144L) - ((int128)tmp_q[4] * 1255024640291L) - ((int128)tmp_q[5] * 3738431663466L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 3738431663466L) - ((int128)tmp_q[1] * 1367763585443L) - ((((int128)tmp_q[2] * 2985528233393L) + ((int128)tmp_q[3] * 1253131150962L) - ((int128)tmp_q[4] * 1363806765144L) - ((int128)tmp_q[5] * 1255024640291L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 1255024640291L) - ((int128)tmp_q[1] * 3738431663466L) - ((int128)tmp_q[2] * 1367763585443L) - ((((int128)tmp_q[3] * 2985528233393L) + ((int128)tmp_q[4] * 1253131150962L) - ((int128)tmp_q[5] * 1363806765144L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 1363806765144L) - ((int128)tmp_q[1] * 1255024640291L) - ((int128)tmp_q[2] * 3738431663466L) - ((int128)tmp_q[3] * 1367763585443L) - ((((int128)tmp_q[4] * 2985528233393L) + ((int128)tmp_q[5] * 1253131150962L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 1253131150962L) - ((int128)tmp_q[1] * 1363806765144L) - ((int128)tmp_q[2] * 1255024640291L) - ((int128)tmp_q[3] * 3738431663466L) - ((int128)tmp_q[4] * 1367763585443L) - ((int128)tmp_q[5] * 20898697633751L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2985528233393L) + ((int128)tmp_q[1] * 1253131150962L) - ((int128)tmp_q[2] * 1363806765144L) - ((int128)tmp_q[3] * 1255024640291L) - ((int128)tmp_q[4] * 3738431663466L) - ((int128)tmp_q[5] * 1367763585443L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

