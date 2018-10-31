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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12258654117409965393UL) + ((((uint64_t)op[1] * 806852976784714608UL) + ((uint64_t)op[2] * 12409398074779057359UL) + ((uint64_t)op[3] * 12777014871664759163UL) + ((uint64_t)op[4] * 8827993689755004339UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 8827993689755004339UL) + ((uint64_t)op[1] * 12258654117409965393UL) + ((((uint64_t)op[2] * 806852976784714608UL) + ((uint64_t)op[3] * 12409398074779057359UL) + ((uint64_t)op[4] * 12777014871664759163UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 12777014871664759163UL) + ((uint64_t)op[1] * 8827993689755004339UL) + ((uint64_t)op[2] * 12258654117409965393UL) + ((((uint64_t)op[3] * 806852976784714608UL) + ((uint64_t)op[4] * 12409398074779057359UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 12409398074779057359UL) + ((uint64_t)op[1] * 12777014871664759163UL) + ((uint64_t)op[2] * 8827993689755004339UL) + ((uint64_t)op[3] * 12258654117409965393UL) + ((uint64_t)op[4] * 11991920259431834752UL);
	tmp_q[4] = ((uint64_t)op[0] * 806852976784714608UL) + ((uint64_t)op[1] * 12409398074779057359UL) + ((uint64_t)op[2] * 12777014871664759163UL) + ((uint64_t)op[3] * 8827993689755004339UL) + ((uint64_t)op[4] * 12258654117409965393UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 403199675913L) - ((((int128)tmp_q[1] * 115898106573L) + ((int128)tmp_q[2] * 142158424384L) + ((int128)tmp_q[3] * 62085230274L) - ((int128)tmp_q[4] * 38734367405L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 38734367405L) - ((int128)tmp_q[1] * 403199675913L) - ((((int128)tmp_q[2] * 115898106573L) + ((int128)tmp_q[3] * 142158424384L) + ((int128)tmp_q[4] * 62085230274L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 62085230274L) - ((int128)tmp_q[1] * 38734367405L) - ((int128)tmp_q[2] * 403199675913L) - ((((int128)tmp_q[3] * 115898106573L) + ((int128)tmp_q[4] * 142158424384L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 142158424384L) + ((int128)tmp_q[1] * 62085230274L) - ((int128)tmp_q[2] * 38734367405L) - ((int128)tmp_q[3] * 403199675913L) - ((int128)tmp_q[4] * 927184852584L);
	tmp_zero[4] = ((int128)tmp_q[0] * 115898106573L) + ((int128)tmp_q[1] * 142158424384L) + ((int128)tmp_q[2] * 62085230274L) - ((int128)tmp_q[3] * 38734367405L) - ((int128)tmp_q[4] * 403199675913L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

