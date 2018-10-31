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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6038036756291059145UL) + ((((uint64_t)op[1] * 1115878358952112494UL) + ((uint64_t)op[2] * 6913763821158272612UL) + ((uint64_t)op[3] * 6904687245259253910UL) + ((uint64_t)op[4] * 6521205463818712610UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 6521205463818712610UL) + ((uint64_t)op[1] * 6038036756291059145UL) + ((((uint64_t)op[2] * 1115878358952112494UL) + ((uint64_t)op[3] * 6913763821158272612UL) + ((uint64_t)op[4] * 6904687245259253910UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 6904687245259253910UL) + ((uint64_t)op[1] * 6521205463818712610UL) + ((uint64_t)op[2] * 6038036756291059145UL) + ((((uint64_t)op[3] * 1115878358952112494UL) + ((uint64_t)op[4] * 6913763821158272612UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 6913763821158272612UL) + ((uint64_t)op[1] * 6904687245259253910UL) + ((uint64_t)op[2] * 6521205463818712610UL) + ((uint64_t)op[3] * 6038036756291059145UL) + ((uint64_t)op[4] * 5579391794760562470UL);
	tmp_q[4] = ((uint64_t)op[0] * 1115878358952112494UL) + ((uint64_t)op[1] * 6913763821158272612UL) + ((uint64_t)op[2] * 6904687245259253910UL) + ((uint64_t)op[3] * 6521205463818712610UL) + ((uint64_t)op[4] * 6038036756291059145UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 110624857871L) + ((-((int128)tmp_q[1] * 75001024854L) - ((int128)tmp_q[2] * 156224839976L) - ((int128)tmp_q[3] * 64403408070L) - ((int128)tmp_q[4] * 179754153302L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 179754153302L) + ((int128)tmp_q[1] * 110624857871L) + ((-((int128)tmp_q[2] * 75001024854L) - ((int128)tmp_q[3] * 156224839976L) - ((int128)tmp_q[4] * 64403408070L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 64403408070L) - ((int128)tmp_q[1] * 179754153302L) + ((int128)tmp_q[2] * 110624857871L) + ((-((int128)tmp_q[3] * 75001024854L) - ((int128)tmp_q[4] * 156224839976L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 156224839976L) - ((int128)tmp_q[1] * 64403408070L) - ((int128)tmp_q[2] * 179754153302L) + ((int128)tmp_q[3] * 110624857871L) - ((int128)tmp_q[4] * 375005124270L);
	tmp_zero[4] = -((int128)tmp_q[0] * 75001024854L) - ((int128)tmp_q[1] * 156224839976L) - ((int128)tmp_q[2] * 64403408070L) - ((int128)tmp_q[3] * 179754153302L) + ((int128)tmp_q[4] * 110624857871L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

