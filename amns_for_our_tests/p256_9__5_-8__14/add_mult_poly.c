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
	tmp_q[0] = ((uint64_t)op[0] * 6010684834938645567UL) + ((((uint64_t)op[1] * 3234024420638131673UL) + ((uint64_t)op[2] * 15735943292792400711UL) + ((uint64_t)op[3] * 16891988692935411520UL) + ((uint64_t)op[4] * 7127449372918822186UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 7127449372918822186UL) + ((uint64_t)op[1] * 6010684834938645567UL) + ((((uint64_t)op[2] * 3234024420638131673UL) + ((uint64_t)op[3] * 15735943292792400711UL) + ((uint64_t)op[4] * 16891988692935411520UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 16891988692935411520UL) + ((uint64_t)op[1] * 7127449372918822186UL) + ((uint64_t)op[2] * 6010684834938645567UL) + ((((uint64_t)op[3] * 3234024420638131673UL) + ((uint64_t)op[4] * 15735943292792400711UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 15735943292792400711UL) + ((uint64_t)op[1] * 16891988692935411520UL) + ((uint64_t)op[2] * 7127449372918822186UL) + ((uint64_t)op[3] * 6010684834938645567UL) + ((uint64_t)op[4] * 11021292782314049848UL);
	tmp_q[4] = ((uint64_t)op[0] * 3234024420638131673UL) + ((uint64_t)op[1] * 15735943292792400711UL) + ((uint64_t)op[2] * 16891988692935411520UL) + ((uint64_t)op[3] * 7127449372918822186UL) + ((uint64_t)op[4] * 6010684834938645567UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 875000193895231L) - ((((int128)tmp_q[1] * 694004850144397L) - ((int128)tmp_q[2] * 1383707940694585L) - ((int128)tmp_q[3] * 824888735225532L) - ((int128)tmp_q[4] * 374379730078590L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 374379730078590L) - ((int128)tmp_q[1] * 875000193895231L) - ((((int128)tmp_q[2] * 694004850144397L) - ((int128)tmp_q[3] * 1383707940694585L) - ((int128)tmp_q[4] * 824888735225532L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 824888735225532L) - ((int128)tmp_q[1] * 374379730078590L) - ((int128)tmp_q[2] * 875000193895231L) - ((((int128)tmp_q[3] * 694004850144397L) - ((int128)tmp_q[4] * 1383707940694585L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 1383707940694585L) - ((int128)tmp_q[1] * 824888735225532L) - ((int128)tmp_q[2] * 374379730078590L) - ((int128)tmp_q[3] * 875000193895231L) - ((int128)tmp_q[4] * 5552038801155176L);
	tmp_zero[4] = ((int128)tmp_q[0] * 694004850144397L) - ((int128)tmp_q[1] * 1383707940694585L) - ((int128)tmp_q[2] * 824888735225532L) - ((int128)tmp_q[3] * 374379730078590L) - ((int128)tmp_q[4] * 875000193895231L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

