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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9308000454732092363UL) + ((((uint64_t)op[1] * 9324754700691716704UL) + ((uint64_t)op[2] * 13526086868177862572UL) + ((uint64_t)op[3] * 17539165994719966075UL) + ((uint64_t)op[4] * 17586993231496953867UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 17586993231496953867UL) + ((uint64_t)op[1] * 9308000454732092363UL) + ((((uint64_t)op[2] * 9324754700691716704UL) + ((uint64_t)op[3] * 13526086868177862572UL) + ((uint64_t)op[4] * 17539165994719966075UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 17539165994719966075UL) + ((uint64_t)op[1] * 17586993231496953867UL) + ((uint64_t)op[2] * 9308000454732092363UL) + ((((uint64_t)op[3] * 9324754700691716704UL) + ((uint64_t)op[4] * 13526086868177862572UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 13526086868177862572UL) + ((uint64_t)op[1] * 17539165994719966075UL) + ((uint64_t)op[2] * 17586993231496953867UL) + ((uint64_t)op[3] * 9308000454732092363UL) + ((uint64_t)op[4] * 202765327673881792UL);
	tmp_q[4] = ((uint64_t)op[0] * 9324754700691716704UL) + ((uint64_t)op[1] * 13526086868177862572UL) + ((uint64_t)op[2] * 17539165994719966075UL) + ((uint64_t)op[3] * 17586993231496953867UL) + ((uint64_t)op[4] * 9308000454732092363UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 14513402996919L) + ((((int128)tmp_q[1] * 4921742133381L) - ((int128)tmp_q[2] * 13435201241011L) + ((int128)tmp_q[3] * 11786839148698L) - ((int128)tmp_q[4] * 15931427721627L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 15931427721627L) - ((int128)tmp_q[1] * 14513402996919L) + ((((int128)tmp_q[2] * 4921742133381L) - ((int128)tmp_q[3] * 13435201241011L) + ((int128)tmp_q[4] * 11786839148698L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 11786839148698L) - ((int128)tmp_q[1] * 15931427721627L) - ((int128)tmp_q[2] * 14513402996919L) + ((((int128)tmp_q[3] * 4921742133381L) - ((int128)tmp_q[4] * 13435201241011L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 13435201241011L) + ((int128)tmp_q[1] * 11786839148698L) - ((int128)tmp_q[2] * 15931427721627L) - ((int128)tmp_q[3] * 14513402996919L) + ((int128)tmp_q[4] * 9843484266762L);
	tmp_zero[4] = ((int128)tmp_q[0] * 4921742133381L) - ((int128)tmp_q[1] * 13435201241011L) + ((int128)tmp_q[2] * 11786839148698L) - ((int128)tmp_q[3] * 15931427721627L) - ((int128)tmp_q[4] * 14513402996919L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

