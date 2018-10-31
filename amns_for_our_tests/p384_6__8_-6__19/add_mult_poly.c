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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17658620449639048605UL) + ((((uint64_t)op[1] * 1154505795420835177UL) + ((uint64_t)op[2] * 14599244127102105038UL) + ((uint64_t)op[3] * 18432363850243359809UL) + ((uint64_t)op[4] * 16354839702589602450UL) + ((uint64_t)op[5] * 1842640072758191403UL) + ((uint64_t)op[6] * 9819296124003355144UL) + ((uint64_t)op[7] * 16359183155090143332UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16359183155090143332UL) + ((uint64_t)op[1] * 17658620449639048605UL) + ((((uint64_t)op[2] * 1154505795420835177UL) + ((uint64_t)op[3] * 14599244127102105038UL) + ((uint64_t)op[4] * 18432363850243359809UL) + ((uint64_t)op[5] * 16354839702589602450UL) + ((uint64_t)op[6] * 1842640072758191403UL) + ((uint64_t)op[7] * 9819296124003355144UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 9819296124003355144UL) + ((uint64_t)op[1] * 16359183155090143332UL) + ((uint64_t)op[2] * 17658620449639048605UL) + ((((uint64_t)op[3] * 1154505795420835177UL) + ((uint64_t)op[4] * 14599244127102105038UL) + ((uint64_t)op[5] * 18432363850243359809UL) + ((uint64_t)op[6] * 16354839702589602450UL) + ((uint64_t)op[7] * 1842640072758191403UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 1842640072758191403UL) + ((uint64_t)op[1] * 9819296124003355144UL) + ((uint64_t)op[2] * 16359183155090143332UL) + ((uint64_t)op[3] * 17658620449639048605UL) + ((((uint64_t)op[4] * 1154505795420835177UL) + ((uint64_t)op[5] * 14599244127102105038UL) + ((uint64_t)op[6] * 18432363850243359809UL) + ((uint64_t)op[7] * 16354839702589602450UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 16354839702589602450UL) + ((uint64_t)op[1] * 1842640072758191403UL) + ((uint64_t)op[2] * 9819296124003355144UL) + ((uint64_t)op[3] * 16359183155090143332UL) + ((uint64_t)op[4] * 17658620449639048605UL) + ((((uint64_t)op[5] * 1154505795420835177UL) + ((uint64_t)op[6] * 14599244127102105038UL) + ((uint64_t)op[7] * 18432363850243359809UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 18432363850243359809UL) + ((uint64_t)op[1] * 16354839702589602450UL) + ((uint64_t)op[2] * 1842640072758191403UL) + ((uint64_t)op[3] * 9819296124003355144UL) + ((uint64_t)op[4] * 16359183155090143332UL) + ((uint64_t)op[5] * 17658620449639048605UL) + ((((uint64_t)op[6] * 1154505795420835177UL) + ((uint64_t)op[7] * 14599244127102105038UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 14599244127102105038UL) + ((uint64_t)op[1] * 18432363850243359809UL) + ((uint64_t)op[2] * 16354839702589602450UL) + ((uint64_t)op[3] * 1842640072758191403UL) + ((uint64_t)op[4] * 9819296124003355144UL) + ((uint64_t)op[5] * 16359183155090143332UL) + ((uint64_t)op[6] * 17658620449639048605UL) + ((uint64_t)op[7] * 11519709301184540554UL);
	tmp_q[7] = ((uint64_t)op[0] * 1154505795420835177UL) + ((uint64_t)op[1] * 14599244127102105038UL) + ((uint64_t)op[2] * 18432363850243359809UL) + ((uint64_t)op[3] * 16354839702589602450UL) + ((uint64_t)op[4] * 1842640072758191403UL) + ((uint64_t)op[5] * 9819296124003355144UL) + ((uint64_t)op[6] * 16359183155090143332UL) + ((uint64_t)op[7] * 17658620449639048605UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 111096538306537L) - ((((int128)tmp_q[1] * 86012065123405L) + ((int128)tmp_q[2] * 43057064728199L) - ((int128)tmp_q[3] * 121568310341231L) - ((int128)tmp_q[4] * 97173856657628L) - ((int128)tmp_q[5] * 103756214773535L) - ((int128)tmp_q[6] * 173397557039726L) - ((int128)tmp_q[7] * 77845180647202L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 77845180647202L) - ((int128)tmp_q[1] * 111096538306537L) - ((((int128)tmp_q[2] * 86012065123405L) + ((int128)tmp_q[3] * 43057064728199L) - ((int128)tmp_q[4] * 121568310341231L) - ((int128)tmp_q[5] * 97173856657628L) - ((int128)tmp_q[6] * 103756214773535L) - ((int128)tmp_q[7] * 173397557039726L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 173397557039726L) - ((int128)tmp_q[1] * 77845180647202L) - ((int128)tmp_q[2] * 111096538306537L) - ((((int128)tmp_q[3] * 86012065123405L) + ((int128)tmp_q[4] * 43057064728199L) - ((int128)tmp_q[5] * 121568310341231L) - ((int128)tmp_q[6] * 97173856657628L) - ((int128)tmp_q[7] * 103756214773535L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 103756214773535L) - ((int128)tmp_q[1] * 173397557039726L) - ((int128)tmp_q[2] * 77845180647202L) - ((int128)tmp_q[3] * 111096538306537L) - ((((int128)tmp_q[4] * 86012065123405L) + ((int128)tmp_q[5] * 43057064728199L) - ((int128)tmp_q[6] * 121568310341231L) - ((int128)tmp_q[7] * 97173856657628L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 97173856657628L) - ((int128)tmp_q[1] * 103756214773535L) - ((int128)tmp_q[2] * 173397557039726L) - ((int128)tmp_q[3] * 77845180647202L) - ((int128)tmp_q[4] * 111096538306537L) - ((((int128)tmp_q[5] * 86012065123405L) + ((int128)tmp_q[6] * 43057064728199L) - ((int128)tmp_q[7] * 121568310341231L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 121568310341231L) - ((int128)tmp_q[1] * 97173856657628L) - ((int128)tmp_q[2] * 103756214773535L) - ((int128)tmp_q[3] * 173397557039726L) - ((int128)tmp_q[4] * 77845180647202L) - ((int128)tmp_q[5] * 111096538306537L) - ((((int128)tmp_q[6] * 86012065123405L) + ((int128)tmp_q[7] * 43057064728199L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 43057064728199L) - ((int128)tmp_q[1] * 121568310341231L) - ((int128)tmp_q[2] * 97173856657628L) - ((int128)tmp_q[3] * 103756214773535L) - ((int128)tmp_q[4] * 173397557039726L) - ((int128)tmp_q[5] * 77845180647202L) - ((int128)tmp_q[6] * 111096538306537L) - ((int128)tmp_q[7] * 516072390740430L);
	tmp_zero[7] = ((int128)tmp_q[0] * 86012065123405L) + ((int128)tmp_q[1] * 43057064728199L) - ((int128)tmp_q[2] * 121568310341231L) - ((int128)tmp_q[3] * 97173856657628L) - ((int128)tmp_q[4] * 103756214773535L) - ((int128)tmp_q[5] * 173397557039726L) - ((int128)tmp_q[6] * 77845180647202L) - ((int128)tmp_q[7] * 111096538306537L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

