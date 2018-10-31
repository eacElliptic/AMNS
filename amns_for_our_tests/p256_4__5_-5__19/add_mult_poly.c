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
	tmp_q[0] = ((uint64_t)op[0] * 11557084361408393936UL) + ((((uint64_t)op[1] * 4591722739246774224UL) + ((uint64_t)op[2] * 2414553948990768442UL) + ((uint64_t)op[3] * 1088102615938305516UL) + ((uint64_t)op[4] * 9673188614892517957UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 9673188614892517957UL) + ((uint64_t)op[1] * 11557084361408393936UL) + ((((uint64_t)op[2] * 4591722739246774224UL) + ((uint64_t)op[3] * 2414553948990768442UL) + ((uint64_t)op[4] * 1088102615938305516UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1088102615938305516UL) + ((uint64_t)op[1] * 9673188614892517957UL) + ((uint64_t)op[2] * 11557084361408393936UL) + ((((uint64_t)op[3] * 4591722739246774224UL) + ((uint64_t)op[4] * 2414553948990768442UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 2414553948990768442UL) + ((uint64_t)op[1] * 1088102615938305516UL) + ((uint64_t)op[2] * 9673188614892517957UL) + ((uint64_t)op[3] * 11557084361408393936UL) + ((uint64_t)op[4] * 13934874451185232112UL);
	tmp_q[4] = ((uint64_t)op[0] * 4591722739246774224UL) + ((uint64_t)op[1] * 2414553948990768442UL) + ((uint64_t)op[2] * 1088102615938305516UL) + ((uint64_t)op[3] * 9673188614892517957UL) + ((uint64_t)op[4] * 11557084361408393936UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 546538529534524L) - ((-((int128)tmp_q[1] * 1498132611765607L) + ((int128)tmp_q[2] * 1519005533038108L) - ((int128)tmp_q[3] * 312851932602832L) - ((int128)tmp_q[4] * 241903426641414L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 241903426641414L) - ((int128)tmp_q[1] * 546538529534524L) - ((-((int128)tmp_q[2] * 1498132611765607L) + ((int128)tmp_q[3] * 1519005533038108L) - ((int128)tmp_q[4] * 312851932602832L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 312851932602832L) - ((int128)tmp_q[1] * 241903426641414L) - ((int128)tmp_q[2] * 546538529534524L) - ((-((int128)tmp_q[3] * 1498132611765607L) + ((int128)tmp_q[4] * 1519005533038108L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1519005533038108L) - ((int128)tmp_q[1] * 312851932602832L) - ((int128)tmp_q[2] * 241903426641414L) - ((int128)tmp_q[3] * 546538529534524L) + ((int128)tmp_q[4] * 7490663058828035L);
	tmp_zero[4] = -((int128)tmp_q[0] * 1498132611765607L) + ((int128)tmp_q[1] * 1519005533038108L) - ((int128)tmp_q[2] * 312851932602832L) - ((int128)tmp_q[3] * 241903426641414L) - ((int128)tmp_q[4] * 546538529534524L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

