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
	tmp_q[0] = ((uint64_t)op[0] * 7926503085318079461UL) + ((((uint64_t)op[1] * 13359108520798040496UL) + ((uint64_t)op[2] * 4436675563331104775UL) + ((uint64_t)op[3] * 2391540371389154979UL) + ((uint64_t)op[4] * 7310819339796767594UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 7310819339796767594UL) + ((uint64_t)op[1] * 7926503085318079461UL) + ((((uint64_t)op[2] * 13359108520798040496UL) + ((uint64_t)op[3] * 4436675563331104775UL) + ((uint64_t)op[4] * 2391540371389154979UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 2391540371389154979UL) + ((uint64_t)op[1] * 7310819339796767594UL) + ((uint64_t)op[2] * 7926503085318079461UL) + ((((uint64_t)op[3] * 13359108520798040496UL) + ((uint64_t)op[4] * 4436675563331104775UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 4436675563331104775UL) + ((uint64_t)op[1] * 2391540371389154979UL) + ((uint64_t)op[2] * 7310819339796767594UL) + ((uint64_t)op[3] * 7926503085318079461UL) + ((uint64_t)op[4] * 11455310382861547632UL);
	tmp_q[4] = ((uint64_t)op[0] * 13359108520798040496UL) + ((uint64_t)op[1] * 4436675563331104775UL) + ((uint64_t)op[2] * 2391540371389154979UL) + ((uint64_t)op[3] * 7310819339796767594UL) + ((uint64_t)op[4] * 7926503085318079461UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 243772644441L) + ((((int128)tmp_q[1] * 211395158063L) + ((int128)tmp_q[2] * 106978645498L) + ((int128)tmp_q[3] * 48955900232L) + ((int128)tmp_q[4] * 22157008559L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 22157008559L) - ((int128)tmp_q[1] * 243772644441L) + ((((int128)tmp_q[2] * 211395158063L) + ((int128)tmp_q[3] * 106978645498L) + ((int128)tmp_q[4] * 48955900232L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 48955900232L) + ((int128)tmp_q[1] * 22157008559L) - ((int128)tmp_q[2] * 243772644441L) + ((((int128)tmp_q[3] * 211395158063L) + ((int128)tmp_q[4] * 106978645498L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 106978645498L) + ((int128)tmp_q[1] * 48955900232L) + ((int128)tmp_q[2] * 22157008559L) - ((int128)tmp_q[3] * 243772644441L) + ((int128)tmp_q[4] * 1056975790315L);
	tmp_zero[4] = ((int128)tmp_q[0] * 211395158063L) + ((int128)tmp_q[1] * 106978645498L) + ((int128)tmp_q[2] * 48955900232L) + ((int128)tmp_q[3] * 22157008559L) - ((int128)tmp_q[4] * 243772644441L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

