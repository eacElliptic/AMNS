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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 966825271350861258UL) + ((((uint64_t)op[1] * 10422173979270161457UL) + ((uint64_t)op[2] * 1959723463719813945UL) + ((uint64_t)op[3] * 40448040033536837UL) + ((uint64_t)op[4] * 15982502796102554461UL) + ((uint64_t)op[5] * 5835019411523540573UL) + ((uint64_t)op[6] * 5836919486189460400UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 5836919486189460400UL) + ((uint64_t)op[1] * 966825271350861258UL) + ((((uint64_t)op[2] * 10422173979270161457UL) + ((uint64_t)op[3] * 1959723463719813945UL) + ((uint64_t)op[4] * 40448040033536837UL) + ((uint64_t)op[5] * 15982502796102554461UL) + ((uint64_t)op[6] * 5835019411523540573UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 5835019411523540573UL) + ((uint64_t)op[1] * 5836919486189460400UL) + ((uint64_t)op[2] * 966825271350861258UL) + ((((uint64_t)op[3] * 10422173979270161457UL) + ((uint64_t)op[4] * 1959723463719813945UL) + ((uint64_t)op[5] * 40448040033536837UL) + ((uint64_t)op[6] * 15982502796102554461UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 15982502796102554461UL) + ((uint64_t)op[1] * 5835019411523540573UL) + ((uint64_t)op[2] * 5836919486189460400UL) + ((uint64_t)op[3] * 966825271350861258UL) + ((((uint64_t)op[4] * 10422173979270161457UL) + ((uint64_t)op[5] * 1959723463719813945UL) + ((uint64_t)op[6] * 40448040033536837UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 40448040033536837UL) + ((uint64_t)op[1] * 15982502796102554461UL) + ((uint64_t)op[2] * 5835019411523540573UL) + ((uint64_t)op[3] * 5836919486189460400UL) + ((uint64_t)op[4] * 966825271350861258UL) + ((((uint64_t)op[5] * 10422173979270161457UL) + ((uint64_t)op[6] * 1959723463719813945UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 1959723463719813945UL) + ((uint64_t)op[1] * 40448040033536837UL) + ((uint64_t)op[2] * 15982502796102554461UL) + ((uint64_t)op[3] * 5835019411523540573UL) + ((uint64_t)op[4] * 5836919486189460400UL) + ((uint64_t)op[5] * 966825271350861258UL) + ((uint64_t)op[6] * 15217381748931704053UL);
	tmp_q[6] = ((uint64_t)op[0] * 10422173979270161457UL) + ((uint64_t)op[1] * 1959723463719813945UL) + ((uint64_t)op[2] * 40448040033536837UL) + ((uint64_t)op[3] * 15982502796102554461UL) + ((uint64_t)op[4] * 5835019411523540573UL) + ((uint64_t)op[5] * 5836919486189460400UL) + ((uint64_t)op[6] * 966825271350861258UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 53832001454L) + ((-((int128)tmp_q[1] * 53860333424L) + ((int128)tmp_q[2] * 46155887765L) + ((int128)tmp_q[3] * 13929704214L) - ((int128)tmp_q[4] * 68730651329L) - ((int128)tmp_q[5] * 8262570736L) + ((int128)tmp_q[6] * 2960700277L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 2960700277L) - ((int128)tmp_q[1] * 53832001454L) + ((-((int128)tmp_q[2] * 53860333424L) + ((int128)tmp_q[3] * 46155887765L) + ((int128)tmp_q[4] * 13929704214L) - ((int128)tmp_q[5] * 68730651329L) - ((int128)tmp_q[6] * 8262570736L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 8262570736L) + ((int128)tmp_q[1] * 2960700277L) - ((int128)tmp_q[2] * 53832001454L) + ((-((int128)tmp_q[3] * 53860333424L) + ((int128)tmp_q[4] * 46155887765L) + ((int128)tmp_q[5] * 13929704214L) - ((int128)tmp_q[6] * 68730651329L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 68730651329L) - ((int128)tmp_q[1] * 8262570736L) + ((int128)tmp_q[2] * 2960700277L) - ((int128)tmp_q[3] * 53832001454L) + ((-((int128)tmp_q[4] * 53860333424L) + ((int128)tmp_q[5] * 46155887765L) + ((int128)tmp_q[6] * 13929704214L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 13929704214L) - ((int128)tmp_q[1] * 68730651329L) - ((int128)tmp_q[2] * 8262570736L) + ((int128)tmp_q[3] * 2960700277L) - ((int128)tmp_q[4] * 53832001454L) + ((-((int128)tmp_q[5] * 53860333424L) + ((int128)tmp_q[6] * 46155887765L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 46155887765L) + ((int128)tmp_q[1] * 13929704214L) - ((int128)tmp_q[2] * 68730651329L) - ((int128)tmp_q[3] * 8262570736L) + ((int128)tmp_q[4] * 2960700277L) - ((int128)tmp_q[5] * 53832001454L) - ((int128)tmp_q[6] * 269301667120L);
	tmp_zero[6] = -((int128)tmp_q[0] * 53860333424L) + ((int128)tmp_q[1] * 46155887765L) + ((int128)tmp_q[2] * 13929704214L) - ((int128)tmp_q[3] * 68730651329L) - ((int128)tmp_q[4] * 8262570736L) + ((int128)tmp_q[5] * 2960700277L) - ((int128)tmp_q[6] * 53832001454L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

