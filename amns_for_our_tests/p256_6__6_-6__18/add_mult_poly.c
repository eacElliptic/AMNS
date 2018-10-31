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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15003863097387989097UL) + ((((uint64_t)op[1] * 16101278054621039687UL) + ((uint64_t)op[2] * 16164434102061198870UL) + ((uint64_t)op[3] * 15337498197876698615UL) + ((uint64_t)op[4] * 393685086721435272UL) + ((uint64_t)op[5] * 12457753832409935727UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 12457753832409935727UL) + ((uint64_t)op[1] * 15003863097387989097UL) + ((((uint64_t)op[2] * 16101278054621039687UL) + ((uint64_t)op[3] * 16164434102061198870UL) + ((uint64_t)op[4] * 15337498197876698615UL) + ((uint64_t)op[5] * 393685086721435272UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 393685086721435272UL) + ((uint64_t)op[1] * 12457753832409935727UL) + ((uint64_t)op[2] * 15003863097387989097UL) + ((((uint64_t)op[3] * 16101278054621039687UL) + ((uint64_t)op[4] * 16164434102061198870UL) + ((uint64_t)op[5] * 15337498197876698615UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 15337498197876698615UL) + ((uint64_t)op[1] * 393685086721435272UL) + ((uint64_t)op[2] * 12457753832409935727UL) + ((uint64_t)op[3] * 15003863097387989097UL) + ((((uint64_t)op[4] * 16101278054621039687UL) + ((uint64_t)op[5] * 16164434102061198870UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 16164434102061198870UL) + ((uint64_t)op[1] * 15337498197876698615UL) + ((uint64_t)op[2] * 393685086721435272UL) + ((uint64_t)op[3] * 12457753832409935727UL) + ((uint64_t)op[4] * 15003863097387989097UL) + ((uint64_t)op[5] * 14072796114531071574UL);
	tmp_q[5] = ((uint64_t)op[0] * 16101278054621039687UL) + ((uint64_t)op[1] * 16164434102061198870UL) + ((uint64_t)op[2] * 15337498197876698615UL) + ((uint64_t)op[3] * 393685086721435272UL) + ((uint64_t)op[4] * 12457753832409935727UL) + ((uint64_t)op[5] * 15003863097387989097UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1309044636845L) - ((-((int128)tmp_q[1] * 4815380636465L) - ((int128)tmp_q[2] * 1328959036295L) - ((int128)tmp_q[3] * 1504186266242L) + ((int128)tmp_q[4] * 2939497897813L) - ((int128)tmp_q[5] * 305820628961L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 305820628961L) - ((int128)tmp_q[1] * 1309044636845L) - ((-((int128)tmp_q[2] * 4815380636465L) - ((int128)tmp_q[3] * 1328959036295L) - ((int128)tmp_q[4] * 1504186266242L) + ((int128)tmp_q[5] * 2939497897813L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 2939497897813L) - ((int128)tmp_q[1] * 305820628961L) - ((int128)tmp_q[2] * 1309044636845L) - ((-((int128)tmp_q[3] * 4815380636465L) - ((int128)tmp_q[4] * 1328959036295L) - ((int128)tmp_q[5] * 1504186266242L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 1504186266242L) + ((int128)tmp_q[1] * 2939497897813L) - ((int128)tmp_q[2] * 305820628961L) - ((int128)tmp_q[3] * 1309044636845L) - ((-((int128)tmp_q[4] * 4815380636465L) - ((int128)tmp_q[5] * 1328959036295L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 1328959036295L) - ((int128)tmp_q[1] * 1504186266242L) + ((int128)tmp_q[2] * 2939497897813L) - ((int128)tmp_q[3] * 305820628961L) - ((int128)tmp_q[4] * 1309044636845L) + ((int128)tmp_q[5] * 28892283818790L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4815380636465L) - ((int128)tmp_q[1] * 1328959036295L) - ((int128)tmp_q[2] * 1504186266242L) + ((int128)tmp_q[3] * 2939497897813L) - ((int128)tmp_q[4] * 305820628961L) - ((int128)tmp_q[5] * 1309044636845L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

