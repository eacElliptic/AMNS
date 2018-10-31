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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18240263364850266193UL) + ((((uint64_t)op[1] * 9593591566218194774UL) + ((uint64_t)op[2] * 12912761266374106896UL) + ((uint64_t)op[3] * 162061107381468571UL) + ((uint64_t)op[4] * 5153996530699554606UL) + ((uint64_t)op[5] * 5067123049366437342UL) + ((uint64_t)op[6] * 11754358658287707833UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 11754358658287707833UL) + ((uint64_t)op[1] * 18240263364850266193UL) + ((((uint64_t)op[2] * 9593591566218194774UL) + ((uint64_t)op[3] * 12912761266374106896UL) + ((uint64_t)op[4] * 162061107381468571UL) + ((uint64_t)op[5] * 5153996530699554606UL) + ((uint64_t)op[6] * 5067123049366437342UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 5067123049366437342UL) + ((uint64_t)op[1] * 11754358658287707833UL) + ((uint64_t)op[2] * 18240263364850266193UL) + ((((uint64_t)op[3] * 9593591566218194774UL) + ((uint64_t)op[4] * 12912761266374106896UL) + ((uint64_t)op[5] * 162061107381468571UL) + ((uint64_t)op[6] * 5153996530699554606UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 5153996530699554606UL) + ((uint64_t)op[1] * 5067123049366437342UL) + ((uint64_t)op[2] * 11754358658287707833UL) + ((uint64_t)op[3] * 18240263364850266193UL) + ((((uint64_t)op[4] * 9593591566218194774UL) + ((uint64_t)op[5] * 12912761266374106896UL) + ((uint64_t)op[6] * 162061107381468571UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 162061107381468571UL) + ((uint64_t)op[1] * 5153996530699554606UL) + ((uint64_t)op[2] * 5067123049366437342UL) + ((uint64_t)op[3] * 11754358658287707833UL) + ((uint64_t)op[4] * 18240263364850266193UL) + ((((uint64_t)op[5] * 9593591566218194774UL) + ((uint64_t)op[6] * 12912761266374106896UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 12912761266374106896UL) + ((uint64_t)op[1] * 162061107381468571UL) + ((uint64_t)op[2] * 5153996530699554606UL) + ((uint64_t)op[3] * 5067123049366437342UL) + ((uint64_t)op[4] * 11754358658287707833UL) + ((uint64_t)op[5] * 18240263364850266193UL) + ((uint64_t)op[6] * 17706305014982713684UL);
	tmp_q[6] = ((uint64_t)op[0] * 9593591566218194774UL) + ((uint64_t)op[1] * 12912761266374106896UL) + ((uint64_t)op[2] * 162061107381468571UL) + ((uint64_t)op[3] * 5153996530699554606UL) + ((uint64_t)op[4] * 5067123049366437342UL) + ((uint64_t)op[5] * 11754358658287707833UL) + ((uint64_t)op[6] * 18240263364850266193UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 19666876112725529L) - ((-((int128)tmp_q[1] * 6076871028556700L) + ((int128)tmp_q[2] * 16398784503070493L) + ((int128)tmp_q[3] * 12232861060948242L) + ((int128)tmp_q[4] * 12991932905546343L) - ((int128)tmp_q[5] * 5305659222713911L) - ((int128)tmp_q[6] * 6266450507518573L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 6266450507518573L) + ((int128)tmp_q[1] * 19666876112725529L) - ((-((int128)tmp_q[2] * 6076871028556700L) + ((int128)tmp_q[3] * 16398784503070493L) + ((int128)tmp_q[4] * 12232861060948242L) + ((int128)tmp_q[5] * 12991932905546343L) - ((int128)tmp_q[6] * 5305659222713911L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 5305659222713911L) - ((int128)tmp_q[1] * 6266450507518573L) + ((int128)tmp_q[2] * 19666876112725529L) - ((-((int128)tmp_q[3] * 6076871028556700L) + ((int128)tmp_q[4] * 16398784503070493L) + ((int128)tmp_q[5] * 12232861060948242L) + ((int128)tmp_q[6] * 12991932905546343L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 12991932905546343L) - ((int128)tmp_q[1] * 5305659222713911L) - ((int128)tmp_q[2] * 6266450507518573L) + ((int128)tmp_q[3] * 19666876112725529L) - ((-((int128)tmp_q[4] * 6076871028556700L) + ((int128)tmp_q[5] * 16398784503070493L) + ((int128)tmp_q[6] * 12232861060948242L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 12232861060948242L) + ((int128)tmp_q[1] * 12991932905546343L) - ((int128)tmp_q[2] * 5305659222713911L) - ((int128)tmp_q[3] * 6266450507518573L) + ((int128)tmp_q[4] * 19666876112725529L) - ((-((int128)tmp_q[5] * 6076871028556700L) + ((int128)tmp_q[6] * 16398784503070493L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 16398784503070493L) + ((int128)tmp_q[1] * 12232861060948242L) + ((int128)tmp_q[2] * 12991932905546343L) - ((int128)tmp_q[3] * 5305659222713911L) - ((int128)tmp_q[4] * 6266450507518573L) + ((int128)tmp_q[5] * 19666876112725529L) + ((int128)tmp_q[6] * 12153742057113400L);
	tmp_zero[6] = -((int128)tmp_q[0] * 6076871028556700L) + ((int128)tmp_q[1] * 16398784503070493L) + ((int128)tmp_q[2] * 12232861060948242L) + ((int128)tmp_q[3] * 12991932905546343L) - ((int128)tmp_q[4] * 5305659222713911L) - ((int128)tmp_q[5] * 6266450507518573L) + ((int128)tmp_q[6] * 19666876112725529L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

