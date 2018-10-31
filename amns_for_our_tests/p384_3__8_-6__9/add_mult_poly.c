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
	tmp_q[0] = ((uint64_t)op[0] * 9198394018944246749UL) + ((((uint64_t)op[1] * 17836239075739092487UL) + ((uint64_t)op[2] * 5135557155526424713UL) + ((uint64_t)op[3] * 4299112259328869557UL) + ((uint64_t)op[4] * 11302597151422636938UL) + ((uint64_t)op[5] * 18026991668138535797UL) + ((uint64_t)op[6] * 6170616996543674640UL) + ((uint64_t)op[7] * 14696779666980174747UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 14696779666980174747UL) + ((uint64_t)op[1] * 9198394018944246749UL) + ((((uint64_t)op[2] * 17836239075739092487UL) + ((uint64_t)op[3] * 5135557155526424713UL) + ((uint64_t)op[4] * 4299112259328869557UL) + ((uint64_t)op[5] * 11302597151422636938UL) + ((uint64_t)op[6] * 18026991668138535797UL) + ((uint64_t)op[7] * 6170616996543674640UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 6170616996543674640UL) + ((uint64_t)op[1] * 14696779666980174747UL) + ((uint64_t)op[2] * 9198394018944246749UL) + ((((uint64_t)op[3] * 17836239075739092487UL) + ((uint64_t)op[4] * 5135557155526424713UL) + ((uint64_t)op[5] * 4299112259328869557UL) + ((uint64_t)op[6] * 11302597151422636938UL) + ((uint64_t)op[7] * 18026991668138535797UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 18026991668138535797UL) + ((uint64_t)op[1] * 6170616996543674640UL) + ((uint64_t)op[2] * 14696779666980174747UL) + ((uint64_t)op[3] * 9198394018944246749UL) + ((((uint64_t)op[4] * 17836239075739092487UL) + ((uint64_t)op[5] * 5135557155526424713UL) + ((uint64_t)op[6] * 4299112259328869557UL) + ((uint64_t)op[7] * 11302597151422636938UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 11302597151422636938UL) + ((uint64_t)op[1] * 18026991668138535797UL) + ((uint64_t)op[2] * 6170616996543674640UL) + ((uint64_t)op[3] * 14696779666980174747UL) + ((uint64_t)op[4] * 9198394018944246749UL) + ((((uint64_t)op[5] * 17836239075739092487UL) + ((uint64_t)op[6] * 5135557155526424713UL) + ((uint64_t)op[7] * 4299112259328869557UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 4299112259328869557UL) + ((uint64_t)op[1] * 11302597151422636938UL) + ((uint64_t)op[2] * 18026991668138535797UL) + ((uint64_t)op[3] * 6170616996543674640UL) + ((uint64_t)op[4] * 14696779666980174747UL) + ((uint64_t)op[5] * 9198394018944246749UL) + ((((uint64_t)op[6] * 17836239075739092487UL) + ((uint64_t)op[7] * 5135557155526424713UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 5135557155526424713UL) + ((uint64_t)op[1] * 4299112259328869557UL) + ((uint64_t)op[2] * 11302597151422636938UL) + ((uint64_t)op[3] * 18026991668138535797UL) + ((uint64_t)op[4] * 6170616996543674640UL) + ((uint64_t)op[5] * 14696779666980174747UL) + ((uint64_t)op[6] * 9198394018944246749UL) + ((uint64_t)op[7] * 3663029987822754774UL);
	tmp_q[7] = ((uint64_t)op[0] * 17836239075739092487UL) + ((uint64_t)op[1] * 5135557155526424713UL) + ((uint64_t)op[2] * 4299112259328869557UL) + ((uint64_t)op[3] * 11302597151422636938UL) + ((uint64_t)op[4] * 18026991668138535797UL) + ((uint64_t)op[5] * 6170616996543674640UL) + ((uint64_t)op[6] * 14696779666980174747UL) + ((uint64_t)op[7] * 9198394018944246749UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 110670459176425L) - ((((int128)tmp_q[1] * 13973312924425L) - ((int128)tmp_q[2] * 21320003533977L) - ((int128)tmp_q[3] * 29789980903975L) - ((int128)tmp_q[4] * 144029163700891L) + ((int128)tmp_q[5] * 137787058787894L) - ((int128)tmp_q[6] * 100273142494801L) - ((int128)tmp_q[7] * 38626406850043L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 38626406850043L) - ((int128)tmp_q[1] * 110670459176425L) - ((((int128)tmp_q[2] * 13973312924425L) - ((int128)tmp_q[3] * 21320003533977L) - ((int128)tmp_q[4] * 29789980903975L) - ((int128)tmp_q[5] * 144029163700891L) + ((int128)tmp_q[6] * 137787058787894L) - ((int128)tmp_q[7] * 100273142494801L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 100273142494801L) - ((int128)tmp_q[1] * 38626406850043L) - ((int128)tmp_q[2] * 110670459176425L) - ((((int128)tmp_q[3] * 13973312924425L) - ((int128)tmp_q[4] * 21320003533977L) - ((int128)tmp_q[5] * 29789980903975L) - ((int128)tmp_q[6] * 144029163700891L) + ((int128)tmp_q[7] * 137787058787894L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 137787058787894L) - ((int128)tmp_q[1] * 100273142494801L) - ((int128)tmp_q[2] * 38626406850043L) - ((int128)tmp_q[3] * 110670459176425L) - ((((int128)tmp_q[4] * 13973312924425L) - ((int128)tmp_q[5] * 21320003533977L) - ((int128)tmp_q[6] * 29789980903975L) - ((int128)tmp_q[7] * 144029163700891L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 144029163700891L) + ((int128)tmp_q[1] * 137787058787894L) - ((int128)tmp_q[2] * 100273142494801L) - ((int128)tmp_q[3] * 38626406850043L) - ((int128)tmp_q[4] * 110670459176425L) - ((((int128)tmp_q[5] * 13973312924425L) - ((int128)tmp_q[6] * 21320003533977L) - ((int128)tmp_q[7] * 29789980903975L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 29789980903975L) - ((int128)tmp_q[1] * 144029163700891L) + ((int128)tmp_q[2] * 137787058787894L) - ((int128)tmp_q[3] * 100273142494801L) - ((int128)tmp_q[4] * 38626406850043L) - ((int128)tmp_q[5] * 110670459176425L) - ((((int128)tmp_q[6] * 13973312924425L) - ((int128)tmp_q[7] * 21320003533977L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 21320003533977L) - ((int128)tmp_q[1] * 29789980903975L) - ((int128)tmp_q[2] * 144029163700891L) + ((int128)tmp_q[3] * 137787058787894L) - ((int128)tmp_q[4] * 100273142494801L) - ((int128)tmp_q[5] * 38626406850043L) - ((int128)tmp_q[6] * 110670459176425L) - ((int128)tmp_q[7] * 83839877546550L);
	tmp_zero[7] = ((int128)tmp_q[0] * 13973312924425L) - ((int128)tmp_q[1] * 21320003533977L) - ((int128)tmp_q[2] * 29789980903975L) - ((int128)tmp_q[3] * 144029163700891L) + ((int128)tmp_q[4] * 137787058787894L) - ((int128)tmp_q[5] * 100273142494801L) - ((int128)tmp_q[6] * 38626406850043L) - ((int128)tmp_q[7] * 110670459176425L);

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

