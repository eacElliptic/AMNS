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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 164464549250277005UL) + ((((uint64_t)op[1] * 564817265252134471UL) + ((uint64_t)op[2] * 10795913474799505701UL) + ((uint64_t)op[3] * 5574569450421752936UL) + ((uint64_t)op[4] * 3669235315003780175UL) + ((uint64_t)op[5] * 647927251720974164UL) + ((uint64_t)op[6] * 853864020104213801UL) + ((uint64_t)op[7] * 10193507159190974111UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 10193507159190974111UL) + ((uint64_t)op[1] * 164464549250277005UL) + ((((uint64_t)op[2] * 564817265252134471UL) + ((uint64_t)op[3] * 10795913474799505701UL) + ((uint64_t)op[4] * 5574569450421752936UL) + ((uint64_t)op[5] * 3669235315003780175UL) + ((uint64_t)op[6] * 647927251720974164UL) + ((uint64_t)op[7] * 853864020104213801UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 853864020104213801UL) + ((uint64_t)op[1] * 10193507159190974111UL) + ((uint64_t)op[2] * 164464549250277005UL) + ((((uint64_t)op[3] * 564817265252134471UL) + ((uint64_t)op[4] * 10795913474799505701UL) + ((uint64_t)op[5] * 5574569450421752936UL) + ((uint64_t)op[6] * 3669235315003780175UL) + ((uint64_t)op[7] * 647927251720974164UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 647927251720974164UL) + ((uint64_t)op[1] * 853864020104213801UL) + ((uint64_t)op[2] * 10193507159190974111UL) + ((uint64_t)op[3] * 164464549250277005UL) + ((((uint64_t)op[4] * 564817265252134471UL) + ((uint64_t)op[5] * 10795913474799505701UL) + ((uint64_t)op[6] * 5574569450421752936UL) + ((uint64_t)op[7] * 3669235315003780175UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 3669235315003780175UL) + ((uint64_t)op[1] * 647927251720974164UL) + ((uint64_t)op[2] * 853864020104213801UL) + ((uint64_t)op[3] * 10193507159190974111UL) + ((uint64_t)op[4] * 164464549250277005UL) + ((((uint64_t)op[5] * 564817265252134471UL) + ((uint64_t)op[6] * 10795913474799505701UL) + ((uint64_t)op[7] * 5574569450421752936UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 5574569450421752936UL) + ((uint64_t)op[1] * 3669235315003780175UL) + ((uint64_t)op[2] * 647927251720974164UL) + ((uint64_t)op[3] * 853864020104213801UL) + ((uint64_t)op[4] * 10193507159190974111UL) + ((uint64_t)op[5] * 164464549250277005UL) + ((((uint64_t)op[6] * 564817265252134471UL) + ((uint64_t)op[7] * 10795913474799505701UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 10795913474799505701UL) + ((uint64_t)op[1] * 5574569450421752936UL) + ((uint64_t)op[2] * 3669235315003780175UL) + ((uint64_t)op[3] * 647927251720974164UL) + ((uint64_t)op[4] * 853864020104213801UL) + ((uint64_t)op[5] * 10193507159190974111UL) + ((uint64_t)op[6] * 164464549250277005UL) + ((uint64_t)op[7] * 1129634530504268942UL);
	tmp_q[7] = ((uint64_t)op[0] * 564817265252134471UL) + ((uint64_t)op[1] * 10795913474799505701UL) + ((uint64_t)op[2] * 5574569450421752936UL) + ((uint64_t)op[3] * 3669235315003780175UL) + ((uint64_t)op[4] * 647927251720974164UL) + ((uint64_t)op[5] * 853864020104213801UL) + ((uint64_t)op[6] * 10193507159190974111UL) + ((uint64_t)op[7] * 164464549250277005UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 114881566643099L) + ((-((int128)tmp_q[1] * 6121094554234L) + ((int128)tmp_q[2] * 128631044531645L) - ((int128)tmp_q[3] * 61372339356354L) + ((int128)tmp_q[4] * 65839855589808L) - ((int128)tmp_q[5] * 131141844930735L) - ((int128)tmp_q[6] * 14006476165982L) + ((int128)tmp_q[7] * 58742003937649L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 58742003937649L) + ((int128)tmp_q[1] * 114881566643099L) + ((-((int128)tmp_q[2] * 6121094554234L) + ((int128)tmp_q[3] * 128631044531645L) - ((int128)tmp_q[4] * 61372339356354L) + ((int128)tmp_q[5] * 65839855589808L) - ((int128)tmp_q[6] * 131141844930735L) - ((int128)tmp_q[7] * 14006476165982L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 14006476165982L) + ((int128)tmp_q[1] * 58742003937649L) + ((int128)tmp_q[2] * 114881566643099L) + ((-((int128)tmp_q[3] * 6121094554234L) + ((int128)tmp_q[4] * 128631044531645L) - ((int128)tmp_q[5] * 61372339356354L) + ((int128)tmp_q[6] * 65839855589808L) - ((int128)tmp_q[7] * 131141844930735L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 131141844930735L) - ((int128)tmp_q[1] * 14006476165982L) + ((int128)tmp_q[2] * 58742003937649L) + ((int128)tmp_q[3] * 114881566643099L) + ((-((int128)tmp_q[4] * 6121094554234L) + ((int128)tmp_q[5] * 128631044531645L) - ((int128)tmp_q[6] * 61372339356354L) + ((int128)tmp_q[7] * 65839855589808L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 65839855589808L) - ((int128)tmp_q[1] * 131141844930735L) - ((int128)tmp_q[2] * 14006476165982L) + ((int128)tmp_q[3] * 58742003937649L) + ((int128)tmp_q[4] * 114881566643099L) + ((-((int128)tmp_q[5] * 6121094554234L) + ((int128)tmp_q[6] * 128631044531645L) - ((int128)tmp_q[7] * 61372339356354L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 61372339356354L) + ((int128)tmp_q[1] * 65839855589808L) - ((int128)tmp_q[2] * 131141844930735L) - ((int128)tmp_q[3] * 14006476165982L) + ((int128)tmp_q[4] * 58742003937649L) + ((int128)tmp_q[5] * 114881566643099L) + ((-((int128)tmp_q[6] * 6121094554234L) + ((int128)tmp_q[7] * 128631044531645L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 128631044531645L) - ((int128)tmp_q[1] * 61372339356354L) + ((int128)tmp_q[2] * 65839855589808L) - ((int128)tmp_q[3] * 131141844930735L) - ((int128)tmp_q[4] * 14006476165982L) + ((int128)tmp_q[5] * 58742003937649L) + ((int128)tmp_q[6] * 114881566643099L) - ((int128)tmp_q[7] * 12242189108468L);
	tmp_zero[7] = -((int128)tmp_q[0] * 6121094554234L) + ((int128)tmp_q[1] * 128631044531645L) - ((int128)tmp_q[2] * 61372339356354L) + ((int128)tmp_q[3] * 65839855589808L) - ((int128)tmp_q[4] * 131141844930735L) - ((int128)tmp_q[5] * 14006476165982L) + ((int128)tmp_q[6] * 58742003937649L) + ((int128)tmp_q[7] * 114881566643099L);

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

