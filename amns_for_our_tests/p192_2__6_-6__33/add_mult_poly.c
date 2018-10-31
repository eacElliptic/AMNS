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
	tmp_q[0] = ((uint64_t)op[0] * 1958471341278944467UL) + ((((uint64_t)op[1] * 14950517439067613568UL) + ((uint64_t)op[2] * 1132898623760754516UL) + ((uint64_t)op[3] * 13916977056006193507UL) + ((uint64_t)op[4] * 1816478208678525164UL) + ((uint64_t)op[5] * 16047624936612213194UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16047624936612213194UL) + ((uint64_t)op[1] * 1958471341278944467UL) + ((((uint64_t)op[2] * 14950517439067613568UL) + ((uint64_t)op[3] * 1132898623760754516UL) + ((uint64_t)op[4] * 13916977056006193507UL) + ((uint64_t)op[5] * 1816478208678525164UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 1816478208678525164UL) + ((uint64_t)op[1] * 16047624936612213194UL) + ((uint64_t)op[2] * 1958471341278944467UL) + ((((uint64_t)op[3] * 14950517439067613568UL) + ((uint64_t)op[4] * 1132898623760754516UL) + ((uint64_t)op[5] * 13916977056006193507UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 13916977056006193507UL) + ((uint64_t)op[1] * 1816478208678525164UL) + ((uint64_t)op[2] * 16047624936612213194UL) + ((uint64_t)op[3] * 1958471341278944467UL) + ((((uint64_t)op[4] * 14950517439067613568UL) + ((uint64_t)op[5] * 1132898623760754516UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 1132898623760754516UL) + ((uint64_t)op[1] * 13916977056006193507UL) + ((uint64_t)op[2] * 1816478208678525164UL) + ((uint64_t)op[3] * 16047624936612213194UL) + ((uint64_t)op[4] * 1958471341278944467UL) + ((uint64_t)op[5] * 2530615734142076672UL);
	tmp_q[5] = ((uint64_t)op[0] * 14950517439067613568UL) + ((uint64_t)op[1] * 1132898623760754516UL) + ((uint64_t)op[2] * 13916977056006193507UL) + ((uint64_t)op[3] * 1816478208678525164UL) + ((uint64_t)op[4] * 16047624936612213194UL) + ((uint64_t)op[5] * 1958471341278944467UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 400688627L) - ((((int128)tmp_q[1] * 280522044L) - ((int128)tmp_q[2] * 2928953384L) + ((int128)tmp_q[3] * 1481723605L) + ((int128)tmp_q[4] * 2565863352L) - ((int128)tmp_q[5] * 686551186L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 686551186L) + ((int128)tmp_q[1] * 400688627L) - ((((int128)tmp_q[2] * 280522044L) - ((int128)tmp_q[3] * 2928953384L) + ((int128)tmp_q[4] * 1481723605L) + ((int128)tmp_q[5] * 2565863352L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 2565863352L) - ((int128)tmp_q[1] * 686551186L) + ((int128)tmp_q[2] * 400688627L) - ((((int128)tmp_q[3] * 280522044L) - ((int128)tmp_q[4] * 2928953384L) + ((int128)tmp_q[5] * 1481723605L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1481723605L) + ((int128)tmp_q[1] * 2565863352L) - ((int128)tmp_q[2] * 686551186L) + ((int128)tmp_q[3] * 400688627L) - ((((int128)tmp_q[4] * 280522044L) - ((int128)tmp_q[5] * 2928953384L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 2928953384L) + ((int128)tmp_q[1] * 1481723605L) + ((int128)tmp_q[2] * 2565863352L) - ((int128)tmp_q[3] * 686551186L) + ((int128)tmp_q[4] * 400688627L) - ((int128)tmp_q[5] * 1683132264L);
	tmp_zero[5] = ((int128)tmp_q[0] * 280522044L) - ((int128)tmp_q[1] * 2928953384L) + ((int128)tmp_q[2] * 1481723605L) + ((int128)tmp_q[3] * 2565863352L) - ((int128)tmp_q[4] * 686551186L) + ((int128)tmp_q[5] * 400688627L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

