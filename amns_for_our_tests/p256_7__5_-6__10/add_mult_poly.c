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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1657262359104878967UL) + ((((uint64_t)op[1] * 12801301339534847489UL) + ((uint64_t)op[2] * 12335257886762897653UL) + ((uint64_t)op[3] * 16495131515819765716UL) + ((uint64_t)op[4] * 17309770669706671342UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 17309770669706671342UL) + ((uint64_t)op[1] * 1657262359104878967UL) + ((((uint64_t)op[2] * 12801301339534847489UL) + ((uint64_t)op[3] * 12335257886762897653UL) + ((uint64_t)op[4] * 16495131515819765716UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 16495131515819765716UL) + ((uint64_t)op[1] * 17309770669706671342UL) + ((uint64_t)op[2] * 1657262359104878967UL) + ((((uint64_t)op[3] * 12801301339534847489UL) + ((uint64_t)op[4] * 12335257886762897653UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 12335257886762897653UL) + ((uint64_t)op[1] * 16495131515819765716UL) + ((uint64_t)op[2] * 17309770669706671342UL) + ((uint64_t)op[3] * 1657262359104878967UL) + ((uint64_t)op[4] * 15425912331338673146UL);
	tmp_q[4] = ((uint64_t)op[0] * 12801301339534847489UL) + ((uint64_t)op[1] * 12335257886762897653UL) + ((uint64_t)op[2] * 16495131515819765716UL) + ((uint64_t)op[3] * 17309770669706671342UL) + ((uint64_t)op[4] * 1657262359104878967UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 285588126399395L) - ((((int128)tmp_q[1] * 92924448348003L) - ((int128)tmp_q[2] * 819177639074025L) + ((int128)tmp_q[3] * 267463280423080L) + ((int128)tmp_q[4] * 526238486593324L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 526238486593324L) - ((int128)tmp_q[1] * 285588126399395L) - ((((int128)tmp_q[2] * 92924448348003L) - ((int128)tmp_q[3] * 819177639074025L) + ((int128)tmp_q[4] * 267463280423080L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 267463280423080L) + ((int128)tmp_q[1] * 526238486593324L) - ((int128)tmp_q[2] * 285588126399395L) - ((((int128)tmp_q[3] * 92924448348003L) - ((int128)tmp_q[4] * 819177639074025L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 819177639074025L) + ((int128)tmp_q[1] * 267463280423080L) + ((int128)tmp_q[2] * 526238486593324L) - ((int128)tmp_q[3] * 285588126399395L) - ((int128)tmp_q[4] * 557546690088018L);
	tmp_zero[4] = ((int128)tmp_q[0] * 92924448348003L) - ((int128)tmp_q[1] * 819177639074025L) + ((int128)tmp_q[2] * 267463280423080L) + ((int128)tmp_q[3] * 526238486593324L) - ((int128)tmp_q[4] * 285588126399395L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

