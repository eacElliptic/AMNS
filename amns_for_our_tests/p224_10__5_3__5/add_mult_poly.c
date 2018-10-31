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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14598551433453741301UL) + ((((uint64_t)op[1] * 4170843000761427489UL) + ((uint64_t)op[2] * 16332509976683406973UL) + ((uint64_t)op[3] * 10496868331823263852UL) + ((uint64_t)op[4] * 7875968826033548924UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 7875968826033548924UL) + ((uint64_t)op[1] * 14598551433453741301UL) + ((((uint64_t)op[2] * 4170843000761427489UL) + ((uint64_t)op[3] * 16332509976683406973UL) + ((uint64_t)op[4] * 10496868331823263852UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 10496868331823263852UL) + ((uint64_t)op[1] * 7875968826033548924UL) + ((uint64_t)op[2] * 14598551433453741301UL) + ((((uint64_t)op[3] * 4170843000761427489UL) + ((uint64_t)op[4] * 16332509976683406973UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 16332509976683406973UL) + ((uint64_t)op[1] * 10496868331823263852UL) + ((uint64_t)op[2] * 7875968826033548924UL) + ((uint64_t)op[3] * 14598551433453741301UL) + ((uint64_t)op[4] * 12512529002284282467UL);
	tmp_q[4] = ((uint64_t)op[0] * 4170843000761427489UL) + ((uint64_t)op[1] * 16332509976683406973UL) + ((uint64_t)op[2] * 10496868331823263852UL) + ((uint64_t)op[3] * 7875968826033548924UL) + ((uint64_t)op[4] * 14598551433453741301UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 574185854704L) + ((((int128)tmp_q[1] * 21364086778189L) + ((int128)tmp_q[2] * 18324256418073L) - ((int128)tmp_q[3] * 944798473906L) - ((int128)tmp_q[4] * 2148609385961L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 2148609385961L) - ((int128)tmp_q[1] * 574185854704L) + ((((int128)tmp_q[2] * 21364086778189L) + ((int128)tmp_q[3] * 18324256418073L) - ((int128)tmp_q[4] * 944798473906L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 944798473906L) - ((int128)tmp_q[1] * 2148609385961L) - ((int128)tmp_q[2] * 574185854704L) + ((((int128)tmp_q[3] * 21364086778189L) + ((int128)tmp_q[4] * 18324256418073L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 18324256418073L) - ((int128)tmp_q[1] * 944798473906L) - ((int128)tmp_q[2] * 2148609385961L) - ((int128)tmp_q[3] * 574185854704L) + ((int128)tmp_q[4] * 64092260334567L);
	tmp_zero[4] = ((int128)tmp_q[0] * 21364086778189L) + ((int128)tmp_q[1] * 18324256418073L) - ((int128)tmp_q[2] * 944798473906L) - ((int128)tmp_q[3] * 2148609385961L) - ((int128)tmp_q[4] * 574185854704L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

