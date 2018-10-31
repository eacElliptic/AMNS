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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5559193755906915865UL) + ((((uint64_t)op[1] * 15801168389286851165UL) + ((uint64_t)op[2] * 9253439345151828637UL) + ((uint64_t)op[3] * 227132465910925275UL) + ((uint64_t)op[4] * 16577027543877593136UL) + ((uint64_t)op[5] * 10114258534547562201UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 10114258534547562201UL) + ((uint64_t)op[1] * 5559193755906915865UL) + ((((uint64_t)op[2] * 15801168389286851165UL) + ((uint64_t)op[3] * 9253439345151828637UL) + ((uint64_t)op[4] * 227132465910925275UL) + ((uint64_t)op[5] * 16577027543877593136UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 16577027543877593136UL) + ((uint64_t)op[1] * 10114258534547562201UL) + ((uint64_t)op[2] * 5559193755906915865UL) + ((((uint64_t)op[3] * 15801168389286851165UL) + ((uint64_t)op[4] * 9253439345151828637UL) + ((uint64_t)op[5] * 227132465910925275UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 227132465910925275UL) + ((uint64_t)op[1] * 16577027543877593136UL) + ((uint64_t)op[2] * 10114258534547562201UL) + ((uint64_t)op[3] * 5559193755906915865UL) + ((((uint64_t)op[4] * 15801168389286851165UL) + ((uint64_t)op[5] * 9253439345151828637UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 9253439345151828637UL) + ((uint64_t)op[1] * 227132465910925275UL) + ((uint64_t)op[2] * 16577027543877593136UL) + ((uint64_t)op[3] * 10114258534547562201UL) + ((uint64_t)op[4] * 5559193755906915865UL) + ((uint64_t)op[5] * 2573289967173348910UL);
	tmp_q[5] = ((uint64_t)op[0] * 15801168389286851165UL) + ((uint64_t)op[1] * 9253439345151828637UL) + ((uint64_t)op[2] * 227132465910925275UL) + ((uint64_t)op[3] * 16577027543877593136UL) + ((uint64_t)op[4] * 10114258534547562201UL) + ((uint64_t)op[5] * 5559193755906915865UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3205828084867L) + ((-((int128)tmp_q[1] * 741143456125L) + ((int128)tmp_q[2] * 2330359287744L) + ((int128)tmp_q[3] * 101850252210L) - ((int128)tmp_q[4] * 709348580227L) + ((int128)tmp_q[5] * 3736891819937L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 3736891819937L) - ((int128)tmp_q[1] * 3205828084867L) + ((-((int128)tmp_q[2] * 741143456125L) + ((int128)tmp_q[3] * 2330359287744L) + ((int128)tmp_q[4] * 101850252210L) - ((int128)tmp_q[5] * 709348580227L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 709348580227L) + ((int128)tmp_q[1] * 3736891819937L) - ((int128)tmp_q[2] * 3205828084867L) + ((-((int128)tmp_q[3] * 741143456125L) + ((int128)tmp_q[4] * 2330359287744L) + ((int128)tmp_q[5] * 101850252210L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 101850252210L) - ((int128)tmp_q[1] * 709348580227L) + ((int128)tmp_q[2] * 3736891819937L) - ((int128)tmp_q[3] * 3205828084867L) + ((-((int128)tmp_q[4] * 741143456125L) + ((int128)tmp_q[5] * 2330359287744L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 2330359287744L) + ((int128)tmp_q[1] * 101850252210L) - ((int128)tmp_q[2] * 709348580227L) + ((int128)tmp_q[3] * 3736891819937L) - ((int128)tmp_q[4] * 3205828084867L) - ((int128)tmp_q[5] * 4446860736750L);
	tmp_zero[5] = -((int128)tmp_q[0] * 741143456125L) + ((int128)tmp_q[1] * 2330359287744L) + ((int128)tmp_q[2] * 101850252210L) - ((int128)tmp_q[3] * 709348580227L) + ((int128)tmp_q[4] * 3736891819937L) - ((int128)tmp_q[5] * 3205828084867L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

