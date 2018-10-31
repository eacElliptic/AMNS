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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6137523157035667329UL) + ((((uint64_t)op[1] * 11578287944714008788UL) + ((uint64_t)op[2] * 5674838932888607574UL) + ((uint64_t)op[3] * 1225391299913295533UL) + ((uint64_t)op[4] * 13870061575875826475UL) + ((uint64_t)op[5] * 849923260061579818UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 849923260061579818UL) + ((uint64_t)op[1] * 6137523157035667329UL) + ((((uint64_t)op[2] * 11578287944714008788UL) + ((uint64_t)op[3] * 5674838932888607574UL) + ((uint64_t)op[4] * 1225391299913295533UL) + ((uint64_t)op[5] * 13870061575875826475UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 13870061575875826475UL) + ((uint64_t)op[1] * 849923260061579818UL) + ((uint64_t)op[2] * 6137523157035667329UL) + ((((uint64_t)op[3] * 11578287944714008788UL) + ((uint64_t)op[4] * 5674838932888607574UL) + ((uint64_t)op[5] * 1225391299913295533UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 1225391299913295533UL) + ((uint64_t)op[1] * 13870061575875826475UL) + ((uint64_t)op[2] * 849923260061579818UL) + ((uint64_t)op[3] * 6137523157035667329UL) + ((((uint64_t)op[4] * 11578287944714008788UL) + ((uint64_t)op[5] * 5674838932888607574UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 5674838932888607574UL) + ((uint64_t)op[1] * 1225391299913295533UL) + ((uint64_t)op[2] * 13870061575875826475UL) + ((uint64_t)op[3] * 849923260061579818UL) + ((uint64_t)op[4] * 6137523157035667329UL) + ((uint64_t)op[5] * 16288119760432474748UL);
	tmp_q[5] = ((uint64_t)op[0] * 11578287944714008788UL) + ((uint64_t)op[1] * 5674838932888607574UL) + ((uint64_t)op[2] * 1225391299913295533UL) + ((uint64_t)op[3] * 13870061575875826475UL) + ((uint64_t)op[4] * 849923260061579818UL) + ((uint64_t)op[5] * 6137523157035667329UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1324118762L) + ((((int128)tmp_q[1] * 88158803L) - ((int128)tmp_q[2] * 1513056839L) + ((int128)tmp_q[3] * 1250138504L) + ((int128)tmp_q[4] * 1201231435L) - ((int128)tmp_q[5] * 926028514L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 926028514L) - ((int128)tmp_q[1] * 1324118762L) + ((((int128)tmp_q[2] * 88158803L) - ((int128)tmp_q[3] * 1513056839L) + ((int128)tmp_q[4] * 1250138504L) + ((int128)tmp_q[5] * 1201231435L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 1201231435L) - ((int128)tmp_q[1] * 926028514L) - ((int128)tmp_q[2] * 1324118762L) + ((((int128)tmp_q[3] * 88158803L) - ((int128)tmp_q[4] * 1513056839L) + ((int128)tmp_q[5] * 1250138504L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1250138504L) + ((int128)tmp_q[1] * 1201231435L) - ((int128)tmp_q[2] * 926028514L) - ((int128)tmp_q[3] * 1324118762L) + ((((int128)tmp_q[4] * 88158803L) - ((int128)tmp_q[5] * 1513056839L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 1513056839L) + ((int128)tmp_q[1] * 1250138504L) + ((int128)tmp_q[2] * 1201231435L) - ((int128)tmp_q[3] * 926028514L) - ((int128)tmp_q[4] * 1324118762L) + ((int128)tmp_q[5] * 264476409L);
	tmp_zero[5] = ((int128)tmp_q[0] * 88158803L) - ((int128)tmp_q[1] * 1513056839L) + ((int128)tmp_q[2] * 1250138504L) + ((int128)tmp_q[3] * 1201231435L) - ((int128)tmp_q[4] * 926028514L) - ((int128)tmp_q[5] * 1324118762L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

