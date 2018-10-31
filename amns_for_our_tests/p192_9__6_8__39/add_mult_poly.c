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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3605095890718428909UL) + ((((uint64_t)op[1] * 17005299622032109353UL) + ((uint64_t)op[2] * 2990760945840981215UL) + ((uint64_t)op[3] * 14031017329540387614UL) + ((uint64_t)op[4] * 17029721951952173001UL) + ((uint64_t)op[5] * 10367248199616527435UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 10367248199616527435UL) + ((uint64_t)op[1] * 3605095890718428909UL) + ((((uint64_t)op[2] * 17005299622032109353UL) + ((uint64_t)op[3] * 2990760945840981215UL) + ((uint64_t)op[4] * 14031017329540387614UL) + ((uint64_t)op[5] * 17029721951952173001UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 17029721951952173001UL) + ((uint64_t)op[1] * 10367248199616527435UL) + ((uint64_t)op[2] * 3605095890718428909UL) + ((((uint64_t)op[3] * 17005299622032109353UL) + ((uint64_t)op[4] * 2990760945840981215UL) + ((uint64_t)op[5] * 14031017329540387614UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 14031017329540387614UL) + ((uint64_t)op[1] * 17029721951952173001UL) + ((uint64_t)op[2] * 10367248199616527435UL) + ((uint64_t)op[3] * 3605095890718428909UL) + ((((uint64_t)op[4] * 17005299622032109353UL) + ((uint64_t)op[5] * 2990760945840981215UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 2990760945840981215UL) + ((uint64_t)op[1] * 14031017329540387614UL) + ((uint64_t)op[2] * 17029721951952173001UL) + ((uint64_t)op[3] * 10367248199616527435UL) + ((uint64_t)op[4] * 3605095890718428909UL) + ((uint64_t)op[5] * 6915188460290013512UL);
	tmp_q[5] = ((uint64_t)op[0] * 17005299622032109353UL) + ((uint64_t)op[1] * 2990760945840981215UL) + ((uint64_t)op[2] * 14031017329540387614UL) + ((uint64_t)op[3] * 17029721951952173001UL) + ((uint64_t)op[4] * 10367248199616527435UL) + ((uint64_t)op[5] * 3605095890718428909UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 14795596731707L) + ((((int128)tmp_q[1] * 42501301070053L) + ((int128)tmp_q[2] * 26123927650996L) + ((int128)tmp_q[3] * 68913777109251L) - ((int128)tmp_q[4] * 42653157288908L) - ((int128)tmp_q[5] * 45082213822141L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 45082213822141L) + ((int128)tmp_q[1] * 14795596731707L) + ((((int128)tmp_q[2] * 42501301070053L) + ((int128)tmp_q[3] * 26123927650996L) + ((int128)tmp_q[4] * 68913777109251L) - ((int128)tmp_q[5] * 42653157288908L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 42653157288908L) - ((int128)tmp_q[1] * 45082213822141L) + ((int128)tmp_q[2] * 14795596731707L) + ((((int128)tmp_q[3] * 42501301070053L) + ((int128)tmp_q[4] * 26123927650996L) + ((int128)tmp_q[5] * 68913777109251L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 68913777109251L) - ((int128)tmp_q[1] * 42653157288908L) - ((int128)tmp_q[2] * 45082213822141L) + ((int128)tmp_q[3] * 14795596731707L) + ((((int128)tmp_q[4] * 42501301070053L) + ((int128)tmp_q[5] * 26123927650996L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 26123927650996L) + ((int128)tmp_q[1] * 68913777109251L) - ((int128)tmp_q[2] * 42653157288908L) - ((int128)tmp_q[3] * 45082213822141L) + ((int128)tmp_q[4] * 14795596731707L) + ((int128)tmp_q[5] * 340010408560424L);
	tmp_zero[5] = ((int128)tmp_q[0] * 42501301070053L) + ((int128)tmp_q[1] * 26123927650996L) + ((int128)tmp_q[2] * 68913777109251L) - ((int128)tmp_q[3] * 42653157288908L) - ((int128)tmp_q[4] * 45082213822141L) + ((int128)tmp_q[5] * 14795596731707L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

