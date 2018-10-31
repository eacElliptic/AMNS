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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8010740120053338493UL) + ((((uint64_t)op[1] * 10777306374862790436UL) + ((uint64_t)op[2] * 18440018612792091128UL) + ((uint64_t)op[3] * 12809006690139679991UL) + ((uint64_t)op[4] * 16343248817553790019UL) + ((uint64_t)op[5] * 759585215975510714UL) + ((uint64_t)op[6] * 11563348793846729256UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 11563348793846729256UL) + ((uint64_t)op[1] * 8010740120053338493UL) + ((((uint64_t)op[2] * 10777306374862790436UL) + ((uint64_t)op[3] * 18440018612792091128UL) + ((uint64_t)op[4] * 12809006690139679991UL) + ((uint64_t)op[5] * 16343248817553790019UL) + ((uint64_t)op[6] * 759585215975510714UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 759585215975510714UL) + ((uint64_t)op[1] * 11563348793846729256UL) + ((uint64_t)op[2] * 8010740120053338493UL) + ((((uint64_t)op[3] * 10777306374862790436UL) + ((uint64_t)op[4] * 18440018612792091128UL) + ((uint64_t)op[5] * 12809006690139679991UL) + ((uint64_t)op[6] * 16343248817553790019UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 16343248817553790019UL) + ((uint64_t)op[1] * 759585215975510714UL) + ((uint64_t)op[2] * 11563348793846729256UL) + ((uint64_t)op[3] * 8010740120053338493UL) + ((((uint64_t)op[4] * 10777306374862790436UL) + ((uint64_t)op[5] * 18440018612792091128UL) + ((uint64_t)op[6] * 12809006690139679991UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 12809006690139679991UL) + ((uint64_t)op[1] * 16343248817553790019UL) + ((uint64_t)op[2] * 759585215975510714UL) + ((uint64_t)op[3] * 11563348793846729256UL) + ((uint64_t)op[4] * 8010740120053338493UL) + ((((uint64_t)op[5] * 10777306374862790436UL) + ((uint64_t)op[6] * 18440018612792091128UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 18440018612792091128UL) + ((uint64_t)op[1] * 12809006690139679991UL) + ((uint64_t)op[2] * 16343248817553790019UL) + ((uint64_t)op[3] * 759585215975510714UL) + ((uint64_t)op[4] * 11563348793846729256UL) + ((uint64_t)op[5] * 8010740120053338493UL) + ((uint64_t)op[6] * 13885175050878819692UL);
	tmp_q[6] = ((uint64_t)op[0] * 10777306374862790436UL) + ((uint64_t)op[1] * 18440018612792091128UL) + ((uint64_t)op[2] * 12809006690139679991UL) + ((uint64_t)op[3] * 16343248817553790019UL) + ((uint64_t)op[4] * 759585215975510714UL) + ((uint64_t)op[5] * 11563348793846729256UL) + ((uint64_t)op[6] * 8010740120053338493UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6976271545L) + ((((int128)tmp_q[1] * 35932020428L) - ((int128)tmp_q[2] * 20894217805L) + ((int128)tmp_q[3] * 26330390307L) + ((int128)tmp_q[4] * 40520227763L) - ((int128)tmp_q[5] * 56804136625L) + ((int128)tmp_q[6] * 5474560354L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 5474560354L) - ((int128)tmp_q[1] * 6976271545L) + ((((int128)tmp_q[2] * 35932020428L) - ((int128)tmp_q[3] * 20894217805L) + ((int128)tmp_q[4] * 26330390307L) + ((int128)tmp_q[5] * 40520227763L) - ((int128)tmp_q[6] * 56804136625L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 56804136625L) + ((int128)tmp_q[1] * 5474560354L) - ((int128)tmp_q[2] * 6976271545L) + ((((int128)tmp_q[3] * 35932020428L) - ((int128)tmp_q[4] * 20894217805L) + ((int128)tmp_q[5] * 26330390307L) + ((int128)tmp_q[6] * 40520227763L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 40520227763L) - ((int128)tmp_q[1] * 56804136625L) + ((int128)tmp_q[2] * 5474560354L) - ((int128)tmp_q[3] * 6976271545L) + ((((int128)tmp_q[4] * 35932020428L) - ((int128)tmp_q[5] * 20894217805L) + ((int128)tmp_q[6] * 26330390307L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 26330390307L) + ((int128)tmp_q[1] * 40520227763L) - ((int128)tmp_q[2] * 56804136625L) + ((int128)tmp_q[3] * 5474560354L) - ((int128)tmp_q[4] * 6976271545L) + ((((int128)tmp_q[5] * 35932020428L) - ((int128)tmp_q[6] * 20894217805L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 20894217805L) + ((int128)tmp_q[1] * 26330390307L) + ((int128)tmp_q[2] * 40520227763L) - ((int128)tmp_q[3] * 56804136625L) + ((int128)tmp_q[4] * 5474560354L) - ((int128)tmp_q[5] * 6976271545L) + ((int128)tmp_q[6] * 107796061284L);
	tmp_zero[6] = ((int128)tmp_q[0] * 35932020428L) - ((int128)tmp_q[1] * 20894217805L) + ((int128)tmp_q[2] * 26330390307L) + ((int128)tmp_q[3] * 40520227763L) - ((int128)tmp_q[4] * 56804136625L) + ((int128)tmp_q[5] * 5474560354L) - ((int128)tmp_q[6] * 6976271545L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

