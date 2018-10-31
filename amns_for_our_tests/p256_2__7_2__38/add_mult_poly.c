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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11473166836076692457UL) + ((((uint64_t)op[1] * 16338735519545254193UL) + ((uint64_t)op[2] * 6468757196959249622UL) + ((uint64_t)op[3] * 343111811518433375UL) + ((uint64_t)op[4] * 15242717238033716950UL) + ((uint64_t)op[5] * 2200178495038581825UL) + ((uint64_t)op[6] * 7045854275977210730UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 7045854275977210730UL) + ((uint64_t)op[1] * 11473166836076692457UL) + ((((uint64_t)op[2] * 16338735519545254193UL) + ((uint64_t)op[3] * 6468757196959249622UL) + ((uint64_t)op[4] * 343111811518433375UL) + ((uint64_t)op[5] * 15242717238033716950UL) + ((uint64_t)op[6] * 2200178495038581825UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 2200178495038581825UL) + ((uint64_t)op[1] * 7045854275977210730UL) + ((uint64_t)op[2] * 11473166836076692457UL) + ((((uint64_t)op[3] * 16338735519545254193UL) + ((uint64_t)op[4] * 6468757196959249622UL) + ((uint64_t)op[5] * 343111811518433375UL) + ((uint64_t)op[6] * 15242717238033716950UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 15242717238033716950UL) + ((uint64_t)op[1] * 2200178495038581825UL) + ((uint64_t)op[2] * 7045854275977210730UL) + ((uint64_t)op[3] * 11473166836076692457UL) + ((((uint64_t)op[4] * 16338735519545254193UL) + ((uint64_t)op[5] * 6468757196959249622UL) + ((uint64_t)op[6] * 343111811518433375UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 343111811518433375UL) + ((uint64_t)op[1] * 15242717238033716950UL) + ((uint64_t)op[2] * 2200178495038581825UL) + ((uint64_t)op[3] * 7045854275977210730UL) + ((uint64_t)op[4] * 11473166836076692457UL) + ((((uint64_t)op[5] * 16338735519545254193UL) + ((uint64_t)op[6] * 6468757196959249622UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 6468757196959249622UL) + ((uint64_t)op[1] * 343111811518433375UL) + ((uint64_t)op[2] * 15242717238033716950UL) + ((uint64_t)op[3] * 2200178495038581825UL) + ((uint64_t)op[4] * 7045854275977210730UL) + ((uint64_t)op[5] * 11473166836076692457UL) + ((uint64_t)op[6] * 14230726965380956770UL);
	tmp_q[6] = ((uint64_t)op[0] * 16338735519545254193UL) + ((uint64_t)op[1] * 6468757196959249622UL) + ((uint64_t)op[2] * 343111811518433375UL) + ((uint64_t)op[3] * 15242717238033716950UL) + ((uint64_t)op[4] * 2200178495038581825UL) + ((uint64_t)op[5] * 7045854275977210730UL) + ((uint64_t)op[6] * 11473166836076692457UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 33824520987L) + ((-((int128)tmp_q[1] * 33913392156L) - ((int128)tmp_q[2] * 18070021172L) + ((int128)tmp_q[3] * 24749858646L) + ((int128)tmp_q[4] * 58197840636L) + ((int128)tmp_q[5] * 28803167405L) - ((int128)tmp_q[6] * 43266996932L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 43266996932L) + ((int128)tmp_q[1] * 33824520987L) + ((-((int128)tmp_q[2] * 33913392156L) - ((int128)tmp_q[3] * 18070021172L) + ((int128)tmp_q[4] * 24749858646L) + ((int128)tmp_q[5] * 58197840636L) + ((int128)tmp_q[6] * 28803167405L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 28803167405L) - ((int128)tmp_q[1] * 43266996932L) + ((int128)tmp_q[2] * 33824520987L) + ((-((int128)tmp_q[3] * 33913392156L) - ((int128)tmp_q[4] * 18070021172L) + ((int128)tmp_q[5] * 24749858646L) + ((int128)tmp_q[6] * 58197840636L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 58197840636L) + ((int128)tmp_q[1] * 28803167405L) - ((int128)tmp_q[2] * 43266996932L) + ((int128)tmp_q[3] * 33824520987L) + ((-((int128)tmp_q[4] * 33913392156L) - ((int128)tmp_q[5] * 18070021172L) + ((int128)tmp_q[6] * 24749858646L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 24749858646L) + ((int128)tmp_q[1] * 58197840636L) + ((int128)tmp_q[2] * 28803167405L) - ((int128)tmp_q[3] * 43266996932L) + ((int128)tmp_q[4] * 33824520987L) + ((-((int128)tmp_q[5] * 33913392156L) - ((int128)tmp_q[6] * 18070021172L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 18070021172L) + ((int128)tmp_q[1] * 24749858646L) + ((int128)tmp_q[2] * 58197840636L) + ((int128)tmp_q[3] * 28803167405L) - ((int128)tmp_q[4] * 43266996932L) + ((int128)tmp_q[5] * 33824520987L) - ((int128)tmp_q[6] * 67826784312L);
	tmp_zero[6] = -((int128)tmp_q[0] * 33913392156L) - ((int128)tmp_q[1] * 18070021172L) + ((int128)tmp_q[2] * 24749858646L) + ((int128)tmp_q[3] * 58197840636L) + ((int128)tmp_q[4] * 28803167405L) - ((int128)tmp_q[5] * 43266996932L) + ((int128)tmp_q[6] * 33824520987L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

