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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8833850911549565499UL) + ((((uint64_t)op[1] * 17083743434838883046UL) + ((uint64_t)op[2] * 17309937286263666853UL) + ((uint64_t)op[3] * 11017668884960018154UL) + ((uint64_t)op[4] * 217455878140450876UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 217455878140450876UL) + ((uint64_t)op[1] * 8833850911549565499UL) + ((((uint64_t)op[2] * 17083743434838883046UL) + ((uint64_t)op[3] * 17309937286263666853UL) + ((uint64_t)op[4] * 11017668884960018154UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 11017668884960018154UL) + ((uint64_t)op[1] * 217455878140450876UL) + ((uint64_t)op[2] * 8833850911549565499UL) + ((((uint64_t)op[3] * 17083743434838883046UL) + ((uint64_t)op[4] * 17309937286263666853UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 17309937286263666853UL) + ((uint64_t)op[1] * 11017668884960018154UL) + ((uint64_t)op[2] * 217455878140450876UL) + ((uint64_t)op[3] * 8833850911549565499UL) + ((uint64_t)op[4] * 15720742795968214476UL);
	tmp_q[4] = ((uint64_t)op[0] * 17083743434838883046UL) + ((uint64_t)op[1] * 17309937286263666853UL) + ((uint64_t)op[2] * 11017668884960018154UL) + ((uint64_t)op[3] * 217455878140450876UL) + ((uint64_t)op[4] * 8833850911549565499UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 10019927770469L) + ((((int128)tmp_q[1] * 18709892912644L) - ((int128)tmp_q[2] * 18375212348535L) + ((int128)tmp_q[3] * 12288753780662L) - ((int128)tmp_q[4] * 5370685210106L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 5370685210106L) + ((int128)tmp_q[1] * 10019927770469L) + ((((int128)tmp_q[2] * 18709892912644L) - ((int128)tmp_q[3] * 18375212348535L) + ((int128)tmp_q[4] * 12288753780662L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 12288753780662L) - ((int128)tmp_q[1] * 5370685210106L) + ((int128)tmp_q[2] * 10019927770469L) + ((((int128)tmp_q[3] * 18709892912644L) - ((int128)tmp_q[4] * 18375212348535L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 18375212348535L) + ((int128)tmp_q[1] * 12288753780662L) - ((int128)tmp_q[2] * 5370685210106L) + ((int128)tmp_q[3] * 10019927770469L) + ((int128)tmp_q[4] * 37419785825288L);
	tmp_zero[4] = ((int128)tmp_q[0] * 18709892912644L) - ((int128)tmp_q[1] * 18375212348535L) + ((int128)tmp_q[2] * 12288753780662L) - ((int128)tmp_q[3] * 5370685210106L) + ((int128)tmp_q[4] * 10019927770469L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

