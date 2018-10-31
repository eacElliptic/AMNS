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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11689827302559183030UL) + ((((uint64_t)op[1] * 438603557140406501UL) + ((uint64_t)op[2] * 1853429344320218400UL) + ((uint64_t)op[3] * 14583443681624311143UL) + ((uint64_t)op[4] * 10532152271113566967UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 10532152271113566967UL) + ((uint64_t)op[1] * 11689827302559183030UL) + ((((uint64_t)op[2] * 438603557140406501UL) + ((uint64_t)op[3] * 1853429344320218400UL) + ((uint64_t)op[4] * 14583443681624311143UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 14583443681624311143UL) + ((uint64_t)op[1] * 10532152271113566967UL) + ((uint64_t)op[2] * 11689827302559183030UL) + ((((uint64_t)op[3] * 438603557140406501UL) + ((uint64_t)op[4] * 1853429344320218400UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 1853429344320218400UL) + ((uint64_t)op[1] * 14583443681624311143UL) + ((uint64_t)op[2] * 10532152271113566967UL) + ((uint64_t)op[3] * 11689827302559183030UL) + ((uint64_t)op[4] * 3070224899982845507UL);
	tmp_q[4] = ((uint64_t)op[0] * 438603557140406501UL) + ((uint64_t)op[1] * 1853429344320218400UL) + ((uint64_t)op[2] * 14583443681624311143UL) + ((uint64_t)op[3] * 10532152271113566967UL) + ((uint64_t)op[4] * 11689827302559183030UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 17654082797019L) + ((((int128)tmp_q[1] * 3941002212374L) + ((int128)tmp_q[2] * 13343671108774L) - ((int128)tmp_q[3] * 8832533083041L) + ((int128)tmp_q[4] * 15828171055079L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 15828171055079L) - ((int128)tmp_q[1] * 17654082797019L) + ((((int128)tmp_q[2] * 3941002212374L) + ((int128)tmp_q[3] * 13343671108774L) - ((int128)tmp_q[4] * 8832533083041L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 8832533083041L) + ((int128)tmp_q[1] * 15828171055079L) - ((int128)tmp_q[2] * 17654082797019L) + ((((int128)tmp_q[3] * 3941002212374L) + ((int128)tmp_q[4] * 13343671108774L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 13343671108774L) - ((int128)tmp_q[1] * 8832533083041L) + ((int128)tmp_q[2] * 15828171055079L) - ((int128)tmp_q[3] * 17654082797019L) + ((int128)tmp_q[4] * 27587015486618L);
	tmp_zero[4] = ((int128)tmp_q[0] * 3941002212374L) + ((int128)tmp_q[1] * 13343671108774L) - ((int128)tmp_q[2] * 8832533083041L) + ((int128)tmp_q[3] * 15828171055079L) - ((int128)tmp_q[4] * 17654082797019L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

