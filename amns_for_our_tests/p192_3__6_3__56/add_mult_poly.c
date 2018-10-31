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
	tmp_q[0] = ((uint64_t)op[0] * 13420998225721242092UL) + ((((uint64_t)op[1] * 8069167330477353880UL) + ((uint64_t)op[2] * 13042173382400082827UL) + ((uint64_t)op[3] * 13460961250822672142UL) + ((uint64_t)op[4] * 16039177128360791268UL) + ((uint64_t)op[5] * 8591013385688375416UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 8591013385688375416UL) + ((uint64_t)op[1] * 13420998225721242092UL) + ((((uint64_t)op[2] * 8069167330477353880UL) + ((uint64_t)op[3] * 13042173382400082827UL) + ((uint64_t)op[4] * 13460961250822672142UL) + ((uint64_t)op[5] * 16039177128360791268UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 16039177128360791268UL) + ((uint64_t)op[1] * 8591013385688375416UL) + ((uint64_t)op[2] * 13420998225721242092UL) + ((((uint64_t)op[3] * 8069167330477353880UL) + ((uint64_t)op[4] * 13042173382400082827UL) + ((uint64_t)op[5] * 13460961250822672142UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 13460961250822672142UL) + ((uint64_t)op[1] * 16039177128360791268UL) + ((uint64_t)op[2] * 8591013385688375416UL) + ((uint64_t)op[3] * 13420998225721242092UL) + ((((uint64_t)op[4] * 8069167330477353880UL) + ((uint64_t)op[5] * 13042173382400082827UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 13042173382400082827UL) + ((uint64_t)op[1] * 13460961250822672142UL) + ((uint64_t)op[2] * 16039177128360791268UL) + ((uint64_t)op[3] * 8591013385688375416UL) + ((uint64_t)op[4] * 13420998225721242092UL) + ((uint64_t)op[5] * 5760757917722510024UL);
	tmp_q[5] = ((uint64_t)op[0] * 8069167330477353880UL) + ((uint64_t)op[1] * 13042173382400082827UL) + ((uint64_t)op[2] * 13460961250822672142UL) + ((uint64_t)op[3] * 16039177128360791268UL) + ((uint64_t)op[4] * 8591013385688375416UL) + ((uint64_t)op[5] * 13420998225721242092UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1267736L) + ((((int128)tmp_q[1] * 786320704L) + ((int128)tmp_q[2] * 1012683228L) + ((int128)tmp_q[3] * 337964984L) - ((int128)tmp_q[4] * 3223485777L) - ((int128)tmp_q[5] * 2241858230L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 2241858230L) + ((int128)tmp_q[1] * 1267736L) + ((((int128)tmp_q[2] * 786320704L) + ((int128)tmp_q[3] * 1012683228L) + ((int128)tmp_q[4] * 337964984L) - ((int128)tmp_q[5] * 3223485777L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 3223485777L) - ((int128)tmp_q[1] * 2241858230L) + ((int128)tmp_q[2] * 1267736L) + ((((int128)tmp_q[3] * 786320704L) + ((int128)tmp_q[4] * 1012683228L) + ((int128)tmp_q[5] * 337964984L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 337964984L) - ((int128)tmp_q[1] * 3223485777L) - ((int128)tmp_q[2] * 2241858230L) + ((int128)tmp_q[3] * 1267736L) + ((((int128)tmp_q[4] * 786320704L) + ((int128)tmp_q[5] * 1012683228L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1012683228L) + ((int128)tmp_q[1] * 337964984L) - ((int128)tmp_q[2] * 3223485777L) - ((int128)tmp_q[3] * 2241858230L) + ((int128)tmp_q[4] * 1267736L) + ((int128)tmp_q[5] * 2358962112L);
	tmp_zero[5] = ((int128)tmp_q[0] * 786320704L) + ((int128)tmp_q[1] * 1012683228L) + ((int128)tmp_q[2] * 337964984L) - ((int128)tmp_q[3] * 3223485777L) - ((int128)tmp_q[4] * 2241858230L) + ((int128)tmp_q[5] * 1267736L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

