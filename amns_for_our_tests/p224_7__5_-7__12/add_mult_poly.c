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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1378242949293080813UL) + ((((uint64_t)op[1] * 11814375741789320001UL) + ((uint64_t)op[2] * 6038599845425625591UL) + ((uint64_t)op[3] * 12301918979697084274UL) + ((uint64_t)op[4] * 631124116883667160UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 631124116883667160UL) + ((uint64_t)op[1] * 1378242949293080813UL) + ((((uint64_t)op[2] * 11814375741789320001UL) + ((uint64_t)op[3] * 6038599845425625591UL) + ((uint64_t)op[4] * 12301918979697084274UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 12301918979697084274UL) + ((uint64_t)op[1] * 631124116883667160UL) + ((uint64_t)op[2] * 1378242949293080813UL) + ((((uint64_t)op[3] * 11814375741789320001UL) + ((uint64_t)op[4] * 6038599845425625591UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 6038599845425625591UL) + ((uint64_t)op[1] * 12301918979697084274UL) + ((uint64_t)op[2] * 631124116883667160UL) + ((uint64_t)op[3] * 1378242949293080813UL) + ((uint64_t)op[4] * 9533090176022518073UL);
	tmp_q[4] = ((uint64_t)op[0] * 11814375741789320001UL) + ((uint64_t)op[1] * 6038599845425625591UL) + ((uint64_t)op[2] * 12301918979697084274UL) + ((uint64_t)op[3] * 631124116883667160UL) + ((uint64_t)op[4] * 1378242949293080813UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5955067730926L) - ((((int128)tmp_q[1] * 1717977415651L) - ((int128)tmp_q[2] * 5519943139675L) + ((int128)tmp_q[3] * 12601189079898L) + ((int128)tmp_q[4] * 11325719064345L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 11325719064345L) + ((int128)tmp_q[1] * 5955067730926L) - ((((int128)tmp_q[2] * 1717977415651L) - ((int128)tmp_q[3] * 5519943139675L) + ((int128)tmp_q[4] * 12601189079898L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 12601189079898L) + ((int128)tmp_q[1] * 11325719064345L) + ((int128)tmp_q[2] * 5955067730926L) - ((((int128)tmp_q[3] * 1717977415651L) - ((int128)tmp_q[4] * 5519943139675L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 5519943139675L) + ((int128)tmp_q[1] * 12601189079898L) + ((int128)tmp_q[2] * 11325719064345L) + ((int128)tmp_q[3] * 5955067730926L) - ((int128)tmp_q[4] * 12025841909557L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1717977415651L) - ((int128)tmp_q[1] * 5519943139675L) + ((int128)tmp_q[2] * 12601189079898L) + ((int128)tmp_q[3] * 11325719064345L) + ((int128)tmp_q[4] * 5955067730926L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

