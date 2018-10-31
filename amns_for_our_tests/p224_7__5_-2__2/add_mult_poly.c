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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16181873410774310829UL) + ((((uint64_t)op[1] * 168919791232449045UL) + ((uint64_t)op[2] * 6315499984030923453UL) + ((uint64_t)op[3] * 2813606731364194461UL) + ((uint64_t)op[4] * 10640821038796936179UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 10640821038796936179UL) + ((uint64_t)op[1] * 16181873410774310829UL) + ((((uint64_t)op[2] * 168919791232449045UL) + ((uint64_t)op[3] * 6315499984030923453UL) + ((uint64_t)op[4] * 2813606731364194461UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2813606731364194461UL) + ((uint64_t)op[1] * 10640821038796936179UL) + ((uint64_t)op[2] * 16181873410774310829UL) + ((((uint64_t)op[3] * 168919791232449045UL) + ((uint64_t)op[4] * 6315499984030923453UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 6315499984030923453UL) + ((uint64_t)op[1] * 2813606731364194461UL) + ((uint64_t)op[2] * 10640821038796936179UL) + ((uint64_t)op[3] * 16181873410774310829UL) + ((uint64_t)op[4] * 18108904491244653526UL);
	tmp_q[4] = ((uint64_t)op[0] * 168919791232449045UL) + ((uint64_t)op[1] * 6315499984030923453UL) + ((uint64_t)op[2] * 2813606731364194461UL) + ((uint64_t)op[3] * 10640821038796936179UL) + ((uint64_t)op[4] * 16181873410774310829UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 23696030476079L) - ((((int128)tmp_q[1] * 7614709056276L) + ((int128)tmp_q[2] * 8295603458306L) + ((int128)tmp_q[3] * 5937304702008L) + ((int128)tmp_q[4] * 10334861225357L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 10334861225357L) - ((int128)tmp_q[1] * 23696030476079L) - ((((int128)tmp_q[2] * 7614709056276L) + ((int128)tmp_q[3] * 8295603458306L) + ((int128)tmp_q[4] * 5937304702008L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 5937304702008L) + ((int128)tmp_q[1] * 10334861225357L) - ((int128)tmp_q[2] * 23696030476079L) - ((((int128)tmp_q[3] * 7614709056276L) + ((int128)tmp_q[4] * 8295603458306L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 8295603458306L) + ((int128)tmp_q[1] * 5937304702008L) + ((int128)tmp_q[2] * 10334861225357L) - ((int128)tmp_q[3] * 23696030476079L) - ((int128)tmp_q[4] * 15229418112552L);
	tmp_zero[4] = ((int128)tmp_q[0] * 7614709056276L) + ((int128)tmp_q[1] * 8295603458306L) + ((int128)tmp_q[2] * 5937304702008L) + ((int128)tmp_q[3] * 10334861225357L) - ((int128)tmp_q[4] * 23696030476079L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

