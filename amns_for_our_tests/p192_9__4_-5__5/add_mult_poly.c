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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8532055964494855364UL) + ((((uint64_t)op[1] * 11384822163364284484UL) + ((uint64_t)op[2] * 15983345233586727290UL) + ((uint64_t)op[3] * 11213112120873248287UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 11213112120873248287UL) + ((uint64_t)op[1] * 8532055964494855364UL) + ((((uint64_t)op[2] * 11384822163364284484UL) + ((uint64_t)op[3] * 15983345233586727290UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 15983345233586727290UL) + ((uint64_t)op[1] * 11213112120873248287UL) + ((uint64_t)op[2] * 8532055964494855364UL) + ((uint64_t)op[3] * 16862865478016784044UL);
	tmp_q[3] = ((uint64_t)op[0] * 11384822163364284484UL) + ((uint64_t)op[1] * 15983345233586727290UL) + ((uint64_t)op[2] * 11213112120873248287UL) + ((uint64_t)op[3] * 8532055964494855364UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 37557818640522L) - ((((int128)tmp_q[1] * 93278953400723L) + ((int128)tmp_q[2] * 98678710670116L) + ((int128)tmp_q[3] * 153458630559512L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 153458630559512L) + ((int128)tmp_q[1] * 37557818640522L) - ((((int128)tmp_q[2] * 93278953400723L) + ((int128)tmp_q[3] * 98678710670116L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 98678710670116L) + ((int128)tmp_q[1] * 153458630559512L) + ((int128)tmp_q[2] * 37557818640522L) - ((int128)tmp_q[3] * 466394767003615L);
	tmp_zero[3] = ((int128)tmp_q[0] * 93278953400723L) + ((int128)tmp_q[1] * 98678710670116L) + ((int128)tmp_q[2] * 153458630559512L) + ((int128)tmp_q[3] * 37557818640522L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

