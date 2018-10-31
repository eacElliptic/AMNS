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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6686285337222573649UL) + ((((uint64_t)op[1] * 13105431857223574936UL) + ((uint64_t)op[2] * 2892743781305468912UL) + ((uint64_t)op[3] * 4251217178017135812UL) + ((uint64_t)op[4] * 3376674554405468996UL) + ((uint64_t)op[5] * 1061228279311039498UL) + ((uint64_t)op[6] * 13119670598373045434UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 13119670598373045434UL) + ((uint64_t)op[1] * 6686285337222573649UL) + ((((uint64_t)op[2] * 13105431857223574936UL) + ((uint64_t)op[3] * 2892743781305468912UL) + ((uint64_t)op[4] * 4251217178017135812UL) + ((uint64_t)op[5] * 3376674554405468996UL) + ((uint64_t)op[6] * 1061228279311039498UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1061228279311039498UL) + ((uint64_t)op[1] * 13119670598373045434UL) + ((uint64_t)op[2] * 6686285337222573649UL) + ((((uint64_t)op[3] * 13105431857223574936UL) + ((uint64_t)op[4] * 2892743781305468912UL) + ((uint64_t)op[5] * 4251217178017135812UL) + ((uint64_t)op[6] * 3376674554405468996UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 3376674554405468996UL) + ((uint64_t)op[1] * 1061228279311039498UL) + ((uint64_t)op[2] * 13119670598373045434UL) + ((uint64_t)op[3] * 6686285337222573649UL) + ((((uint64_t)op[4] * 13105431857223574936UL) + ((uint64_t)op[5] * 2892743781305468912UL) + ((uint64_t)op[6] * 4251217178017135812UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 4251217178017135812UL) + ((uint64_t)op[1] * 3376674554405468996UL) + ((uint64_t)op[2] * 1061228279311039498UL) + ((uint64_t)op[3] * 13119670598373045434UL) + ((uint64_t)op[4] * 6686285337222573649UL) + ((((uint64_t)op[5] * 13105431857223574936UL) + ((uint64_t)op[6] * 2892743781305468912UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 2892743781305468912UL) + ((uint64_t)op[1] * 4251217178017135812UL) + ((uint64_t)op[2] * 3376674554405468996UL) + ((uint64_t)op[3] * 1061228279311039498UL) + ((uint64_t)op[4] * 13119670598373045434UL) + ((uint64_t)op[5] * 6686285337222573649UL) + ((uint64_t)op[6] * 8259817008720331784UL);
	tmp_q[6] = ((uint64_t)op[0] * 13105431857223574936UL) + ((uint64_t)op[1] * 2892743781305468912UL) + ((uint64_t)op[2] * 4251217178017135812UL) + ((uint64_t)op[3] * 3376674554405468996UL) + ((uint64_t)op[4] * 1061228279311039498UL) + ((uint64_t)op[5] * 13119670598373045434UL) + ((uint64_t)op[6] * 6686285337222573649UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 33522299937L) - ((-((int128)tmp_q[1] * 1095256720L) + ((int128)tmp_q[2] * 66352186072L) + ((int128)tmp_q[3] * 1140414200L) - ((int128)tmp_q[4] * 62615976700L) - ((int128)tmp_q[5] * 33353575642L) + ((int128)tmp_q[6] * 26604521098L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 26604521098L) - ((int128)tmp_q[1] * 33522299937L) - ((-((int128)tmp_q[2] * 1095256720L) + ((int128)tmp_q[3] * 66352186072L) + ((int128)tmp_q[4] * 1140414200L) - ((int128)tmp_q[5] * 62615976700L) - ((int128)tmp_q[6] * 33353575642L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 33353575642L) + ((int128)tmp_q[1] * 26604521098L) - ((int128)tmp_q[2] * 33522299937L) - ((-((int128)tmp_q[3] * 1095256720L) + ((int128)tmp_q[4] * 66352186072L) + ((int128)tmp_q[5] * 1140414200L) - ((int128)tmp_q[6] * 62615976700L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 62615976700L) - ((int128)tmp_q[1] * 33353575642L) + ((int128)tmp_q[2] * 26604521098L) - ((int128)tmp_q[3] * 33522299937L) - ((-((int128)tmp_q[4] * 1095256720L) + ((int128)tmp_q[5] * 66352186072L) + ((int128)tmp_q[6] * 1140414200L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1140414200L) - ((int128)tmp_q[1] * 62615976700L) - ((int128)tmp_q[2] * 33353575642L) + ((int128)tmp_q[3] * 26604521098L) - ((int128)tmp_q[4] * 33522299937L) - ((-((int128)tmp_q[5] * 1095256720L) + ((int128)tmp_q[6] * 66352186072L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 66352186072L) + ((int128)tmp_q[1] * 1140414200L) - ((int128)tmp_q[2] * 62615976700L) - ((int128)tmp_q[3] * 33353575642L) + ((int128)tmp_q[4] * 26604521098L) - ((int128)tmp_q[5] * 33522299937L) + ((int128)tmp_q[6] * 5476283600L);
	tmp_zero[6] = -((int128)tmp_q[0] * 1095256720L) + ((int128)tmp_q[1] * 66352186072L) + ((int128)tmp_q[2] * 1140414200L) - ((int128)tmp_q[3] * 62615976700L) - ((int128)tmp_q[4] * 33353575642L) + ((int128)tmp_q[5] * 26604521098L) - ((int128)tmp_q[6] * 33522299937L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

