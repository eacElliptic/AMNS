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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4203915704000889295UL) + ((((uint64_t)op[1] * 8573924981176199562UL) + ((uint64_t)op[2] * 11322267856134404221UL) + ((uint64_t)op[3] * 4161970180214934385UL) + ((uint64_t)op[4] * 4619632370718191504UL) + ((uint64_t)op[5] * 1398547629167344129UL) + ((uint64_t)op[6] * 3758513694242011902UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 3758513694242011902UL) + ((uint64_t)op[1] * 4203915704000889295UL) + ((((uint64_t)op[2] * 8573924981176199562UL) + ((uint64_t)op[3] * 11322267856134404221UL) + ((uint64_t)op[4] * 4161970180214934385UL) + ((uint64_t)op[5] * 4619632370718191504UL) + ((uint64_t)op[6] * 1398547629167344129UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 1398547629167344129UL) + ((uint64_t)op[1] * 3758513694242011902UL) + ((uint64_t)op[2] * 4203915704000889295UL) + ((((uint64_t)op[3] * 8573924981176199562UL) + ((uint64_t)op[4] * 11322267856134404221UL) + ((uint64_t)op[5] * 4161970180214934385UL) + ((uint64_t)op[6] * 4619632370718191504UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 4619632370718191504UL) + ((uint64_t)op[1] * 1398547629167344129UL) + ((uint64_t)op[2] * 3758513694242011902UL) + ((uint64_t)op[3] * 4203915704000889295UL) + ((((uint64_t)op[4] * 8573924981176199562UL) + ((uint64_t)op[5] * 11322267856134404221UL) + ((uint64_t)op[6] * 4161970180214934385UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 4161970180214934385UL) + ((uint64_t)op[1] * 4619632370718191504UL) + ((uint64_t)op[2] * 1398547629167344129UL) + ((uint64_t)op[3] * 3758513694242011902UL) + ((uint64_t)op[4] * 4203915704000889295UL) + ((((uint64_t)op[5] * 8573924981176199562UL) + ((uint64_t)op[6] * 11322267856134404221UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 11322267856134404221UL) + ((uint64_t)op[1] * 4161970180214934385UL) + ((uint64_t)op[2] * 4619632370718191504UL) + ((uint64_t)op[3] * 1398547629167344129UL) + ((uint64_t)op[4] * 3758513694242011902UL) + ((uint64_t)op[5] * 4203915704000889295UL) + ((uint64_t)op[6] * 15848955850995246632UL);
	tmp_q[6] = ((uint64_t)op[0] * 8573924981176199562UL) + ((uint64_t)op[1] * 11322267856134404221UL) + ((uint64_t)op[2] * 4161970180214934385UL) + ((uint64_t)op[3] * 4619632370718191504UL) + ((uint64_t)op[4] * 1398547629167344129UL) + ((uint64_t)op[5] * 3758513694242011902UL) + ((uint64_t)op[6] * 4203915704000889295UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 31483720473L) + ((-((int128)tmp_q[1] * 24399545883L) + ((int128)tmp_q[2] * 39197671779L) - ((int128)tmp_q[3] * 6805594930L) + ((int128)tmp_q[4] * 29280925888L) + ((int128)tmp_q[5] * 17138643977L) + ((int128)tmp_q[6] * 75427326690L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 75427326690L) + ((int128)tmp_q[1] * 31483720473L) + ((-((int128)tmp_q[2] * 24399545883L) + ((int128)tmp_q[3] * 39197671779L) - ((int128)tmp_q[4] * 6805594930L) + ((int128)tmp_q[5] * 29280925888L) + ((int128)tmp_q[6] * 17138643977L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 17138643977L) + ((int128)tmp_q[1] * 75427326690L) + ((int128)tmp_q[2] * 31483720473L) + ((-((int128)tmp_q[3] * 24399545883L) + ((int128)tmp_q[4] * 39197671779L) - ((int128)tmp_q[5] * 6805594930L) + ((int128)tmp_q[6] * 29280925888L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 29280925888L) + ((int128)tmp_q[1] * 17138643977L) + ((int128)tmp_q[2] * 75427326690L) + ((int128)tmp_q[3] * 31483720473L) + ((-((int128)tmp_q[4] * 24399545883L) + ((int128)tmp_q[5] * 39197671779L) - ((int128)tmp_q[6] * 6805594930L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 6805594930L) + ((int128)tmp_q[1] * 29280925888L) + ((int128)tmp_q[2] * 17138643977L) + ((int128)tmp_q[3] * 75427326690L) + ((int128)tmp_q[4] * 31483720473L) + ((-((int128)tmp_q[5] * 24399545883L) + ((int128)tmp_q[6] * 39197671779L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 39197671779L) - ((int128)tmp_q[1] * 6805594930L) + ((int128)tmp_q[2] * 29280925888L) + ((int128)tmp_q[3] * 17138643977L) + ((int128)tmp_q[4] * 75427326690L) + ((int128)tmp_q[5] * 31483720473L) - ((int128)tmp_q[6] * 97598183532L);
	tmp_zero[6] = -((int128)tmp_q[0] * 24399545883L) + ((int128)tmp_q[1] * 39197671779L) - ((int128)tmp_q[2] * 6805594930L) + ((int128)tmp_q[3] * 29280925888L) + ((int128)tmp_q[4] * 17138643977L) + ((int128)tmp_q[5] * 75427326690L) + ((int128)tmp_q[6] * 31483720473L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

