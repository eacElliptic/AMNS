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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13531828907252823013UL) + ((((uint64_t)op[1] * 10931078004031736854UL) + ((uint64_t)op[2] * 5181377733757576196UL) + ((uint64_t)op[3] * 6617354547718721513UL) + ((uint64_t)op[4] * 11118918869899696302UL) + ((uint64_t)op[5] * 4440114343429287723UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 4440114343429287723UL) + ((uint64_t)op[1] * 13531828907252823013UL) + ((((uint64_t)op[2] * 10931078004031736854UL) + ((uint64_t)op[3] * 5181377733757576196UL) + ((uint64_t)op[4] * 6617354547718721513UL) + ((uint64_t)op[5] * 11118918869899696302UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 11118918869899696302UL) + ((uint64_t)op[1] * 4440114343429287723UL) + ((uint64_t)op[2] * 13531828907252823013UL) + ((((uint64_t)op[3] * 10931078004031736854UL) + ((uint64_t)op[4] * 5181377733757576196UL) + ((uint64_t)op[5] * 6617354547718721513UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 6617354547718721513UL) + ((uint64_t)op[1] * 11118918869899696302UL) + ((uint64_t)op[2] * 4440114343429287723UL) + ((uint64_t)op[3] * 13531828907252823013UL) + ((((uint64_t)op[4] * 10931078004031736854UL) + ((uint64_t)op[5] * 5181377733757576196UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 5181377733757576196UL) + ((uint64_t)op[1] * 6617354547718721513UL) + ((uint64_t)op[2] * 11118918869899696302UL) + ((uint64_t)op[3] * 4440114343429287723UL) + ((uint64_t)op[4] * 13531828907252823013UL) + ((uint64_t)op[5] * 8200508270647785340UL);
	tmp_q[5] = ((uint64_t)op[0] * 10931078004031736854UL) + ((uint64_t)op[1] * 5181377733757576196UL) + ((uint64_t)op[2] * 6617354547718721513UL) + ((uint64_t)op[3] * 11118918869899696302UL) + ((uint64_t)op[4] * 4440114343429287723UL) + ((uint64_t)op[5] * 13531828907252823013UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 31628681709L) - ((((int128)tmp_q[1] * 105388972606L) - ((int128)tmp_q[2] * 11652683377L) - ((int128)tmp_q[3] * 28020689950L) - ((int128)tmp_q[4] * 39128029781L) - ((int128)tmp_q[5] * 81456812951L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 81456812951L) - ((int128)tmp_q[1] * 31628681709L) - ((((int128)tmp_q[2] * 105388972606L) - ((int128)tmp_q[3] * 11652683377L) - ((int128)tmp_q[4] * 28020689950L) - ((int128)tmp_q[5] * 39128029781L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 39128029781L) - ((int128)tmp_q[1] * 81456812951L) - ((int128)tmp_q[2] * 31628681709L) - ((((int128)tmp_q[3] * 105388972606L) - ((int128)tmp_q[4] * 11652683377L) - ((int128)tmp_q[5] * 28020689950L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 28020689950L) - ((int128)tmp_q[1] * 39128029781L) - ((int128)tmp_q[2] * 81456812951L) - ((int128)tmp_q[3] * 31628681709L) - ((((int128)tmp_q[4] * 105388972606L) - ((int128)tmp_q[5] * 11652683377L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 11652683377L) - ((int128)tmp_q[1] * 28020689950L) - ((int128)tmp_q[2] * 39128029781L) - ((int128)tmp_q[3] * 81456812951L) - ((int128)tmp_q[4] * 31628681709L) - ((int128)tmp_q[5] * 632333835636L);
	tmp_zero[5] = ((int128)tmp_q[0] * 105388972606L) - ((int128)tmp_q[1] * 11652683377L) - ((int128)tmp_q[2] * 28020689950L) - ((int128)tmp_q[3] * 39128029781L) - ((int128)tmp_q[4] * 81456812951L) - ((int128)tmp_q[5] * 31628681709L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

