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
	tmp_q[0] = ((uint64_t)op[0] * 12421240620220707599UL) + ((((uint64_t)op[1] * 18039524544705255274UL) + ((uint64_t)op[2] * 13877558409615725329UL) + ((uint64_t)op[3] * 363601807860198527UL) + ((uint64_t)op[4] * 13240685104086641776UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 13240685104086641776UL) + ((uint64_t)op[1] * 12421240620220707599UL) + ((((uint64_t)op[2] * 18039524544705255274UL) + ((uint64_t)op[3] * 13877558409615725329UL) + ((uint64_t)op[4] * 363601807860198527UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 363601807860198527UL) + ((uint64_t)op[1] * 13240685104086641776UL) + ((uint64_t)op[2] * 12421240620220707599UL) + ((((uint64_t)op[3] * 18039524544705255274UL) + ((uint64_t)op[4] * 13877558409615725329UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 13877558409615725329UL) + ((uint64_t)op[1] * 363601807860198527UL) + ((uint64_t)op[2] * 13240685104086641776UL) + ((uint64_t)op[3] * 12421240620220707599UL) + ((uint64_t)op[4] * 17632305015700958932UL);
	tmp_q[4] = ((uint64_t)op[0] * 18039524544705255274UL) + ((uint64_t)op[1] * 13877558409615725329UL) + ((uint64_t)op[2] * 363601807860198527UL) + ((uint64_t)op[3] * 13240685104086641776UL) + ((uint64_t)op[4] * 12421240620220707599UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 74662383497L) + ((-((int128)tmp_q[1] * 16107695551L) + ((int128)tmp_q[2] * 109312170953L) + ((int128)tmp_q[3] * 149929971869L) - ((int128)tmp_q[4] * 266278875804L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 266278875804L) + ((int128)tmp_q[1] * 74662383497L) + ((-((int128)tmp_q[2] * 16107695551L) + ((int128)tmp_q[3] * 109312170953L) + ((int128)tmp_q[4] * 149929971869L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 149929971869L) - ((int128)tmp_q[1] * 266278875804L) + ((int128)tmp_q[2] * 74662383497L) + ((-((int128)tmp_q[3] * 16107695551L) + ((int128)tmp_q[4] * 109312170953L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 109312170953L) + ((int128)tmp_q[1] * 149929971869L) - ((int128)tmp_q[2] * 266278875804L) + ((int128)tmp_q[3] * 74662383497L) - ((int128)tmp_q[4] * 32215391102L);
	tmp_zero[4] = -((int128)tmp_q[0] * 16107695551L) + ((int128)tmp_q[1] * 109312170953L) + ((int128)tmp_q[2] * 149929971869L) - ((int128)tmp_q[3] * 266278875804L) + ((int128)tmp_q[4] * 74662383497L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

