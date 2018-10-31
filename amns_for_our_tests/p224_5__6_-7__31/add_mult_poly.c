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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3141692149071644925UL) + ((((uint64_t)op[1] * 9354326408620187970UL) + ((uint64_t)op[2] * 12296337649898852326UL) + ((uint64_t)op[3] * 12940046079435614499UL) + ((uint64_t)op[4] * 17491480678112287984UL) + ((uint64_t)op[5] * 17163200724027140505UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 17163200724027140505UL) + ((uint64_t)op[1] * 3141692149071644925UL) + ((((uint64_t)op[2] * 9354326408620187970UL) + ((uint64_t)op[3] * 12296337649898852326UL) + ((uint64_t)op[4] * 12940046079435614499UL) + ((uint64_t)op[5] * 17491480678112287984UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 17491480678112287984UL) + ((uint64_t)op[1] * 17163200724027140505UL) + ((uint64_t)op[2] * 3141692149071644925UL) + ((((uint64_t)op[3] * 9354326408620187970UL) + ((uint64_t)op[4] * 12296337649898852326UL) + ((uint64_t)op[5] * 12940046079435614499UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 12940046079435614499UL) + ((uint64_t)op[1] * 17491480678112287984UL) + ((uint64_t)op[2] * 17163200724027140505UL) + ((uint64_t)op[3] * 3141692149071644925UL) + ((((uint64_t)op[4] * 9354326408620187970UL) + ((uint64_t)op[5] * 12296337649898852326UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 12296337649898852326UL) + ((uint64_t)op[1] * 12940046079435614499UL) + ((uint64_t)op[2] * 17491480678112287984UL) + ((uint64_t)op[3] * 17163200724027140505UL) + ((uint64_t)op[4] * 3141692149071644925UL) + ((uint64_t)op[5] * 8306691434496890674UL);
	tmp_q[5] = ((uint64_t)op[0] * 9354326408620187970UL) + ((uint64_t)op[1] * 12296337649898852326UL) + ((uint64_t)op[2] * 12940046079435614499UL) + ((uint64_t)op[3] * 17491480678112287984UL) + ((uint64_t)op[4] * 17163200724027140505UL) + ((uint64_t)op[5] * 3141692149071644925UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 748819158L) - ((((int128)tmp_q[1] * 55938528067L) + ((int128)tmp_q[2] * 29275008497L) + ((int128)tmp_q[3] * 27659701056L) - ((int128)tmp_q[4] * 80838540370L) - ((int128)tmp_q[5] * 65466936469L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 65466936469L) - ((int128)tmp_q[1] * 748819158L) - ((((int128)tmp_q[2] * 55938528067L) + ((int128)tmp_q[3] * 29275008497L) + ((int128)tmp_q[4] * 27659701056L) - ((int128)tmp_q[5] * 80838540370L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 80838540370L) - ((int128)tmp_q[1] * 65466936469L) - ((int128)tmp_q[2] * 748819158L) - ((((int128)tmp_q[3] * 55938528067L) + ((int128)tmp_q[4] * 29275008497L) + ((int128)tmp_q[5] * 27659701056L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 27659701056L) - ((int128)tmp_q[1] * 80838540370L) - ((int128)tmp_q[2] * 65466936469L) - ((int128)tmp_q[3] * 748819158L) - ((((int128)tmp_q[4] * 55938528067L) + ((int128)tmp_q[5] * 29275008497L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 29275008497L) + ((int128)tmp_q[1] * 27659701056L) - ((int128)tmp_q[2] * 80838540370L) - ((int128)tmp_q[3] * 65466936469L) - ((int128)tmp_q[4] * 748819158L) - ((int128)tmp_q[5] * 391569696469L);
	tmp_zero[5] = ((int128)tmp_q[0] * 55938528067L) + ((int128)tmp_q[1] * 29275008497L) + ((int128)tmp_q[2] * 27659701056L) - ((int128)tmp_q[3] * 80838540370L) - ((int128)tmp_q[4] * 65466936469L) - ((int128)tmp_q[5] * 748819158L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

