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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8593608149801795759UL) + ((((uint64_t)op[1] * 15326476944043422688UL) + ((uint64_t)op[2] * 14846800643766699816UL) + ((uint64_t)op[3] * 6440405290091284611UL) + ((uint64_t)op[4] * 5196151286101031891UL) + ((uint64_t)op[5] * 9524218545273481160UL) + ((uint64_t)op[6] * 6587503043414919951UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 6587503043414919951UL) + ((uint64_t)op[1] * 8593608149801795759UL) + ((((uint64_t)op[2] * 15326476944043422688UL) + ((uint64_t)op[3] * 14846800643766699816UL) + ((uint64_t)op[4] * 6440405290091284611UL) + ((uint64_t)op[5] * 5196151286101031891UL) + ((uint64_t)op[6] * 9524218545273481160UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 9524218545273481160UL) + ((uint64_t)op[1] * 6587503043414919951UL) + ((uint64_t)op[2] * 8593608149801795759UL) + ((((uint64_t)op[3] * 15326476944043422688UL) + ((uint64_t)op[4] * 14846800643766699816UL) + ((uint64_t)op[5] * 6440405290091284611UL) + ((uint64_t)op[6] * 5196151286101031891UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 5196151286101031891UL) + ((uint64_t)op[1] * 9524218545273481160UL) + ((uint64_t)op[2] * 6587503043414919951UL) + ((uint64_t)op[3] * 8593608149801795759UL) + ((((uint64_t)op[4] * 15326476944043422688UL) + ((uint64_t)op[5] * 14846800643766699816UL) + ((uint64_t)op[6] * 6440405290091284611UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 6440405290091284611UL) + ((uint64_t)op[1] * 5196151286101031891UL) + ((uint64_t)op[2] * 9524218545273481160UL) + ((uint64_t)op[3] * 6587503043414919951UL) + ((uint64_t)op[4] * 8593608149801795759UL) + ((((uint64_t)op[5] * 15326476944043422688UL) + ((uint64_t)op[6] * 14846800643766699816UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 14846800643766699816UL) + ((uint64_t)op[1] * 6440405290091284611UL) + ((uint64_t)op[2] * 5196151286101031891UL) + ((uint64_t)op[3] * 9524218545273481160UL) + ((uint64_t)op[4] * 6587503043414919951UL) + ((uint64_t)op[5] * 8593608149801795759UL) + ((uint64_t)op[6] * 12206209814377293760UL);
	tmp_q[6] = ((uint64_t)op[0] * 15326476944043422688UL) + ((uint64_t)op[1] * 14846800643766699816UL) + ((uint64_t)op[2] * 6440405290091284611UL) + ((uint64_t)op[3] * 5196151286101031891UL) + ((uint64_t)op[4] * 9524218545273481160UL) + ((uint64_t)op[5] * 6587503043414919951UL) + ((uint64_t)op[6] * 8593608149801795759UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 26899500953L) + ((-((int128)tmp_q[1] * 4788000319L) - ((int128)tmp_q[2] * 57423085104L) + ((int128)tmp_q[3] * 37621224866L) - ((int128)tmp_q[4] * 8900227782L) - ((int128)tmp_q[5] * 53150852215L) + ((int128)tmp_q[6] * 15329040789L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 15329040789L) - ((int128)tmp_q[1] * 26899500953L) + ((-((int128)tmp_q[2] * 4788000319L) - ((int128)tmp_q[3] * 57423085104L) + ((int128)tmp_q[4] * 37621224866L) - ((int128)tmp_q[5] * 8900227782L) - ((int128)tmp_q[6] * 53150852215L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 53150852215L) + ((int128)tmp_q[1] * 15329040789L) - ((int128)tmp_q[2] * 26899500953L) + ((-((int128)tmp_q[3] * 4788000319L) - ((int128)tmp_q[4] * 57423085104L) + ((int128)tmp_q[5] * 37621224866L) - ((int128)tmp_q[6] * 8900227782L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 8900227782L) - ((int128)tmp_q[1] * 53150852215L) + ((int128)tmp_q[2] * 15329040789L) - ((int128)tmp_q[3] * 26899500953L) + ((-((int128)tmp_q[4] * 4788000319L) - ((int128)tmp_q[5] * 57423085104L) + ((int128)tmp_q[6] * 37621224866L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 37621224866L) - ((int128)tmp_q[1] * 8900227782L) - ((int128)tmp_q[2] * 53150852215L) + ((int128)tmp_q[3] * 15329040789L) - ((int128)tmp_q[4] * 26899500953L) + ((-((int128)tmp_q[5] * 4788000319L) - ((int128)tmp_q[6] * 57423085104L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 57423085104L) + ((int128)tmp_q[1] * 37621224866L) - ((int128)tmp_q[2] * 8900227782L) - ((int128)tmp_q[3] * 53150852215L) + ((int128)tmp_q[4] * 15329040789L) - ((int128)tmp_q[5] * 26899500953L) - ((int128)tmp_q[6] * 9576000638L);
	tmp_zero[6] = -((int128)tmp_q[0] * 4788000319L) - ((int128)tmp_q[1] * 57423085104L) + ((int128)tmp_q[2] * 37621224866L) - ((int128)tmp_q[3] * 8900227782L) - ((int128)tmp_q[4] * 53150852215L) + ((int128)tmp_q[5] * 15329040789L) - ((int128)tmp_q[6] * 26899500953L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

