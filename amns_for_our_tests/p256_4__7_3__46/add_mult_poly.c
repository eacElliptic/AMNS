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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8255209416359351951UL) + ((((uint64_t)op[1] * 8020707734171428624UL) + ((uint64_t)op[2] * 12265687480391165057UL) + ((uint64_t)op[3] * 5080000625035297799UL) + ((uint64_t)op[4] * 16222046545862132401UL) + ((uint64_t)op[5] * 5506633726557625598UL) + ((uint64_t)op[6] * 17783441207941129907UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 17783441207941129907UL) + ((uint64_t)op[1] * 8255209416359351951UL) + ((((uint64_t)op[2] * 8020707734171428624UL) + ((uint64_t)op[3] * 12265687480391165057UL) + ((uint64_t)op[4] * 5080000625035297799UL) + ((uint64_t)op[5] * 16222046545862132401UL) + ((uint64_t)op[6] * 5506633726557625598UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 5506633726557625598UL) + ((uint64_t)op[1] * 17783441207941129907UL) + ((uint64_t)op[2] * 8255209416359351951UL) + ((((uint64_t)op[3] * 8020707734171428624UL) + ((uint64_t)op[4] * 12265687480391165057UL) + ((uint64_t)op[5] * 5080000625035297799UL) + ((uint64_t)op[6] * 16222046545862132401UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 16222046545862132401UL) + ((uint64_t)op[1] * 5506633726557625598UL) + ((uint64_t)op[2] * 17783441207941129907UL) + ((uint64_t)op[3] * 8255209416359351951UL) + ((((uint64_t)op[4] * 8020707734171428624UL) + ((uint64_t)op[5] * 12265687480391165057UL) + ((uint64_t)op[6] * 5080000625035297799UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 5080000625035297799UL) + ((uint64_t)op[1] * 16222046545862132401UL) + ((uint64_t)op[2] * 5506633726557625598UL) + ((uint64_t)op[3] * 17783441207941129907UL) + ((uint64_t)op[4] * 8255209416359351951UL) + ((((uint64_t)op[5] * 8020707734171428624UL) + ((uint64_t)op[6] * 12265687480391165057UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 12265687480391165057UL) + ((uint64_t)op[1] * 5080000625035297799UL) + ((uint64_t)op[2] * 16222046545862132401UL) + ((uint64_t)op[3] * 5506633726557625598UL) + ((uint64_t)op[4] * 17783441207941129907UL) + ((uint64_t)op[5] * 8255209416359351951UL) + ((uint64_t)op[6] * 5615379128804734256UL);
	tmp_q[6] = ((uint64_t)op[0] * 8020707734171428624UL) + ((uint64_t)op[1] * 12265687480391165057UL) + ((uint64_t)op[2] * 5080000625035297799UL) + ((uint64_t)op[3] * 16222046545862132401UL) + ((uint64_t)op[4] * 5506633726557625598UL) + ((uint64_t)op[5] * 17783441207941129907UL) + ((uint64_t)op[6] * 8255209416359351951UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 21138205472L) + ((-((int128)tmp_q[1] * 7835447332L) - ((int128)tmp_q[2] * 42483751152L) - ((int128)tmp_q[3] * 39738421145L) + ((int128)tmp_q[4] * 5694180899L) + ((int128)tmp_q[5] * 49553861205L) + ((int128)tmp_q[6] * 48993739060L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 48993739060L) - ((int128)tmp_q[1] * 21138205472L) + ((-((int128)tmp_q[2] * 7835447332L) - ((int128)tmp_q[3] * 42483751152L) - ((int128)tmp_q[4] * 39738421145L) + ((int128)tmp_q[5] * 5694180899L) + ((int128)tmp_q[6] * 49553861205L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 49553861205L) + ((int128)tmp_q[1] * 48993739060L) - ((int128)tmp_q[2] * 21138205472L) + ((-((int128)tmp_q[3] * 7835447332L) - ((int128)tmp_q[4] * 42483751152L) - ((int128)tmp_q[5] * 39738421145L) + ((int128)tmp_q[6] * 5694180899L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 5694180899L) + ((int128)tmp_q[1] * 49553861205L) + ((int128)tmp_q[2] * 48993739060L) - ((int128)tmp_q[3] * 21138205472L) + ((-((int128)tmp_q[4] * 7835447332L) - ((int128)tmp_q[5] * 42483751152L) - ((int128)tmp_q[6] * 39738421145L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 39738421145L) + ((int128)tmp_q[1] * 5694180899L) + ((int128)tmp_q[2] * 49553861205L) + ((int128)tmp_q[3] * 48993739060L) - ((int128)tmp_q[4] * 21138205472L) + ((-((int128)tmp_q[5] * 7835447332L) - ((int128)tmp_q[6] * 42483751152L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 42483751152L) - ((int128)tmp_q[1] * 39738421145L) + ((int128)tmp_q[2] * 5694180899L) + ((int128)tmp_q[3] * 49553861205L) + ((int128)tmp_q[4] * 48993739060L) - ((int128)tmp_q[5] * 21138205472L) - ((int128)tmp_q[6] * 23506341996L);
	tmp_zero[6] = -((int128)tmp_q[0] * 7835447332L) - ((int128)tmp_q[1] * 42483751152L) - ((int128)tmp_q[2] * 39738421145L) + ((int128)tmp_q[3] * 5694180899L) + ((int128)tmp_q[4] * 49553861205L) + ((int128)tmp_q[5] * 48993739060L) - ((int128)tmp_q[6] * 21138205472L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

