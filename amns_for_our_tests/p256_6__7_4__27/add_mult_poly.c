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
	tmp_q[0] = ((uint64_t)op[0] * 11166999560050885621UL) + ((((uint64_t)op[1] * 11181835378205512552UL) + ((uint64_t)op[2] * 17579151192417079834UL) + ((uint64_t)op[3] * 10317624453454712108UL) + ((uint64_t)op[4] * 1286619325004665618UL) + ((uint64_t)op[5] * 8401789891454390295UL) + ((uint64_t)op[6] * 12048883510852937423UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 12048883510852937423UL) + ((uint64_t)op[1] * 11166999560050885621UL) + ((((uint64_t)op[2] * 11181835378205512552UL) + ((uint64_t)op[3] * 17579151192417079834UL) + ((uint64_t)op[4] * 10317624453454712108UL) + ((uint64_t)op[5] * 1286619325004665618UL) + ((uint64_t)op[6] * 8401789891454390295UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 8401789891454390295UL) + ((uint64_t)op[1] * 12048883510852937423UL) + ((uint64_t)op[2] * 11166999560050885621UL) + ((((uint64_t)op[3] * 11181835378205512552UL) + ((uint64_t)op[4] * 17579151192417079834UL) + ((uint64_t)op[5] * 10317624453454712108UL) + ((uint64_t)op[6] * 1286619325004665618UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 1286619325004665618UL) + ((uint64_t)op[1] * 8401789891454390295UL) + ((uint64_t)op[2] * 12048883510852937423UL) + ((uint64_t)op[3] * 11166999560050885621UL) + ((((uint64_t)op[4] * 11181835378205512552UL) + ((uint64_t)op[5] * 17579151192417079834UL) + ((uint64_t)op[6] * 10317624453454712108UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 10317624453454712108UL) + ((uint64_t)op[1] * 1286619325004665618UL) + ((uint64_t)op[2] * 8401789891454390295UL) + ((uint64_t)op[3] * 12048883510852937423UL) + ((uint64_t)op[4] * 11166999560050885621UL) + ((((uint64_t)op[5] * 11181835378205512552UL) + ((uint64_t)op[6] * 17579151192417079834UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 17579151192417079834UL) + ((uint64_t)op[1] * 10317624453454712108UL) + ((uint64_t)op[2] * 1286619325004665618UL) + ((uint64_t)op[3] * 8401789891454390295UL) + ((uint64_t)op[4] * 12048883510852937423UL) + ((uint64_t)op[5] * 11166999560050885621UL) + ((uint64_t)op[6] * 7833853365402946976UL);
	tmp_q[6] = ((uint64_t)op[0] * 11181835378205512552UL) + ((uint64_t)op[1] * 17579151192417079834UL) + ((uint64_t)op[2] * 10317624453454712108UL) + ((uint64_t)op[3] * 1286619325004665618UL) + ((uint64_t)op[4] * 8401789891454390295UL) + ((uint64_t)op[5] * 12048883510852937423UL) + ((uint64_t)op[6] * 11166999560050885621UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 22990264551L) + ((-((int128)tmp_q[1] * 3558153077L) + ((int128)tmp_q[2] * 1357532504L) - ((int128)tmp_q[3] * 63760122253L) + ((int128)tmp_q[4] * 66175310579L) + ((int128)tmp_q[5] * 4384714774L) + ((int128)tmp_q[6] * 34527633711L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 34527633711L) + ((int128)tmp_q[1] * 22990264551L) + ((-((int128)tmp_q[2] * 3558153077L) + ((int128)tmp_q[3] * 1357532504L) - ((int128)tmp_q[4] * 63760122253L) + ((int128)tmp_q[5] * 66175310579L) + ((int128)tmp_q[6] * 4384714774L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 4384714774L) + ((int128)tmp_q[1] * 34527633711L) + ((int128)tmp_q[2] * 22990264551L) + ((-((int128)tmp_q[3] * 3558153077L) + ((int128)tmp_q[4] * 1357532504L) - ((int128)tmp_q[5] * 63760122253L) + ((int128)tmp_q[6] * 66175310579L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 66175310579L) + ((int128)tmp_q[1] * 4384714774L) + ((int128)tmp_q[2] * 34527633711L) + ((int128)tmp_q[3] * 22990264551L) + ((-((int128)tmp_q[4] * 3558153077L) + ((int128)tmp_q[5] * 1357532504L) - ((int128)tmp_q[6] * 63760122253L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 63760122253L) + ((int128)tmp_q[1] * 66175310579L) + ((int128)tmp_q[2] * 4384714774L) + ((int128)tmp_q[3] * 34527633711L) + ((int128)tmp_q[4] * 22990264551L) + ((-((int128)tmp_q[5] * 3558153077L) + ((int128)tmp_q[6] * 1357532504L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 1357532504L) - ((int128)tmp_q[1] * 63760122253L) + ((int128)tmp_q[2] * 66175310579L) + ((int128)tmp_q[3] * 4384714774L) + ((int128)tmp_q[4] * 34527633711L) + ((int128)tmp_q[5] * 22990264551L) - ((int128)tmp_q[6] * 14232612308L);
	tmp_zero[6] = -((int128)tmp_q[0] * 3558153077L) + ((int128)tmp_q[1] * 1357532504L) - ((int128)tmp_q[2] * 63760122253L) + ((int128)tmp_q[3] * 66175310579L) + ((int128)tmp_q[4] * 4384714774L) + ((int128)tmp_q[5] * 34527633711L) + ((int128)tmp_q[6] * 22990264551L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

