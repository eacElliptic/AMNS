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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[3] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14556129444783801505UL) + ((((uint64_t)op[1] * 17874413178479774009UL) + ((uint64_t)op[2] * 9231673381939475530UL) + ((uint64_t)op[3] * 1930717150984222010UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 1930717150984222010UL) + ((uint64_t)op[1] * 14556129444783801505UL) + ((((uint64_t)op[2] * 17874413178479774009UL) + ((uint64_t)op[3] * 9231673381939475530UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 9231673381939475530UL) + ((uint64_t)op[1] * 1930717150984222010UL) + ((uint64_t)op[2] * 14556129444783801505UL) + ((uint64_t)op[3] * 13868096911871330760UL);
	tmp_q[3] = ((uint64_t)op[0] * 17874413178479774009UL) + ((uint64_t)op[1] * 9231673381939475530UL) + ((uint64_t)op[2] * 1930717150984222010UL) + ((uint64_t)op[3] * 14556129444783801505UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 74855486069663L) + ((-((int128)tmp_q[1] * 92404884657175L) + ((int128)tmp_q[2] * 34058729700510L) + ((int128)tmp_q[3] * 144622032079482L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 144622032079482L) + ((int128)tmp_q[1] * 74855486069663L) + ((-((int128)tmp_q[2] * 92404884657175L) + ((int128)tmp_q[3] * 34058729700510L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 34058729700510L) + ((int128)tmp_q[1] * 144622032079482L) + ((int128)tmp_q[2] * 74855486069663L) - ((int128)tmp_q[3] * 739239077257400L);
	tmp_zero[3] = -((int128)tmp_q[0] * 92404884657175L) + ((int128)tmp_q[1] * 34058729700510L) + ((int128)tmp_q[2] * 144622032079482L) + ((int128)tmp_q[3] * 74855486069663L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

