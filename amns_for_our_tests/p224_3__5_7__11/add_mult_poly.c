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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14158267772612618175UL) + ((((uint64_t)op[1] * 4595174815985665637UL) + ((uint64_t)op[2] * 18172432519245305739UL) + ((uint64_t)op[3] * 13573567646286451258UL) + ((uint64_t)op[4] * 17889402388516551054UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 17889402388516551054UL) + ((uint64_t)op[1] * 14158267772612618175UL) + ((((uint64_t)op[2] * 4595174815985665637UL) + ((uint64_t)op[3] * 18172432519245305739UL) + ((uint64_t)op[4] * 13573567646286451258UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 13573567646286451258UL) + ((uint64_t)op[1] * 17889402388516551054UL) + ((uint64_t)op[2] * 14158267772612618175UL) + ((((uint64_t)op[3] * 4595174815985665637UL) + ((uint64_t)op[4] * 18172432519245305739UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 18172432519245305739UL) + ((uint64_t)op[1] * 13573567646286451258UL) + ((uint64_t)op[2] * 17889402388516551054UL) + ((uint64_t)op[3] * 14158267772612618175UL) + ((uint64_t)op[4] * 13719479638190107843UL);
	tmp_q[4] = ((uint64_t)op[0] * 4595174815985665637UL) + ((uint64_t)op[1] * 18172432519245305739UL) + ((uint64_t)op[2] * 13573567646286451258UL) + ((uint64_t)op[3] * 17889402388516551054UL) + ((uint64_t)op[4] * 14158267772612618175UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 8073500747408L) + ((-((int128)tmp_q[1] * 250587414357L) - ((int128)tmp_q[2] * 9905877218925L) + ((int128)tmp_q[3] * 677842970986L) - ((int128)tmp_q[4] * 16067130217017L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 16067130217017L) + ((int128)tmp_q[1] * 8073500747408L) + ((-((int128)tmp_q[2] * 250587414357L) - ((int128)tmp_q[3] * 9905877218925L) + ((int128)tmp_q[4] * 677842970986L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 677842970986L) - ((int128)tmp_q[1] * 16067130217017L) + ((int128)tmp_q[2] * 8073500747408L) + ((-((int128)tmp_q[3] * 250587414357L) - ((int128)tmp_q[4] * 9905877218925L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 9905877218925L) + ((int128)tmp_q[1] * 677842970986L) - ((int128)tmp_q[2] * 16067130217017L) + ((int128)tmp_q[3] * 8073500747408L) - ((int128)tmp_q[4] * 1754111900499L);
	tmp_zero[4] = -((int128)tmp_q[0] * 250587414357L) - ((int128)tmp_q[1] * 9905877218925L) + ((int128)tmp_q[2] * 677842970986L) - ((int128)tmp_q[3] * 16067130217017L) + ((int128)tmp_q[4] * 8073500747408L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

