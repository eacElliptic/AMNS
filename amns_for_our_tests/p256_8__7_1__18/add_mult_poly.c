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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + ((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + ((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + ((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + ((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + ((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + ((int128)pa[6] * pb[6]);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + ((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + ((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + ((int128)pa[6] * pa[6]);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17130256066536578734UL) + ((uint64_t)op[1] * 1617787785200538278UL) + ((uint64_t)op[2] * 376288571705508531UL) + ((uint64_t)op[3] * 4882460954289740432UL) + ((uint64_t)op[4] * 259251780415428890UL) + ((uint64_t)op[5] * 15650751168344772518UL) + ((uint64_t)op[6] * 7517688434474851058UL);
	tmp_q[1] = ((uint64_t)op[0] * 7517688434474851058UL) + ((uint64_t)op[1] * 17130256066536578734UL) + ((uint64_t)op[2] * 1617787785200538278UL) + ((uint64_t)op[3] * 376288571705508531UL) + ((uint64_t)op[4] * 4882460954289740432UL) + ((uint64_t)op[5] * 259251780415428890UL) + ((uint64_t)op[6] * 15650751168344772518UL);
	tmp_q[2] = ((uint64_t)op[0] * 15650751168344772518UL) + ((uint64_t)op[1] * 7517688434474851058UL) + ((uint64_t)op[2] * 17130256066536578734UL) + ((uint64_t)op[3] * 1617787785200538278UL) + ((uint64_t)op[4] * 376288571705508531UL) + ((uint64_t)op[5] * 4882460954289740432UL) + ((uint64_t)op[6] * 259251780415428890UL);
	tmp_q[3] = ((uint64_t)op[0] * 259251780415428890UL) + ((uint64_t)op[1] * 15650751168344772518UL) + ((uint64_t)op[2] * 7517688434474851058UL) + ((uint64_t)op[3] * 17130256066536578734UL) + ((uint64_t)op[4] * 1617787785200538278UL) + ((uint64_t)op[5] * 376288571705508531UL) + ((uint64_t)op[6] * 4882460954289740432UL);
	tmp_q[4] = ((uint64_t)op[0] * 4882460954289740432UL) + ((uint64_t)op[1] * 259251780415428890UL) + ((uint64_t)op[2] * 15650751168344772518UL) + ((uint64_t)op[3] * 7517688434474851058UL) + ((uint64_t)op[4] * 17130256066536578734UL) + ((uint64_t)op[5] * 1617787785200538278UL) + ((uint64_t)op[6] * 376288571705508531UL);
	tmp_q[5] = ((uint64_t)op[0] * 376288571705508531UL) + ((uint64_t)op[1] * 4882460954289740432UL) + ((uint64_t)op[2] * 259251780415428890UL) + ((uint64_t)op[3] * 15650751168344772518UL) + ((uint64_t)op[4] * 7517688434474851058UL) + ((uint64_t)op[5] * 17130256066536578734UL) + ((uint64_t)op[6] * 1617787785200538278UL);
	tmp_q[6] = ((uint64_t)op[0] * 1617787785200538278UL) + ((uint64_t)op[1] * 376288571705508531UL) + ((uint64_t)op[2] * 4882460954289740432UL) + ((uint64_t)op[3] * 259251780415428890UL) + ((uint64_t)op[4] * 15650751168344772518UL) + ((uint64_t)op[5] * 7517688434474851058UL) + ((uint64_t)op[6] * 17130256066536578734UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4467354304726L) + (-((int128)tmp_q[1] * 3501616326238L) + ((int128)tmp_q[2] * 3695037944646L) + ((int128)tmp_q[3] * 4204982265442L) + ((int128)tmp_q[4] * 1203312297178L) - ((int128)tmp_q[5] * 3494274517851L) + ((int128)tmp_q[6] * 2359912641556L));
	tmp_zero[1] = ((int128)tmp_q[0] * 2359912641556L) - ((int128)tmp_q[1] * 4467354304726L) + (-((int128)tmp_q[2] * 3501616326238L) + ((int128)tmp_q[3] * 3695037944646L) + ((int128)tmp_q[4] * 4204982265442L) + ((int128)tmp_q[5] * 1203312297178L) - ((int128)tmp_q[6] * 3494274517851L));
	tmp_zero[2] = -((int128)tmp_q[0] * 3494274517851L) + ((int128)tmp_q[1] * 2359912641556L) - ((int128)tmp_q[2] * 4467354304726L) + (-((int128)tmp_q[3] * 3501616326238L) + ((int128)tmp_q[4] * 3695037944646L) + ((int128)tmp_q[5] * 4204982265442L) + ((int128)tmp_q[6] * 1203312297178L));
	tmp_zero[3] = ((int128)tmp_q[0] * 1203312297178L) - ((int128)tmp_q[1] * 3494274517851L) + ((int128)tmp_q[2] * 2359912641556L) - ((int128)tmp_q[3] * 4467354304726L) + (-((int128)tmp_q[4] * 3501616326238L) + ((int128)tmp_q[5] * 3695037944646L) + ((int128)tmp_q[6] * 4204982265442L));
	tmp_zero[4] = ((int128)tmp_q[0] * 4204982265442L) + ((int128)tmp_q[1] * 1203312297178L) - ((int128)tmp_q[2] * 3494274517851L) + ((int128)tmp_q[3] * 2359912641556L) - ((int128)tmp_q[4] * 4467354304726L) + (-((int128)tmp_q[5] * 3501616326238L) + ((int128)tmp_q[6] * 3695037944646L));
	tmp_zero[5] = ((int128)tmp_q[0] * 3695037944646L) + ((int128)tmp_q[1] * 4204982265442L) + ((int128)tmp_q[2] * 1203312297178L) - ((int128)tmp_q[3] * 3494274517851L) + ((int128)tmp_q[4] * 2359912641556L) - ((int128)tmp_q[5] * 4467354304726L) - ((int128)tmp_q[6] * 3501616326238L);
	tmp_zero[6] = -((int128)tmp_q[0] * 3501616326238L) + ((int128)tmp_q[1] * 3695037944646L) + ((int128)tmp_q[2] * 4204982265442L) + ((int128)tmp_q[3] * 1203312297178L) - ((int128)tmp_q[4] * 3494274517851L) + ((int128)tmp_q[5] * 2359912641556L) - ((int128)tmp_q[6] * 4467354304726L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

