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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14294782742149001872UL) + ((((uint64_t)op[1] * 16601090425410334486UL) + ((uint64_t)op[2] * 4975286061077153486UL) + ((uint64_t)op[3] * 5585813703495424764UL) + ((uint64_t)op[4] * 12327544035597330695UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 12327544035597330695UL) + ((uint64_t)op[1] * 14294782742149001872UL) + ((((uint64_t)op[2] * 16601090425410334486UL) + ((uint64_t)op[3] * 4975286061077153486UL) + ((uint64_t)op[4] * 5585813703495424764UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 5585813703495424764UL) + ((uint64_t)op[1] * 12327544035597330695UL) + ((uint64_t)op[2] * 14294782742149001872UL) + ((((uint64_t)op[3] * 16601090425410334486UL) + ((uint64_t)op[4] * 4975286061077153486UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 4975286061077153486UL) + ((uint64_t)op[1] * 5585813703495424764UL) + ((uint64_t)op[2] * 12327544035597330695UL) + ((uint64_t)op[3] * 14294782742149001872UL) + ((uint64_t)op[4] * 9218475832213465966UL);
	tmp_q[4] = ((uint64_t)op[0] * 16601090425410334486UL) + ((uint64_t)op[1] * 4975286061077153486UL) + ((uint64_t)op[2] * 5585813703495424764UL) + ((uint64_t)op[3] * 12327544035597330695UL) + ((uint64_t)op[4] * 14294782742149001872UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 260064863768L) + ((-((int128)tmp_q[1] * 45765687827L) + ((int128)tmp_q[2] * 120879890684L) + ((int128)tmp_q[3] * 234788721950L) - ((int128)tmp_q[4] * 15437595834L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 15437595834L) - ((int128)tmp_q[1] * 260064863768L) + ((-((int128)tmp_q[2] * 45765687827L) + ((int128)tmp_q[3] * 120879890684L) + ((int128)tmp_q[4] * 234788721950L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 234788721950L) - ((int128)tmp_q[1] * 15437595834L) - ((int128)tmp_q[2] * 260064863768L) + ((-((int128)tmp_q[3] * 45765687827L) + ((int128)tmp_q[4] * 120879890684L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 120879890684L) + ((int128)tmp_q[1] * 234788721950L) - ((int128)tmp_q[2] * 15437595834L) - ((int128)tmp_q[3] * 260064863768L) - ((int128)tmp_q[4] * 228828439135L);
	tmp_zero[4] = -((int128)tmp_q[0] * 45765687827L) + ((int128)tmp_q[1] * 120879890684L) + ((int128)tmp_q[2] * 234788721950L) - ((int128)tmp_q[3] * 15437595834L) - ((int128)tmp_q[4] * 260064863768L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

