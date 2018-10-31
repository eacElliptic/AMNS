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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16175329874982539360UL) + ((((uint64_t)op[1] * 760557835416487947UL) + ((uint64_t)op[2] * 10804605776845918714UL) + ((uint64_t)op[3] * 17709806165331583967UL) + ((uint64_t)op[4] * 14688734758557400897UL) + ((uint64_t)op[5] * 2663127379175128014UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 2663127379175128014UL) + ((uint64_t)op[1] * 16175329874982539360UL) + ((((uint64_t)op[2] * 760557835416487947UL) + ((uint64_t)op[3] * 10804605776845918714UL) + ((uint64_t)op[4] * 17709806165331583967UL) + ((uint64_t)op[5] * 14688734758557400897UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 14688734758557400897UL) + ((uint64_t)op[1] * 2663127379175128014UL) + ((uint64_t)op[2] * 16175329874982539360UL) + ((((uint64_t)op[3] * 760557835416487947UL) + ((uint64_t)op[4] * 10804605776845918714UL) + ((uint64_t)op[5] * 17709806165331583967UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 17709806165331583967UL) + ((uint64_t)op[1] * 14688734758557400897UL) + ((uint64_t)op[2] * 2663127379175128014UL) + ((uint64_t)op[3] * 16175329874982539360UL) + ((((uint64_t)op[4] * 760557835416487947UL) + ((uint64_t)op[5] * 10804605776845918714UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 10804605776845918714UL) + ((uint64_t)op[1] * 17709806165331583967UL) + ((uint64_t)op[2] * 14688734758557400897UL) + ((uint64_t)op[3] * 2663127379175128014UL) + ((uint64_t)op[4] * 16175329874982539360UL) + ((uint64_t)op[5] * 2281673506249463841UL);
	tmp_q[5] = ((uint64_t)op[0] * 760557835416487947UL) + ((uint64_t)op[1] * 10804605776845918714UL) + ((uint64_t)op[2] * 17709806165331583967UL) + ((uint64_t)op[3] * 14688734758557400897UL) + ((uint64_t)op[4] * 2663127379175128014UL) + ((uint64_t)op[5] * 16175329874982539360UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1064198513242L) + ((((int128)tmp_q[1] * 1442445741681L) + ((int128)tmp_q[2] * 1197115499460L) - ((int128)tmp_q[3] * 3236685154791L) + ((int128)tmp_q[4] * 2518895643207L) + ((int128)tmp_q[5] * 680725315090L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 680725315090L) - ((int128)tmp_q[1] * 1064198513242L) + ((((int128)tmp_q[2] * 1442445741681L) + ((int128)tmp_q[3] * 1197115499460L) - ((int128)tmp_q[4] * 3236685154791L) + ((int128)tmp_q[5] * 2518895643207L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2518895643207L) + ((int128)tmp_q[1] * 680725315090L) - ((int128)tmp_q[2] * 1064198513242L) + ((((int128)tmp_q[3] * 1442445741681L) + ((int128)tmp_q[4] * 1197115499460L) - ((int128)tmp_q[5] * 3236685154791L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 3236685154791L) + ((int128)tmp_q[1] * 2518895643207L) + ((int128)tmp_q[2] * 680725315090L) - ((int128)tmp_q[3] * 1064198513242L) + ((((int128)tmp_q[4] * 1442445741681L) + ((int128)tmp_q[5] * 1197115499460L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1197115499460L) - ((int128)tmp_q[1] * 3236685154791L) + ((int128)tmp_q[2] * 2518895643207L) + ((int128)tmp_q[3] * 680725315090L) - ((int128)tmp_q[4] * 1064198513242L) + ((int128)tmp_q[5] * 4327337225043L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1442445741681L) + ((int128)tmp_q[1] * 1197115499460L) - ((int128)tmp_q[2] * 3236685154791L) + ((int128)tmp_q[3] * 2518895643207L) + ((int128)tmp_q[4] * 680725315090L) - ((int128)tmp_q[5] * 1064198513242L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

