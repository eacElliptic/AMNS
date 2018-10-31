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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 120646286235574767UL) + ((((uint64_t)op[1] * 15559519765447490485UL) + ((uint64_t)op[2] * 14011069347539808682UL) + ((uint64_t)op[3] * 13593648029355976007UL) + ((uint64_t)op[4] * 17616924472544751521UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 17616924472544751521UL) + ((uint64_t)op[1] * 120646286235574767UL) + ((((uint64_t)op[2] * 15559519765447490485UL) + ((uint64_t)op[3] * 14011069347539808682UL) + ((uint64_t)op[4] * 13593648029355976007UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 13593648029355976007UL) + ((uint64_t)op[1] * 17616924472544751521UL) + ((uint64_t)op[2] * 120646286235574767UL) + ((((uint64_t)op[3] * 15559519765447490485UL) + ((uint64_t)op[4] * 14011069347539808682UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 14011069347539808682UL) + ((uint64_t)op[1] * 13593648029355976007UL) + ((uint64_t)op[2] * 17616924472544751521UL) + ((uint64_t)op[3] * 120646286235574767UL) + ((uint64_t)op[4] * 13795693681322614184UL);
	tmp_q[4] = ((uint64_t)op[0] * 15559519765447490485UL) + ((uint64_t)op[1] * 14011069347539808682UL) + ((uint64_t)op[2] * 13593648029355976007UL) + ((uint64_t)op[3] * 17616924472544751521UL) + ((uint64_t)op[4] * 120646286235574767UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2496629561169L) + ((((int128)tmp_q[1] * 3353436180928L) + ((int128)tmp_q[2] * 16467694261313L) - ((int128)tmp_q[3] * 25808585266320L) - ((int128)tmp_q[4] * 3420663450447L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 3420663450447L) + ((int128)tmp_q[1] * 2496629561169L) + ((((int128)tmp_q[2] * 3353436180928L) + ((int128)tmp_q[3] * 16467694261313L) - ((int128)tmp_q[4] * 25808585266320L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 25808585266320L) - ((int128)tmp_q[1] * 3420663450447L) + ((int128)tmp_q[2] * 2496629561169L) + ((((int128)tmp_q[3] * 3353436180928L) + ((int128)tmp_q[4] * 16467694261313L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 16467694261313L) - ((int128)tmp_q[1] * 25808585266320L) - ((int128)tmp_q[2] * 3420663450447L) + ((int128)tmp_q[3] * 2496629561169L) + ((int128)tmp_q[4] * 26827489447424L);
	tmp_zero[4] = ((int128)tmp_q[0] * 3353436180928L) + ((int128)tmp_q[1] * 16467694261313L) - ((int128)tmp_q[2] * 25808585266320L) - ((int128)tmp_q[3] * 3420663450447L) + ((int128)tmp_q[4] * 2496629561169L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

