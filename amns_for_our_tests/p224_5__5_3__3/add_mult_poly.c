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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17354937965552512765UL) + ((((uint64_t)op[1] * 5201819169796114696UL) + ((uint64_t)op[2] * 10610118658184290042UL) + ((uint64_t)op[3] * 7126930884540692456UL) + ((uint64_t)op[4] * 108426929825813072UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 108426929825813072UL) + ((uint64_t)op[1] * 17354937965552512765UL) + ((((uint64_t)op[2] * 5201819169796114696UL) + ((uint64_t)op[3] * 10610118658184290042UL) + ((uint64_t)op[4] * 7126930884540692456UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 7126930884540692456UL) + ((uint64_t)op[1] * 108426929825813072UL) + ((uint64_t)op[2] * 17354937965552512765UL) + ((((uint64_t)op[3] * 5201819169796114696UL) + ((uint64_t)op[4] * 10610118658184290042UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 10610118658184290042UL) + ((uint64_t)op[1] * 7126930884540692456UL) + ((uint64_t)op[2] * 108426929825813072UL) + ((uint64_t)op[3] * 17354937965552512765UL) + ((uint64_t)op[4] * 15605457509388344088UL);
	tmp_q[4] = ((uint64_t)op[0] * 5201819169796114696UL) + ((uint64_t)op[1] * 10610118658184290042UL) + ((uint64_t)op[2] * 7126930884540692456UL) + ((uint64_t)op[3] * 108426929825813072UL) + ((uint64_t)op[4] * 17354937965552512765UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3610934596597L) + ((-((int128)tmp_q[1] * 13158711525696L) - ((int128)tmp_q[2] * 7030613303094L) - ((int128)tmp_q[3] * 13991443173320L) - ((int128)tmp_q[4] * 13093159986476L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 13093159986476L) - ((int128)tmp_q[1] * 3610934596597L) + ((-((int128)tmp_q[2] * 13158711525696L) - ((int128)tmp_q[3] * 7030613303094L) - ((int128)tmp_q[4] * 13991443173320L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 13991443173320L) - ((int128)tmp_q[1] * 13093159986476L) - ((int128)tmp_q[2] * 3610934596597L) + ((-((int128)tmp_q[3] * 13158711525696L) - ((int128)tmp_q[4] * 7030613303094L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 7030613303094L) - ((int128)tmp_q[1] * 13991443173320L) - ((int128)tmp_q[2] * 13093159986476L) - ((int128)tmp_q[3] * 3610934596597L) - ((int128)tmp_q[4] * 39476134577088L);
	tmp_zero[4] = -((int128)tmp_q[0] * 13158711525696L) - ((int128)tmp_q[1] * 7030613303094L) - ((int128)tmp_q[2] * 13991443173320L) - ((int128)tmp_q[3] * 13093159986476L) - ((int128)tmp_q[4] * 3610934596597L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

