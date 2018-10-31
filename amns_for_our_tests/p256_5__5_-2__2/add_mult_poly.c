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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2373467291411748709UL) + ((((uint64_t)op[1] * 13398378273852063970UL) + ((uint64_t)op[2] * 4447930887779008950UL) + ((uint64_t)op[3] * 9591434027058088736UL) + ((uint64_t)op[4] * 6513448271990229040UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 6513448271990229040UL) + ((uint64_t)op[1] * 2373467291411748709UL) + ((((uint64_t)op[2] * 13398378273852063970UL) + ((uint64_t)op[3] * 4447930887779008950UL) + ((uint64_t)op[4] * 9591434027058088736UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 9591434027058088736UL) + ((uint64_t)op[1] * 6513448271990229040UL) + ((uint64_t)op[2] * 2373467291411748709UL) + ((((uint64_t)op[3] * 13398378273852063970UL) + ((uint64_t)op[4] * 4447930887779008950UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 4447930887779008950UL) + ((uint64_t)op[1] * 9591434027058088736UL) + ((uint64_t)op[2] * 6513448271990229040UL) + ((uint64_t)op[3] * 2373467291411748709UL) + ((uint64_t)op[4] * 10096731599714975292UL);
	tmp_q[4] = ((uint64_t)op[0] * 13398378273852063970UL) + ((uint64_t)op[1] * 4447930887779008950UL) + ((uint64_t)op[2] * 9591434027058088736UL) + ((uint64_t)op[3] * 6513448271990229040UL) + ((uint64_t)op[4] * 2373467291411748709UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 260012889325683L) - ((((int128)tmp_q[1] * 398496247683394L) - ((int128)tmp_q[2] * 197340263445682L) + ((int128)tmp_q[3] * 1928193251354608L) - ((int128)tmp_q[4] * 102611367928392L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 102611367928392L) + ((int128)tmp_q[1] * 260012889325683L) - ((((int128)tmp_q[2] * 398496247683394L) - ((int128)tmp_q[3] * 197340263445682L) + ((int128)tmp_q[4] * 1928193251354608L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1928193251354608L) - ((int128)tmp_q[1] * 102611367928392L) + ((int128)tmp_q[2] * 260012889325683L) - ((((int128)tmp_q[3] * 398496247683394L) - ((int128)tmp_q[4] * 197340263445682L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 197340263445682L) + ((int128)tmp_q[1] * 1928193251354608L) - ((int128)tmp_q[2] * 102611367928392L) + ((int128)tmp_q[3] * 260012889325683L) - ((int128)tmp_q[4] * 796992495366788L);
	tmp_zero[4] = ((int128)tmp_q[0] * 398496247683394L) - ((int128)tmp_q[1] * 197340263445682L) + ((int128)tmp_q[2] * 1928193251354608L) - ((int128)tmp_q[3] * 102611367928392L) + ((int128)tmp_q[4] * 260012889325683L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

