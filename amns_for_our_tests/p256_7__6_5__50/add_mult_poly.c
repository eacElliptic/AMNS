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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18277734048191897745UL) + ((((uint64_t)op[1] * 9172271497221285080UL) + ((uint64_t)op[2] * 6148769694441883290UL) + ((uint64_t)op[3] * 17364025244578405991UL) + ((uint64_t)op[4] * 13210088814633206531UL) + ((uint64_t)op[5] * 17619905913632512702UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 17619905913632512702UL) + ((uint64_t)op[1] * 18277734048191897745UL) + ((((uint64_t)op[2] * 9172271497221285080UL) + ((uint64_t)op[3] * 6148769694441883290UL) + ((uint64_t)op[4] * 17364025244578405991UL) + ((uint64_t)op[5] * 13210088814633206531UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13210088814633206531UL) + ((uint64_t)op[1] * 17619905913632512702UL) + ((uint64_t)op[2] * 18277734048191897745UL) + ((((uint64_t)op[3] * 9172271497221285080UL) + ((uint64_t)op[4] * 6148769694441883290UL) + ((uint64_t)op[5] * 17364025244578405991UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 17364025244578405991UL) + ((uint64_t)op[1] * 13210088814633206531UL) + ((uint64_t)op[2] * 17619905913632512702UL) + ((uint64_t)op[3] * 18277734048191897745UL) + ((((uint64_t)op[4] * 9172271497221285080UL) + ((uint64_t)op[5] * 6148769694441883290UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 6148769694441883290UL) + ((uint64_t)op[1] * 17364025244578405991UL) + ((uint64_t)op[2] * 13210088814633206531UL) + ((uint64_t)op[3] * 17619905913632512702UL) + ((uint64_t)op[4] * 18277734048191897745UL) + ((uint64_t)op[5] * 8967869338687322168UL);
	tmp_q[5] = ((uint64_t)op[0] * 9172271497221285080UL) + ((uint64_t)op[1] * 6148769694441883290UL) + ((uint64_t)op[2] * 17364025244578405991UL) + ((uint64_t)op[3] * 13210088814633206531UL) + ((uint64_t)op[4] * 17619905913632512702UL) + ((uint64_t)op[5] * 18277734048191897745UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1812166793488L) + ((-((int128)tmp_q[1] * 4162929281981L) + ((int128)tmp_q[2] * 1125316642889L) + ((int128)tmp_q[3] * 3455774510068L) + ((int128)tmp_q[4] * 935591577589L) + ((int128)tmp_q[5] * 3247794992472L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 3247794992472L) + ((int128)tmp_q[1] * 1812166793488L) + ((-((int128)tmp_q[2] * 4162929281981L) + ((int128)tmp_q[3] * 1125316642889L) + ((int128)tmp_q[4] * 3455774510068L) + ((int128)tmp_q[5] * 935591577589L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 935591577589L) + ((int128)tmp_q[1] * 3247794992472L) + ((int128)tmp_q[2] * 1812166793488L) + ((-((int128)tmp_q[3] * 4162929281981L) + ((int128)tmp_q[4] * 1125316642889L) + ((int128)tmp_q[5] * 3455774510068L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 3455774510068L) + ((int128)tmp_q[1] * 935591577589L) + ((int128)tmp_q[2] * 3247794992472L) + ((int128)tmp_q[3] * 1812166793488L) + ((-((int128)tmp_q[4] * 4162929281981L) + ((int128)tmp_q[5] * 1125316642889L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1125316642889L) + ((int128)tmp_q[1] * 3455774510068L) + ((int128)tmp_q[2] * 935591577589L) + ((int128)tmp_q[3] * 3247794992472L) + ((int128)tmp_q[4] * 1812166793488L) - ((int128)tmp_q[5] * 20814646409905L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4162929281981L) + ((int128)tmp_q[1] * 1125316642889L) + ((int128)tmp_q[2] * 3455774510068L) + ((int128)tmp_q[3] * 935591577589L) + ((int128)tmp_q[4] * 3247794992472L) + ((int128)tmp_q[5] * 1812166793488L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

