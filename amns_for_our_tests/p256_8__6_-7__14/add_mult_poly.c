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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7193667809138476485UL) + ((((uint64_t)op[1] * 17073580320375225296UL) + ((uint64_t)op[2] * 7858519171719075908UL) + ((uint64_t)op[3] * 13982030731897526122UL) + ((uint64_t)op[4] * 16190159106176905830UL) + ((uint64_t)op[5] * 8464029031799416346UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 8464029031799416346UL) + ((uint64_t)op[1] * 7193667809138476485UL) + ((((uint64_t)op[2] * 17073580320375225296UL) + ((uint64_t)op[3] * 7858519171719075908UL) + ((uint64_t)op[4] * 13982030731897526122UL) + ((uint64_t)op[5] * 16190159106176905830UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 16190159106176905830UL) + ((uint64_t)op[1] * 8464029031799416346UL) + ((uint64_t)op[2] * 7193667809138476485UL) + ((((uint64_t)op[3] * 17073580320375225296UL) + ((uint64_t)op[4] * 7858519171719075908UL) + ((uint64_t)op[5] * 13982030731897526122UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 13982030731897526122UL) + ((uint64_t)op[1] * 16190159106176905830UL) + ((uint64_t)op[2] * 8464029031799416346UL) + ((uint64_t)op[3] * 7193667809138476485UL) + ((((uint64_t)op[4] * 17073580320375225296UL) + ((uint64_t)op[5] * 7858519171719075908UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 7858519171719075908UL) + ((uint64_t)op[1] * 13982030731897526122UL) + ((uint64_t)op[2] * 16190159106176905830UL) + ((uint64_t)op[3] * 8464029031799416346UL) + ((uint64_t)op[4] * 7193667809138476485UL) + ((uint64_t)op[5] * 9612146273340284240UL);
	tmp_q[5] = ((uint64_t)op[0] * 17073580320375225296UL) + ((uint64_t)op[1] * 7858519171719075908UL) + ((uint64_t)op[2] * 13982030731897526122UL) + ((uint64_t)op[3] * 16190159106176905830UL) + ((uint64_t)op[4] * 8464029031799416346UL) + ((uint64_t)op[5] * 7193667809138476485UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2191763614071L) - ((-((int128)tmp_q[1] * 667156386056L) - ((int128)tmp_q[2] * 2321050196912L) - ((int128)tmp_q[3] * 2555109548622L) + ((int128)tmp_q[4] * 2010859924314L) + ((int128)tmp_q[5] * 5053910703658L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 5053910703658L) + ((int128)tmp_q[1] * 2191763614071L) - ((-((int128)tmp_q[2] * 667156386056L) - ((int128)tmp_q[3] * 2321050196912L) - ((int128)tmp_q[4] * 2555109548622L) + ((int128)tmp_q[5] * 2010859924314L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 2010859924314L) + ((int128)tmp_q[1] * 5053910703658L) + ((int128)tmp_q[2] * 2191763614071L) - ((-((int128)tmp_q[3] * 667156386056L) - ((int128)tmp_q[4] * 2321050196912L) - ((int128)tmp_q[5] * 2555109548622L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 2555109548622L) + ((int128)tmp_q[1] * 2010859924314L) + ((int128)tmp_q[2] * 5053910703658L) + ((int128)tmp_q[3] * 2191763614071L) - ((-((int128)tmp_q[4] * 667156386056L) - ((int128)tmp_q[5] * 2321050196912L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 2321050196912L) - ((int128)tmp_q[1] * 2555109548622L) + ((int128)tmp_q[2] * 2010859924314L) + ((int128)tmp_q[3] * 5053910703658L) + ((int128)tmp_q[4] * 2191763614071L) + ((int128)tmp_q[5] * 4670094702392L);
	tmp_zero[5] = -((int128)tmp_q[0] * 667156386056L) - ((int128)tmp_q[1] * 2321050196912L) - ((int128)tmp_q[2] * 2555109548622L) + ((int128)tmp_q[3] * 2010859924314L) + ((int128)tmp_q[4] * 5053910703658L) + ((int128)tmp_q[5] * 2191763614071L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

