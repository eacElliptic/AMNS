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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7519537181221375271UL) + ((((uint64_t)op[1] * 4550051389112940595UL) + ((uint64_t)op[2] * 18146727287301539829UL) + ((uint64_t)op[3] * 14382186944855523058UL) + ((uint64_t)op[4] * 14481763294906836951UL) + ((uint64_t)op[5] * 17843442298706223393UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 17843442298706223393UL) + ((uint64_t)op[1] * 7519537181221375271UL) + ((((uint64_t)op[2] * 4550051389112940595UL) + ((uint64_t)op[3] * 18146727287301539829UL) + ((uint64_t)op[4] * 14382186944855523058UL) + ((uint64_t)op[5] * 14481763294906836951UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 14481763294906836951UL) + ((uint64_t)op[1] * 17843442298706223393UL) + ((uint64_t)op[2] * 7519537181221375271UL) + ((((uint64_t)op[3] * 4550051389112940595UL) + ((uint64_t)op[4] * 18146727287301539829UL) + ((uint64_t)op[5] * 14382186944855523058UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 14382186944855523058UL) + ((uint64_t)op[1] * 14481763294906836951UL) + ((uint64_t)op[2] * 17843442298706223393UL) + ((uint64_t)op[3] * 7519537181221375271UL) + ((((uint64_t)op[4] * 4550051389112940595UL) + ((uint64_t)op[5] * 18146727287301539829UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 18146727287301539829UL) + ((uint64_t)op[1] * 14382186944855523058UL) + ((uint64_t)op[2] * 14481763294906836951UL) + ((uint64_t)op[3] * 17843442298706223393UL) + ((uint64_t)op[4] * 7519537181221375271UL) + ((uint64_t)op[5] * 246538517257789236UL);
	tmp_q[5] = ((uint64_t)op[0] * 4550051389112940595UL) + ((uint64_t)op[1] * 18146727287301539829UL) + ((uint64_t)op[2] * 14382186944855523058UL) + ((uint64_t)op[3] * 14481763294906836951UL) + ((uint64_t)op[4] * 17843442298706223393UL) + ((uint64_t)op[5] * 7519537181221375271UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3916861857559L) - ((-((int128)tmp_q[1] * 1038324134537L) - ((int128)tmp_q[2] * 1852897410472L) - ((int128)tmp_q[3] * 4785269274975L) + ((int128)tmp_q[4] * 509287747316L) + ((int128)tmp_q[5] * 678974175017L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 678974175017L) - ((int128)tmp_q[1] * 3916861857559L) - ((-((int128)tmp_q[2] * 1038324134537L) - ((int128)tmp_q[3] * 1852897410472L) - ((int128)tmp_q[4] * 4785269274975L) + ((int128)tmp_q[5] * 509287747316L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 509287747316L) + ((int128)tmp_q[1] * 678974175017L) - ((int128)tmp_q[2] * 3916861857559L) - ((-((int128)tmp_q[3] * 1038324134537L) - ((int128)tmp_q[4] * 1852897410472L) - ((int128)tmp_q[5] * 4785269274975L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 4785269274975L) + ((int128)tmp_q[1] * 509287747316L) + ((int128)tmp_q[2] * 678974175017L) - ((int128)tmp_q[3] * 3916861857559L) - ((-((int128)tmp_q[4] * 1038324134537L) - ((int128)tmp_q[5] * 1852897410472L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 1852897410472L) - ((int128)tmp_q[1] * 4785269274975L) + ((int128)tmp_q[2] * 509287747316L) + ((int128)tmp_q[3] * 678974175017L) - ((int128)tmp_q[4] * 3916861857559L) + ((int128)tmp_q[5] * 4153296538148L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1038324134537L) - ((int128)tmp_q[1] * 1852897410472L) - ((int128)tmp_q[2] * 4785269274975L) + ((int128)tmp_q[3] * 509287747316L) + ((int128)tmp_q[4] * 678974175017L) - ((int128)tmp_q[5] * 3916861857559L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

