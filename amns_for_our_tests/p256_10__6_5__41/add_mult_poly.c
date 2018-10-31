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
	tmp_q[0] = ((uint64_t)op[0] * 6639265928802566489UL) + ((((uint64_t)op[1] * 2342670993397580548UL) + ((uint64_t)op[2] * 14542742547969863887UL) + ((uint64_t)op[3] * 1556123459629021822UL) + ((uint64_t)op[4] * 5330418212788200536UL) + ((uint64_t)op[5] * 13546994818982681979UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 13546994818982681979UL) + ((uint64_t)op[1] * 6639265928802566489UL) + ((((uint64_t)op[2] * 2342670993397580548UL) + ((uint64_t)op[3] * 14542742547969863887UL) + ((uint64_t)op[4] * 1556123459629021822UL) + ((uint64_t)op[5] * 5330418212788200536UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 5330418212788200536UL) + ((uint64_t)op[1] * 13546994818982681979UL) + ((uint64_t)op[2] * 6639265928802566489UL) + ((((uint64_t)op[3] * 2342670993397580548UL) + ((uint64_t)op[4] * 14542742547969863887UL) + ((uint64_t)op[5] * 1556123459629021822UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 1556123459629021822UL) + ((uint64_t)op[1] * 5330418212788200536UL) + ((uint64_t)op[2] * 13546994818982681979UL) + ((uint64_t)op[3] * 6639265928802566489UL) + ((((uint64_t)op[4] * 2342670993397580548UL) + ((uint64_t)op[5] * 14542742547969863887UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 14542742547969863887UL) + ((uint64_t)op[1] * 1556123459629021822UL) + ((uint64_t)op[2] * 5330418212788200536UL) + ((uint64_t)op[3] * 13546994818982681979UL) + ((uint64_t)op[4] * 6639265928802566489UL) + ((uint64_t)op[5] * 11713354966987902740UL);
	tmp_q[5] = ((uint64_t)op[0] * 2342670993397580548UL) + ((uint64_t)op[1] * 14542742547969863887UL) + ((uint64_t)op[2] * 1556123459629021822UL) + ((uint64_t)op[3] * 5330418212788200536UL) + ((uint64_t)op[4] * 13546994818982681979UL) + ((uint64_t)op[5] * 6639265928802566489UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2569087564597L) + ((((int128)tmp_q[1] * 808644078550L) + ((int128)tmp_q[2] * 2949784727567L) - ((int128)tmp_q[3] * 2242411200010L) + ((int128)tmp_q[4] * 3771550552862L) - ((int128)tmp_q[5] * 2725236743157L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 2725236743157L) - ((int128)tmp_q[1] * 2569087564597L) + ((((int128)tmp_q[2] * 808644078550L) + ((int128)tmp_q[3] * 2949784727567L) - ((int128)tmp_q[4] * 2242411200010L) + ((int128)tmp_q[5] * 3771550552862L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 3771550552862L) - ((int128)tmp_q[1] * 2725236743157L) - ((int128)tmp_q[2] * 2569087564597L) + ((((int128)tmp_q[3] * 808644078550L) + ((int128)tmp_q[4] * 2949784727567L) - ((int128)tmp_q[5] * 2242411200010L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2242411200010L) + ((int128)tmp_q[1] * 3771550552862L) - ((int128)tmp_q[2] * 2725236743157L) - ((int128)tmp_q[3] * 2569087564597L) + ((((int128)tmp_q[4] * 808644078550L) + ((int128)tmp_q[5] * 2949784727567L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 2949784727567L) - ((int128)tmp_q[1] * 2242411200010L) + ((int128)tmp_q[2] * 3771550552862L) - ((int128)tmp_q[3] * 2725236743157L) - ((int128)tmp_q[4] * 2569087564597L) + ((int128)tmp_q[5] * 4043220392750L);
	tmp_zero[5] = ((int128)tmp_q[0] * 808644078550L) + ((int128)tmp_q[1] * 2949784727567L) - ((int128)tmp_q[2] * 2242411200010L) + ((int128)tmp_q[3] * 3771550552862L) - ((int128)tmp_q[4] * 2725236743157L) - ((int128)tmp_q[5] * 2569087564597L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

