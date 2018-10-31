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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7694687945627488685UL) + ((((uint64_t)op[1] * 3166470379756018677UL) + ((uint64_t)op[2] * 13842017160969907978UL) + ((uint64_t)op[3] * 7572763781264286794UL) + ((uint64_t)op[4] * 9987016803379322603UL) + ((uint64_t)op[5] * 17085508628062399020UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 17085508628062399020UL) + ((uint64_t)op[1] * 7694687945627488685UL) + ((((uint64_t)op[2] * 3166470379756018677UL) + ((uint64_t)op[3] * 13842017160969907978UL) + ((uint64_t)op[4] * 7572763781264286794UL) + ((uint64_t)op[5] * 9987016803379322603UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 9987016803379322603UL) + ((uint64_t)op[1] * 17085508628062399020UL) + ((uint64_t)op[2] * 7694687945627488685UL) + ((((uint64_t)op[3] * 3166470379756018677UL) + ((uint64_t)op[4] * 13842017160969907978UL) + ((uint64_t)op[5] * 7572763781264286794UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7572763781264286794UL) + ((uint64_t)op[1] * 9987016803379322603UL) + ((uint64_t)op[2] * 17085508628062399020UL) + ((uint64_t)op[3] * 7694687945627488685UL) + ((((uint64_t)op[4] * 3166470379756018677UL) + ((uint64_t)op[5] * 13842017160969907978UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 13842017160969907978UL) + ((uint64_t)op[1] * 7572763781264286794UL) + ((uint64_t)op[2] * 9987016803379322603UL) + ((uint64_t)op[3] * 17085508628062399020UL) + ((uint64_t)op[4] * 7694687945627488685UL) + ((uint64_t)op[5] * 8947332934441495585UL);
	tmp_q[5] = ((uint64_t)op[0] * 3166470379756018677UL) + ((uint64_t)op[1] * 13842017160969907978UL) + ((uint64_t)op[2] * 7572763781264286794UL) + ((uint64_t)op[3] * 9987016803379322603UL) + ((uint64_t)op[4] * 17085508628062399020UL) + ((uint64_t)op[5] * 7694687945627488685UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 905234880197L) - ((-((int128)tmp_q[1] * 545769043655L) - ((int128)tmp_q[2] * 2752504631420L) - ((int128)tmp_q[3] * 1232196310098L) - ((int128)tmp_q[4] * 3562321580221L) - ((int128)tmp_q[5] * 3194110223610L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 3194110223610L) - ((int128)tmp_q[1] * 905234880197L) - ((-((int128)tmp_q[2] * 545769043655L) - ((int128)tmp_q[3] * 2752504631420L) - ((int128)tmp_q[4] * 1232196310098L) - ((int128)tmp_q[5] * 3562321580221L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 3562321580221L) - ((int128)tmp_q[1] * 3194110223610L) - ((int128)tmp_q[2] * 905234880197L) - ((-((int128)tmp_q[3] * 545769043655L) - ((int128)tmp_q[4] * 2752504631420L) - ((int128)tmp_q[5] * 1232196310098L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 1232196310098L) - ((int128)tmp_q[1] * 3562321580221L) - ((int128)tmp_q[2] * 3194110223610L) - ((int128)tmp_q[3] * 905234880197L) - ((-((int128)tmp_q[4] * 545769043655L) - ((int128)tmp_q[5] * 2752504631420L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 2752504631420L) - ((int128)tmp_q[1] * 1232196310098L) - ((int128)tmp_q[2] * 3562321580221L) - ((int128)tmp_q[3] * 3194110223610L) - ((int128)tmp_q[4] * 905234880197L) + ((int128)tmp_q[5] * 1637307130965L);
	tmp_zero[5] = -((int128)tmp_q[0] * 545769043655L) - ((int128)tmp_q[1] * 2752504631420L) - ((int128)tmp_q[2] * 1232196310098L) - ((int128)tmp_q[3] * 3562321580221L) - ((int128)tmp_q[4] * 3194110223610L) - ((int128)tmp_q[5] * 905234880197L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

