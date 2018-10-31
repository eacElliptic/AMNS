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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 148171736175084429UL) + ((((uint64_t)op[1] * 14563528291195670774UL) + ((uint64_t)op[2] * 7218373514717132610UL) + ((uint64_t)op[3] * 17878730741261868319UL) + ((uint64_t)op[4] * 18421128221058774153UL) + ((uint64_t)op[5] * 7524722380641091606UL) + ((uint64_t)op[6] * 13536186683436562779UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 13536186683436562779UL) + ((uint64_t)op[1] * 148171736175084429UL) + ((((uint64_t)op[2] * 14563528291195670774UL) + ((uint64_t)op[3] * 7218373514717132610UL) + ((uint64_t)op[4] * 17878730741261868319UL) + ((uint64_t)op[5] * 18421128221058774153UL) + ((uint64_t)op[6] * 7524722380641091606UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 7524722380641091606UL) + ((uint64_t)op[1] * 13536186683436562779UL) + ((uint64_t)op[2] * 148171736175084429UL) + ((((uint64_t)op[3] * 14563528291195670774UL) + ((uint64_t)op[4] * 7218373514717132610UL) + ((uint64_t)op[5] * 17878730741261868319UL) + ((uint64_t)op[6] * 18421128221058774153UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 18421128221058774153UL) + ((uint64_t)op[1] * 7524722380641091606UL) + ((uint64_t)op[2] * 13536186683436562779UL) + ((uint64_t)op[3] * 148171736175084429UL) + ((((uint64_t)op[4] * 14563528291195670774UL) + ((uint64_t)op[5] * 7218373514717132610UL) + ((uint64_t)op[6] * 17878730741261868319UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 17878730741261868319UL) + ((uint64_t)op[1] * 18421128221058774153UL) + ((uint64_t)op[2] * 7524722380641091606UL) + ((uint64_t)op[3] * 13536186683436562779UL) + ((uint64_t)op[4] * 148171736175084429UL) + ((((uint64_t)op[5] * 14563528291195670774UL) + ((uint64_t)op[6] * 7218373514717132610UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 7218373514717132610UL) + ((uint64_t)op[1] * 17878730741261868319UL) + ((uint64_t)op[2] * 18421128221058774153UL) + ((uint64_t)op[3] * 7524722380641091606UL) + ((uint64_t)op[4] * 13536186683436562779UL) + ((uint64_t)op[5] * 148171736175084429UL) + ((uint64_t)op[6] * 7766431565027761684UL);
	tmp_q[6] = ((uint64_t)op[0] * 14563528291195670774UL) + ((uint64_t)op[1] * 7218373514717132610UL) + ((uint64_t)op[2] * 17878730741261868319UL) + ((uint64_t)op[3] * 18421128221058774153UL) + ((uint64_t)op[4] * 7524722380641091606UL) + ((uint64_t)op[5] * 13536186683436562779UL) + ((uint64_t)op[6] * 148171736175084429UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6140018443L) - ((-((int128)tmp_q[1] * 14229196395L) + ((int128)tmp_q[2] * 45796226980L) - ((int128)tmp_q[3] * 49067978434L) - ((int128)tmp_q[4] * 19316864156L) + ((int128)tmp_q[5] * 24749131621L) + ((int128)tmp_q[6] * 43443214065L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 43443214065L) - ((int128)tmp_q[1] * 6140018443L) - ((-((int128)tmp_q[2] * 14229196395L) + ((int128)tmp_q[3] * 45796226980L) - ((int128)tmp_q[4] * 49067978434L) - ((int128)tmp_q[5] * 19316864156L) + ((int128)tmp_q[6] * 24749131621L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 24749131621L) + ((int128)tmp_q[1] * 43443214065L) - ((int128)tmp_q[2] * 6140018443L) - ((-((int128)tmp_q[3] * 14229196395L) + ((int128)tmp_q[4] * 45796226980L) - ((int128)tmp_q[5] * 49067978434L) - ((int128)tmp_q[6] * 19316864156L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 19316864156L) + ((int128)tmp_q[1] * 24749131621L) + ((int128)tmp_q[2] * 43443214065L) - ((int128)tmp_q[3] * 6140018443L) - ((-((int128)tmp_q[4] * 14229196395L) + ((int128)tmp_q[5] * 45796226980L) - ((int128)tmp_q[6] * 49067978434L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 49067978434L) - ((int128)tmp_q[1] * 19316864156L) + ((int128)tmp_q[2] * 24749131621L) + ((int128)tmp_q[3] * 43443214065L) - ((int128)tmp_q[4] * 6140018443L) - ((-((int128)tmp_q[5] * 14229196395L) + ((int128)tmp_q[6] * 45796226980L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 45796226980L) - ((int128)tmp_q[1] * 49067978434L) - ((int128)tmp_q[2] * 19316864156L) + ((int128)tmp_q[3] * 24749131621L) + ((int128)tmp_q[4] * 43443214065L) - ((int128)tmp_q[5] * 6140018443L) + ((int128)tmp_q[6] * 28458392790L);
	tmp_zero[6] = -((int128)tmp_q[0] * 14229196395L) + ((int128)tmp_q[1] * 45796226980L) - ((int128)tmp_q[2] * 49067978434L) - ((int128)tmp_q[3] * 19316864156L) + ((int128)tmp_q[4] * 24749131621L) + ((int128)tmp_q[5] * 43443214065L) - ((int128)tmp_q[6] * 6140018443L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

