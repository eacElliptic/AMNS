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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10788131371398924611UL) + ((((uint64_t)op[1] * 6842951830157911777UL) + ((uint64_t)op[2] * 12574648594070799155UL) + ((uint64_t)op[3] * 2610384155089415020UL) + ((uint64_t)op[4] * 2054455067071932053UL) + ((uint64_t)op[5] * 12019196794661859017UL) + ((uint64_t)op[6] * 3051511843206794305UL) + ((uint64_t)op[7] * 16596432931194944870UL) + ((uint64_t)op[8] * 1173923864233943876UL) + ((uint64_t)op[9] * 11643937432567285189UL) + ((uint64_t)op[10] * 13492113065484167636UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 13492113065484167636UL) + ((uint64_t)op[1] * 10788131371398924611UL) + ((((uint64_t)op[2] * 6842951830157911777UL) + ((uint64_t)op[3] * 12574648594070799155UL) + ((uint64_t)op[4] * 2610384155089415020UL) + ((uint64_t)op[5] * 2054455067071932053UL) + ((uint64_t)op[6] * 12019196794661859017UL) + ((uint64_t)op[7] * 3051511843206794305UL) + ((uint64_t)op[8] * 16596432931194944870UL) + ((uint64_t)op[9] * 1173923864233943876UL) + ((uint64_t)op[10] * 11643937432567285189UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 11643937432567285189UL) + ((uint64_t)op[1] * 13492113065484167636UL) + ((uint64_t)op[2] * 10788131371398924611UL) + ((((uint64_t)op[3] * 6842951830157911777UL) + ((uint64_t)op[4] * 12574648594070799155UL) + ((uint64_t)op[5] * 2610384155089415020UL) + ((uint64_t)op[6] * 2054455067071932053UL) + ((uint64_t)op[7] * 12019196794661859017UL) + ((uint64_t)op[8] * 3051511843206794305UL) + ((uint64_t)op[9] * 16596432931194944870UL) + ((uint64_t)op[10] * 1173923864233943876UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 1173923864233943876UL) + ((uint64_t)op[1] * 11643937432567285189UL) + ((uint64_t)op[2] * 13492113065484167636UL) + ((uint64_t)op[3] * 10788131371398924611UL) + ((((uint64_t)op[4] * 6842951830157911777UL) + ((uint64_t)op[5] * 12574648594070799155UL) + ((uint64_t)op[6] * 2610384155089415020UL) + ((uint64_t)op[7] * 2054455067071932053UL) + ((uint64_t)op[8] * 12019196794661859017UL) + ((uint64_t)op[9] * 3051511843206794305UL) + ((uint64_t)op[10] * 16596432931194944870UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 16596432931194944870UL) + ((uint64_t)op[1] * 1173923864233943876UL) + ((uint64_t)op[2] * 11643937432567285189UL) + ((uint64_t)op[3] * 13492113065484167636UL) + ((uint64_t)op[4] * 10788131371398924611UL) + ((((uint64_t)op[5] * 6842951830157911777UL) + ((uint64_t)op[6] * 12574648594070799155UL) + ((uint64_t)op[7] * 2610384155089415020UL) + ((uint64_t)op[8] * 2054455067071932053UL) + ((uint64_t)op[9] * 12019196794661859017UL) + ((uint64_t)op[10] * 3051511843206794305UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 3051511843206794305UL) + ((uint64_t)op[1] * 16596432931194944870UL) + ((uint64_t)op[2] * 1173923864233943876UL) + ((uint64_t)op[3] * 11643937432567285189UL) + ((uint64_t)op[4] * 13492113065484167636UL) + ((uint64_t)op[5] * 10788131371398924611UL) + ((((uint64_t)op[6] * 6842951830157911777UL) + ((uint64_t)op[7] * 12574648594070799155UL) + ((uint64_t)op[8] * 2610384155089415020UL) + ((uint64_t)op[9] * 2054455067071932053UL) + ((uint64_t)op[10] * 12019196794661859017UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 12019196794661859017UL) + ((uint64_t)op[1] * 3051511843206794305UL) + ((uint64_t)op[2] * 16596432931194944870UL) + ((uint64_t)op[3] * 1173923864233943876UL) + ((uint64_t)op[4] * 11643937432567285189UL) + ((uint64_t)op[5] * 13492113065484167636UL) + ((uint64_t)op[6] * 10788131371398924611UL) + ((((uint64_t)op[7] * 6842951830157911777UL) + ((uint64_t)op[8] * 12574648594070799155UL) + ((uint64_t)op[9] * 2610384155089415020UL) + ((uint64_t)op[10] * 2054455067071932053UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 2054455067071932053UL) + ((uint64_t)op[1] * 12019196794661859017UL) + ((uint64_t)op[2] * 3051511843206794305UL) + ((uint64_t)op[3] * 16596432931194944870UL) + ((uint64_t)op[4] * 1173923864233943876UL) + ((uint64_t)op[5] * 11643937432567285189UL) + ((uint64_t)op[6] * 13492113065484167636UL) + ((uint64_t)op[7] * 10788131371398924611UL) + ((((uint64_t)op[8] * 6842951830157911777UL) + ((uint64_t)op[9] * 12574648594070799155UL) + ((uint64_t)op[10] * 2610384155089415020UL)) * 18446744073709551610);
	tmp_q[8] = ((uint64_t)op[0] * 2610384155089415020UL) + ((uint64_t)op[1] * 2054455067071932053UL) + ((uint64_t)op[2] * 12019196794661859017UL) + ((uint64_t)op[3] * 3051511843206794305UL) + ((uint64_t)op[4] * 16596432931194944870UL) + ((uint64_t)op[5] * 1173923864233943876UL) + ((uint64_t)op[6] * 11643937432567285189UL) + ((uint64_t)op[7] * 13492113065484167636UL) + ((uint64_t)op[8] * 10788131371398924611UL) + ((((uint64_t)op[9] * 6842951830157911777UL) + ((uint64_t)op[10] * 12574648594070799155UL)) * 18446744073709551610);
	tmp_q[9] = ((uint64_t)op[0] * 12574648594070799155UL) + ((uint64_t)op[1] * 2610384155089415020UL) + ((uint64_t)op[2] * 2054455067071932053UL) + ((uint64_t)op[3] * 12019196794661859017UL) + ((uint64_t)op[4] * 3051511843206794305UL) + ((uint64_t)op[5] * 16596432931194944870UL) + ((uint64_t)op[6] * 1173923864233943876UL) + ((uint64_t)op[7] * 11643937432567285189UL) + ((uint64_t)op[8] * 13492113065484167636UL) + ((uint64_t)op[9] * 10788131371398924611UL) + ((uint64_t)op[10] * 14282521240181184186UL);
	tmp_q[10] = ((uint64_t)op[0] * 6842951830157911777UL) + ((uint64_t)op[1] * 12574648594070799155UL) + ((uint64_t)op[2] * 2610384155089415020UL) + ((uint64_t)op[3] * 2054455067071932053UL) + ((uint64_t)op[4] * 12019196794661859017UL) + ((uint64_t)op[5] * 3051511843206794305UL) + ((uint64_t)op[6] * 16596432931194944870UL) + ((uint64_t)op[7] * 1173923864233943876UL) + ((uint64_t)op[8] * 11643937432567285189UL) + ((uint64_t)op[9] * 13492113065484167636UL) + ((uint64_t)op[10] * 10788131371398924611UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 62758010530901L) - ((((int128)tmp_q[1] * 16819734153378L) + ((int128)tmp_q[2] * 18523334347440L) - ((int128)tmp_q[3] * 67196767386583L) + ((int128)tmp_q[4] * 58890218523927L) - ((int128)tmp_q[5] * 72838850341670L) - ((int128)tmp_q[6] * 79445255805001L) - ((int128)tmp_q[7] * 51014017380057L) + ((int128)tmp_q[8] * 87346431091274L) - ((int128)tmp_q[9] * 44862091927883L) - ((int128)tmp_q[10] * 42760779907230L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 42760779907230L) - ((int128)tmp_q[1] * 62758010530901L) - ((((int128)tmp_q[2] * 16819734153378L) + ((int128)tmp_q[3] * 18523334347440L) - ((int128)tmp_q[4] * 67196767386583L) + ((int128)tmp_q[5] * 58890218523927L) - ((int128)tmp_q[6] * 72838850341670L) - ((int128)tmp_q[7] * 79445255805001L) - ((int128)tmp_q[8] * 51014017380057L) + ((int128)tmp_q[9] * 87346431091274L) - ((int128)tmp_q[10] * 44862091927883L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 44862091927883L) - ((int128)tmp_q[1] * 42760779907230L) - ((int128)tmp_q[2] * 62758010530901L) - ((((int128)tmp_q[3] * 16819734153378L) + ((int128)tmp_q[4] * 18523334347440L) - ((int128)tmp_q[5] * 67196767386583L) + ((int128)tmp_q[6] * 58890218523927L) - ((int128)tmp_q[7] * 72838850341670L) - ((int128)tmp_q[8] * 79445255805001L) - ((int128)tmp_q[9] * 51014017380057L) + ((int128)tmp_q[10] * 87346431091274L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 87346431091274L) - ((int128)tmp_q[1] * 44862091927883L) - ((int128)tmp_q[2] * 42760779907230L) - ((int128)tmp_q[3] * 62758010530901L) - ((((int128)tmp_q[4] * 16819734153378L) + ((int128)tmp_q[5] * 18523334347440L) - ((int128)tmp_q[6] * 67196767386583L) + ((int128)tmp_q[7] * 58890218523927L) - ((int128)tmp_q[8] * 72838850341670L) - ((int128)tmp_q[9] * 79445255805001L) - ((int128)tmp_q[10] * 51014017380057L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 51014017380057L) + ((int128)tmp_q[1] * 87346431091274L) - ((int128)tmp_q[2] * 44862091927883L) - ((int128)tmp_q[3] * 42760779907230L) - ((int128)tmp_q[4] * 62758010530901L) - ((((int128)tmp_q[5] * 16819734153378L) + ((int128)tmp_q[6] * 18523334347440L) - ((int128)tmp_q[7] * 67196767386583L) + ((int128)tmp_q[8] * 58890218523927L) - ((int128)tmp_q[9] * 72838850341670L) - ((int128)tmp_q[10] * 79445255805001L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 79445255805001L) - ((int128)tmp_q[1] * 51014017380057L) + ((int128)tmp_q[2] * 87346431091274L) - ((int128)tmp_q[3] * 44862091927883L) - ((int128)tmp_q[4] * 42760779907230L) - ((int128)tmp_q[5] * 62758010530901L) - ((((int128)tmp_q[6] * 16819734153378L) + ((int128)tmp_q[7] * 18523334347440L) - ((int128)tmp_q[8] * 67196767386583L) + ((int128)tmp_q[9] * 58890218523927L) - ((int128)tmp_q[10] * 72838850341670L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 72838850341670L) - ((int128)tmp_q[1] * 79445255805001L) - ((int128)tmp_q[2] * 51014017380057L) + ((int128)tmp_q[3] * 87346431091274L) - ((int128)tmp_q[4] * 44862091927883L) - ((int128)tmp_q[5] * 42760779907230L) - ((int128)tmp_q[6] * 62758010530901L) - ((((int128)tmp_q[7] * 16819734153378L) + ((int128)tmp_q[8] * 18523334347440L) - ((int128)tmp_q[9] * 67196767386583L) + ((int128)tmp_q[10] * 58890218523927L)) * 6);
	tmp_zero[7] = ((int128)tmp_q[0] * 58890218523927L) - ((int128)tmp_q[1] * 72838850341670L) - ((int128)tmp_q[2] * 79445255805001L) - ((int128)tmp_q[3] * 51014017380057L) + ((int128)tmp_q[4] * 87346431091274L) - ((int128)tmp_q[5] * 44862091927883L) - ((int128)tmp_q[6] * 42760779907230L) - ((int128)tmp_q[7] * 62758010530901L) - ((((int128)tmp_q[8] * 16819734153378L) + ((int128)tmp_q[9] * 18523334347440L) - ((int128)tmp_q[10] * 67196767386583L)) * 6);
	tmp_zero[8] = -((int128)tmp_q[0] * 67196767386583L) + ((int128)tmp_q[1] * 58890218523927L) - ((int128)tmp_q[2] * 72838850341670L) - ((int128)tmp_q[3] * 79445255805001L) - ((int128)tmp_q[4] * 51014017380057L) + ((int128)tmp_q[5] * 87346431091274L) - ((int128)tmp_q[6] * 44862091927883L) - ((int128)tmp_q[7] * 42760779907230L) - ((int128)tmp_q[8] * 62758010530901L) - ((((int128)tmp_q[9] * 16819734153378L) + ((int128)tmp_q[10] * 18523334347440L)) * 6);
	tmp_zero[9] = ((int128)tmp_q[0] * 18523334347440L) - ((int128)tmp_q[1] * 67196767386583L) + ((int128)tmp_q[2] * 58890218523927L) - ((int128)tmp_q[3] * 72838850341670L) - ((int128)tmp_q[4] * 79445255805001L) - ((int128)tmp_q[5] * 51014017380057L) + ((int128)tmp_q[6] * 87346431091274L) - ((int128)tmp_q[7] * 44862091927883L) - ((int128)tmp_q[8] * 42760779907230L) - ((int128)tmp_q[9] * 62758010530901L) - ((int128)tmp_q[10] * 100918404920268L);
	tmp_zero[10] = ((int128)tmp_q[0] * 16819734153378L) + ((int128)tmp_q[1] * 18523334347440L) - ((int128)tmp_q[2] * 67196767386583L) + ((int128)tmp_q[3] * 58890218523927L) - ((int128)tmp_q[4] * 72838850341670L) - ((int128)tmp_q[5] * 79445255805001L) - ((int128)tmp_q[6] * 51014017380057L) + ((int128)tmp_q[7] * 87346431091274L) - ((int128)tmp_q[8] * 44862091927883L) - ((int128)tmp_q[9] * 42760779907230L) - ((int128)tmp_q[10] * 62758010530901L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
}

