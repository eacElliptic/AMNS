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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18415531108464612935UL) + ((((uint64_t)op[1] * 7349181473276032313UL) + ((uint64_t)op[2] * 7868938781831905001UL) + ((uint64_t)op[3] * 14206693913815609646UL) + ((uint64_t)op[4] * 13823895817601966885UL) + ((uint64_t)op[5] * 10182052057624725370UL) + ((uint64_t)op[6] * 18031788028952716542UL) + ((uint64_t)op[7] * 13531614810474583979UL) + ((uint64_t)op[8] * 8664828748238991555UL) + ((uint64_t)op[9] * 14995837233037358424UL) + ((uint64_t)op[10] * 7445765303064700292UL) + ((uint64_t)op[11] * 1860807927596649322UL)) * 14);
	tmp_q[1] = ((uint64_t)op[0] * 1860807927596649322UL) + ((uint64_t)op[1] * 18415531108464612935UL) + ((((uint64_t)op[2] * 7349181473276032313UL) + ((uint64_t)op[3] * 7868938781831905001UL) + ((uint64_t)op[4] * 14206693913815609646UL) + ((uint64_t)op[5] * 13823895817601966885UL) + ((uint64_t)op[6] * 10182052057624725370UL) + ((uint64_t)op[7] * 18031788028952716542UL) + ((uint64_t)op[8] * 13531614810474583979UL) + ((uint64_t)op[9] * 8664828748238991555UL) + ((uint64_t)op[10] * 14995837233037358424UL) + ((uint64_t)op[11] * 7445765303064700292UL)) * 14);
	tmp_q[2] = ((uint64_t)op[0] * 7445765303064700292UL) + ((uint64_t)op[1] * 1860807927596649322UL) + ((uint64_t)op[2] * 18415531108464612935UL) + ((((uint64_t)op[3] * 7349181473276032313UL) + ((uint64_t)op[4] * 7868938781831905001UL) + ((uint64_t)op[5] * 14206693913815609646UL) + ((uint64_t)op[6] * 13823895817601966885UL) + ((uint64_t)op[7] * 10182052057624725370UL) + ((uint64_t)op[8] * 18031788028952716542UL) + ((uint64_t)op[9] * 13531614810474583979UL) + ((uint64_t)op[10] * 8664828748238991555UL) + ((uint64_t)op[11] * 14995837233037358424UL)) * 14);
	tmp_q[3] = ((uint64_t)op[0] * 14995837233037358424UL) + ((uint64_t)op[1] * 7445765303064700292UL) + ((uint64_t)op[2] * 1860807927596649322UL) + ((uint64_t)op[3] * 18415531108464612935UL) + ((((uint64_t)op[4] * 7349181473276032313UL) + ((uint64_t)op[5] * 7868938781831905001UL) + ((uint64_t)op[6] * 14206693913815609646UL) + ((uint64_t)op[7] * 13823895817601966885UL) + ((uint64_t)op[8] * 10182052057624725370UL) + ((uint64_t)op[9] * 18031788028952716542UL) + ((uint64_t)op[10] * 13531614810474583979UL) + ((uint64_t)op[11] * 8664828748238991555UL)) * 14);
	tmp_q[4] = ((uint64_t)op[0] * 8664828748238991555UL) + ((uint64_t)op[1] * 14995837233037358424UL) + ((uint64_t)op[2] * 7445765303064700292UL) + ((uint64_t)op[3] * 1860807927596649322UL) + ((uint64_t)op[4] * 18415531108464612935UL) + ((((uint64_t)op[5] * 7349181473276032313UL) + ((uint64_t)op[6] * 7868938781831905001UL) + ((uint64_t)op[7] * 14206693913815609646UL) + ((uint64_t)op[8] * 13823895817601966885UL) + ((uint64_t)op[9] * 10182052057624725370UL) + ((uint64_t)op[10] * 18031788028952716542UL) + ((uint64_t)op[11] * 13531614810474583979UL)) * 14);
	tmp_q[5] = ((uint64_t)op[0] * 13531614810474583979UL) + ((uint64_t)op[1] * 8664828748238991555UL) + ((uint64_t)op[2] * 14995837233037358424UL) + ((uint64_t)op[3] * 7445765303064700292UL) + ((uint64_t)op[4] * 1860807927596649322UL) + ((uint64_t)op[5] * 18415531108464612935UL) + ((((uint64_t)op[6] * 7349181473276032313UL) + ((uint64_t)op[7] * 7868938781831905001UL) + ((uint64_t)op[8] * 14206693913815609646UL) + ((uint64_t)op[9] * 13823895817601966885UL) + ((uint64_t)op[10] * 10182052057624725370UL) + ((uint64_t)op[11] * 18031788028952716542UL)) * 14);
	tmp_q[6] = ((uint64_t)op[0] * 18031788028952716542UL) + ((uint64_t)op[1] * 13531614810474583979UL) + ((uint64_t)op[2] * 8664828748238991555UL) + ((uint64_t)op[3] * 14995837233037358424UL) + ((uint64_t)op[4] * 7445765303064700292UL) + ((uint64_t)op[5] * 1860807927596649322UL) + ((uint64_t)op[6] * 18415531108464612935UL) + ((((uint64_t)op[7] * 7349181473276032313UL) + ((uint64_t)op[8] * 7868938781831905001UL) + ((uint64_t)op[9] * 14206693913815609646UL) + ((uint64_t)op[10] * 13823895817601966885UL) + ((uint64_t)op[11] * 10182052057624725370UL)) * 14);
	tmp_q[7] = ((uint64_t)op[0] * 10182052057624725370UL) + ((uint64_t)op[1] * 18031788028952716542UL) + ((uint64_t)op[2] * 13531614810474583979UL) + ((uint64_t)op[3] * 8664828748238991555UL) + ((uint64_t)op[4] * 14995837233037358424UL) + ((uint64_t)op[5] * 7445765303064700292UL) + ((uint64_t)op[6] * 1860807927596649322UL) + ((uint64_t)op[7] * 18415531108464612935UL) + ((((uint64_t)op[8] * 7349181473276032313UL) + ((uint64_t)op[9] * 7868938781831905001UL) + ((uint64_t)op[10] * 14206693913815609646UL) + ((uint64_t)op[11] * 13823895817601966885UL)) * 14);
	tmp_q[8] = ((uint64_t)op[0] * 13823895817601966885UL) + ((uint64_t)op[1] * 10182052057624725370UL) + ((uint64_t)op[2] * 18031788028952716542UL) + ((uint64_t)op[3] * 13531614810474583979UL) + ((uint64_t)op[4] * 8664828748238991555UL) + ((uint64_t)op[5] * 14995837233037358424UL) + ((uint64_t)op[6] * 7445765303064700292UL) + ((uint64_t)op[7] * 1860807927596649322UL) + ((uint64_t)op[8] * 18415531108464612935UL) + ((((uint64_t)op[9] * 7349181473276032313UL) + ((uint64_t)op[10] * 7868938781831905001UL) + ((uint64_t)op[11] * 14206693913815609646UL)) * 14);
	tmp_q[9] = ((uint64_t)op[0] * 14206693913815609646UL) + ((uint64_t)op[1] * 13823895817601966885UL) + ((uint64_t)op[2] * 10182052057624725370UL) + ((uint64_t)op[3] * 18031788028952716542UL) + ((uint64_t)op[4] * 13531614810474583979UL) + ((uint64_t)op[5] * 8664828748238991555UL) + ((uint64_t)op[6] * 14995837233037358424UL) + ((uint64_t)op[7] * 7445765303064700292UL) + ((uint64_t)op[8] * 1860807927596649322UL) + ((uint64_t)op[9] * 18415531108464612935UL) + ((((uint64_t)op[10] * 7349181473276032313UL) + ((uint64_t)op[11] * 7868938781831905001UL)) * 14);
	tmp_q[10] = ((uint64_t)op[0] * 7868938781831905001UL) + ((uint64_t)op[1] * 14206693913815609646UL) + ((uint64_t)op[2] * 13823895817601966885UL) + ((uint64_t)op[3] * 10182052057624725370UL) + ((uint64_t)op[4] * 18031788028952716542UL) + ((uint64_t)op[5] * 13531614810474583979UL) + ((uint64_t)op[6] * 8664828748238991555UL) + ((uint64_t)op[7] * 14995837233037358424UL) + ((uint64_t)op[8] * 7445765303064700292UL) + ((uint64_t)op[9] * 1860807927596649322UL) + ((uint64_t)op[10] * 18415531108464612935UL) + ((uint64_t)op[11] * 10654820257316694302UL);
	tmp_q[11] = ((uint64_t)op[0] * 7349181473276032313UL) + ((uint64_t)op[1] * 7868938781831905001UL) + ((uint64_t)op[2] * 14206693913815609646UL) + ((uint64_t)op[3] * 13823895817601966885UL) + ((uint64_t)op[4] * 10182052057624725370UL) + ((uint64_t)op[5] * 18031788028952716542UL) + ((uint64_t)op[6] * 13531614810474583979UL) + ((uint64_t)op[7] * 8664828748238991555UL) + ((uint64_t)op[8] * 14995837233037358424UL) + ((uint64_t)op[9] * 7445765303064700292UL) + ((uint64_t)op[10] * 1860807927596649322UL) + ((uint64_t)op[11] * 18415531108464612935UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3981685863613L) + ((-((int128)tmp_q[1] * 748909423533L) - ((int128)tmp_q[2] * 2044008684108L) - ((int128)tmp_q[3] * 6615804962984L) + ((int128)tmp_q[4] * 6485404615496L) - ((int128)tmp_q[5] * 3270825733120L) - ((int128)tmp_q[6] * 3177140187402L) - ((int128)tmp_q[7] * 339088445809L) + ((int128)tmp_q[8] * 2045806498689L) - ((int128)tmp_q[9] * 1592755831710L) + ((int128)tmp_q[10] * 2557198974550L) - ((int128)tmp_q[11] * 1806717757764L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 1806717757764L) - ((int128)tmp_q[1] * 3981685863613L) + ((-((int128)tmp_q[2] * 748909423533L) - ((int128)tmp_q[3] * 2044008684108L) - ((int128)tmp_q[4] * 6615804962984L) + ((int128)tmp_q[5] * 6485404615496L) - ((int128)tmp_q[6] * 3270825733120L) - ((int128)tmp_q[7] * 3177140187402L) - ((int128)tmp_q[8] * 339088445809L) + ((int128)tmp_q[9] * 2045806498689L) - ((int128)tmp_q[10] * 1592755831710L) + ((int128)tmp_q[11] * 2557198974550L)) * 14);
	tmp_zero[2] = ((int128)tmp_q[0] * 2557198974550L) - ((int128)tmp_q[1] * 1806717757764L) - ((int128)tmp_q[2] * 3981685863613L) + ((-((int128)tmp_q[3] * 748909423533L) - ((int128)tmp_q[4] * 2044008684108L) - ((int128)tmp_q[5] * 6615804962984L) + ((int128)tmp_q[6] * 6485404615496L) - ((int128)tmp_q[7] * 3270825733120L) - ((int128)tmp_q[8] * 3177140187402L) - ((int128)tmp_q[9] * 339088445809L) + ((int128)tmp_q[10] * 2045806498689L) - ((int128)tmp_q[11] * 1592755831710L)) * 14);
	tmp_zero[3] = -((int128)tmp_q[0] * 1592755831710L) + ((int128)tmp_q[1] * 2557198974550L) - ((int128)tmp_q[2] * 1806717757764L) - ((int128)tmp_q[3] * 3981685863613L) + ((-((int128)tmp_q[4] * 748909423533L) - ((int128)tmp_q[5] * 2044008684108L) - ((int128)tmp_q[6] * 6615804962984L) + ((int128)tmp_q[7] * 6485404615496L) - ((int128)tmp_q[8] * 3270825733120L) - ((int128)tmp_q[9] * 3177140187402L) - ((int128)tmp_q[10] * 339088445809L) + ((int128)tmp_q[11] * 2045806498689L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 2045806498689L) - ((int128)tmp_q[1] * 1592755831710L) + ((int128)tmp_q[2] * 2557198974550L) - ((int128)tmp_q[3] * 1806717757764L) - ((int128)tmp_q[4] * 3981685863613L) + ((-((int128)tmp_q[5] * 748909423533L) - ((int128)tmp_q[6] * 2044008684108L) - ((int128)tmp_q[7] * 6615804962984L) + ((int128)tmp_q[8] * 6485404615496L) - ((int128)tmp_q[9] * 3270825733120L) - ((int128)tmp_q[10] * 3177140187402L) - ((int128)tmp_q[11] * 339088445809L)) * 14);
	tmp_zero[5] = -((int128)tmp_q[0] * 339088445809L) + ((int128)tmp_q[1] * 2045806498689L) - ((int128)tmp_q[2] * 1592755831710L) + ((int128)tmp_q[3] * 2557198974550L) - ((int128)tmp_q[4] * 1806717757764L) - ((int128)tmp_q[5] * 3981685863613L) + ((-((int128)tmp_q[6] * 748909423533L) - ((int128)tmp_q[7] * 2044008684108L) - ((int128)tmp_q[8] * 6615804962984L) + ((int128)tmp_q[9] * 6485404615496L) - ((int128)tmp_q[10] * 3270825733120L) - ((int128)tmp_q[11] * 3177140187402L)) * 14);
	tmp_zero[6] = -((int128)tmp_q[0] * 3177140187402L) - ((int128)tmp_q[1] * 339088445809L) + ((int128)tmp_q[2] * 2045806498689L) - ((int128)tmp_q[3] * 1592755831710L) + ((int128)tmp_q[4] * 2557198974550L) - ((int128)tmp_q[5] * 1806717757764L) - ((int128)tmp_q[6] * 3981685863613L) + ((-((int128)tmp_q[7] * 748909423533L) - ((int128)tmp_q[8] * 2044008684108L) - ((int128)tmp_q[9] * 6615804962984L) + ((int128)tmp_q[10] * 6485404615496L) - ((int128)tmp_q[11] * 3270825733120L)) * 14);
	tmp_zero[7] = -((int128)tmp_q[0] * 3270825733120L) - ((int128)tmp_q[1] * 3177140187402L) - ((int128)tmp_q[2] * 339088445809L) + ((int128)tmp_q[3] * 2045806498689L) - ((int128)tmp_q[4] * 1592755831710L) + ((int128)tmp_q[5] * 2557198974550L) - ((int128)tmp_q[6] * 1806717757764L) - ((int128)tmp_q[7] * 3981685863613L) + ((-((int128)tmp_q[8] * 748909423533L) - ((int128)tmp_q[9] * 2044008684108L) - ((int128)tmp_q[10] * 6615804962984L) + ((int128)tmp_q[11] * 6485404615496L)) * 14);
	tmp_zero[8] = ((int128)tmp_q[0] * 6485404615496L) - ((int128)tmp_q[1] * 3270825733120L) - ((int128)tmp_q[2] * 3177140187402L) - ((int128)tmp_q[3] * 339088445809L) + ((int128)tmp_q[4] * 2045806498689L) - ((int128)tmp_q[5] * 1592755831710L) + ((int128)tmp_q[6] * 2557198974550L) - ((int128)tmp_q[7] * 1806717757764L) - ((int128)tmp_q[8] * 3981685863613L) + ((-((int128)tmp_q[9] * 748909423533L) - ((int128)tmp_q[10] * 2044008684108L) - ((int128)tmp_q[11] * 6615804962984L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 6615804962984L) + ((int128)tmp_q[1] * 6485404615496L) - ((int128)tmp_q[2] * 3270825733120L) - ((int128)tmp_q[3] * 3177140187402L) - ((int128)tmp_q[4] * 339088445809L) + ((int128)tmp_q[5] * 2045806498689L) - ((int128)tmp_q[6] * 1592755831710L) + ((int128)tmp_q[7] * 2557198974550L) - ((int128)tmp_q[8] * 1806717757764L) - ((int128)tmp_q[9] * 3981685863613L) + ((-((int128)tmp_q[10] * 748909423533L) - ((int128)tmp_q[11] * 2044008684108L)) * 14);
	tmp_zero[10] = -((int128)tmp_q[0] * 2044008684108L) - ((int128)tmp_q[1] * 6615804962984L) + ((int128)tmp_q[2] * 6485404615496L) - ((int128)tmp_q[3] * 3270825733120L) - ((int128)tmp_q[4] * 3177140187402L) - ((int128)tmp_q[5] * 339088445809L) + ((int128)tmp_q[6] * 2045806498689L) - ((int128)tmp_q[7] * 1592755831710L) + ((int128)tmp_q[8] * 2557198974550L) - ((int128)tmp_q[9] * 1806717757764L) - ((int128)tmp_q[10] * 3981685863613L) - ((int128)tmp_q[11] * 10484731929462L);
	tmp_zero[11] = -((int128)tmp_q[0] * 748909423533L) - ((int128)tmp_q[1] * 2044008684108L) - ((int128)tmp_q[2] * 6615804962984L) + ((int128)tmp_q[3] * 6485404615496L) - ((int128)tmp_q[4] * 3270825733120L) - ((int128)tmp_q[5] * 3177140187402L) - ((int128)tmp_q[6] * 339088445809L) + ((int128)tmp_q[7] * 2045806498689L) - ((int128)tmp_q[8] * 1592755831710L) + ((int128)tmp_q[9] * 2557198974550L) - ((int128)tmp_q[10] * 1806717757764L) - ((int128)tmp_q[11] * 3981685863613L);

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
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

