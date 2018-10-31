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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14543758820756763305UL) + ((((uint64_t)op[1] * 8370661783949653883UL) + ((uint64_t)op[2] * 13013691345033838675UL) + ((uint64_t)op[3] * 12901735865263496257UL) + ((uint64_t)op[4] * 10195808684652200459UL) + ((uint64_t)op[5] * 190663711029241994UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 190663711029241994UL) + ((uint64_t)op[1] * 14543758820756763305UL) + ((((uint64_t)op[2] * 8370661783949653883UL) + ((uint64_t)op[3] * 13013691345033838675UL) + ((uint64_t)op[4] * 12901735865263496257UL) + ((uint64_t)op[5] * 10195808684652200459UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 10195808684652200459UL) + ((uint64_t)op[1] * 190663711029241994UL) + ((uint64_t)op[2] * 14543758820756763305UL) + ((((uint64_t)op[3] * 8370661783949653883UL) + ((uint64_t)op[4] * 13013691345033838675UL) + ((uint64_t)op[5] * 12901735865263496257UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 12901735865263496257UL) + ((uint64_t)op[1] * 10195808684652200459UL) + ((uint64_t)op[2] * 190663711029241994UL) + ((uint64_t)op[3] * 14543758820756763305UL) + ((((uint64_t)op[4] * 8370661783949653883UL) + ((uint64_t)op[5] * 13013691345033838675UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 13013691345033838675UL) + ((uint64_t)op[1] * 12901735865263496257UL) + ((uint64_t)op[2] * 10195808684652200459UL) + ((uint64_t)op[3] * 190663711029241994UL) + ((uint64_t)op[4] * 14543758820756763305UL) + ((uint64_t)op[5] * 10076082289759897733UL);
	tmp_q[5] = ((uint64_t)op[0] * 8370661783949653883UL) + ((uint64_t)op[1] * 13013691345033838675UL) + ((uint64_t)op[2] * 12901735865263496257UL) + ((uint64_t)op[3] * 10195808684652200459UL) + ((uint64_t)op[4] * 190663711029241994UL) + ((uint64_t)op[5] * 14543758820756763305UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5886194106221L) - (((int128)tmp_q[1] * 41180980500154L) + ((int128)tmp_q[2] * 61866182682047L) - ((int128)tmp_q[3] * 113990481779095L) + ((int128)tmp_q[4] * 55979988575827L) - ((int128)tmp_q[5] * 155171462279245L));
	tmp_zero[1] = -((int128)tmp_q[0] * 155171462279245L) - ((int128)tmp_q[1] * 5886194106221L) - (((int128)tmp_q[2] * 41180980500154L) + ((int128)tmp_q[3] * 61866182682047L) - ((int128)tmp_q[4] * 113990481779095L) + ((int128)tmp_q[5] * 55979988575827L));
	tmp_zero[2] = ((int128)tmp_q[0] * 55979988575827L) - ((int128)tmp_q[1] * 155171462279245L) - ((int128)tmp_q[2] * 5886194106221L) - (((int128)tmp_q[3] * 41180980500154L) + ((int128)tmp_q[4] * 61866182682047L) - ((int128)tmp_q[5] * 113990481779095L));
	tmp_zero[3] = -((int128)tmp_q[0] * 113990481779095L) + ((int128)tmp_q[1] * 55979988575827L) - ((int128)tmp_q[2] * 155171462279245L) - ((int128)tmp_q[3] * 5886194106221L) - (((int128)tmp_q[4] * 41180980500154L) + ((int128)tmp_q[5] * 61866182682047L));
	tmp_zero[4] = ((int128)tmp_q[0] * 61866182682047L) - ((int128)tmp_q[1] * 113990481779095L) + ((int128)tmp_q[2] * 55979988575827L) - ((int128)tmp_q[3] * 155171462279245L) - ((int128)tmp_q[4] * 5886194106221L) - ((int128)tmp_q[5] * 41180980500154L);
	tmp_zero[5] = ((int128)tmp_q[0] * 41180980500154L) + ((int128)tmp_q[1] * 61866182682047L) - ((int128)tmp_q[2] * 113990481779095L) + ((int128)tmp_q[3] * 55979988575827L) - ((int128)tmp_q[4] * 155171462279245L) - ((int128)tmp_q[5] * 5886194106221L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

