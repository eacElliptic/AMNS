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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9282402602135281887UL) + ((((uint64_t)op[1] * 17776312502898389807UL) + ((uint64_t)op[2] * 7290059319876608324UL) + ((uint64_t)op[3] * 12014736450102806765UL) + ((uint64_t)op[4] * 3399037164017689250UL) + ((uint64_t)op[5] * 1682956648982724691UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 1682956648982724691UL) + ((uint64_t)op[1] * 9282402602135281887UL) + ((((uint64_t)op[2] * 17776312502898389807UL) + ((uint64_t)op[3] * 7290059319876608324UL) + ((uint64_t)op[4] * 12014736450102806765UL) + ((uint64_t)op[5] * 3399037164017689250UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 3399037164017689250UL) + ((uint64_t)op[1] * 1682956648982724691UL) + ((uint64_t)op[2] * 9282402602135281887UL) + ((((uint64_t)op[3] * 17776312502898389807UL) + ((uint64_t)op[4] * 7290059319876608324UL) + ((uint64_t)op[5] * 12014736450102806765UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 12014736450102806765UL) + ((uint64_t)op[1] * 3399037164017689250UL) + ((uint64_t)op[2] * 1682956648982724691UL) + ((uint64_t)op[3] * 9282402602135281887UL) + ((((uint64_t)op[4] * 17776312502898389807UL) + ((uint64_t)op[5] * 7290059319876608324UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 7290059319876608324UL) + ((uint64_t)op[1] * 12014736450102806765UL) + ((uint64_t)op[2] * 3399037164017689250UL) + ((uint64_t)op[3] * 1682956648982724691UL) + ((uint64_t)op[4] * 9282402602135281887UL) + ((uint64_t)op[5] * 1340863141622323618UL);
	tmp_q[5] = ((uint64_t)op[0] * 17776312502898389807UL) + ((uint64_t)op[1] * 7290059319876608324UL) + ((uint64_t)op[2] * 12014736450102806765UL) + ((uint64_t)op[3] * 3399037164017689250UL) + ((uint64_t)op[4] * 1682956648982724691UL) + ((uint64_t)op[5] * 9282402602135281887UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 386485634293L) - ((((int128)tmp_q[1] * 964620256965L) - ((int128)tmp_q[2] * 952749757161L) + ((int128)tmp_q[3] * 3659116788996L) - ((int128)tmp_q[4] * 616702038507L) - ((int128)tmp_q[5] * 5472805248305L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 5472805248305L) + ((int128)tmp_q[1] * 386485634293L) - ((((int128)tmp_q[2] * 964620256965L) - ((int128)tmp_q[3] * 952749757161L) + ((int128)tmp_q[4] * 3659116788996L) - ((int128)tmp_q[5] * 616702038507L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 616702038507L) - ((int128)tmp_q[1] * 5472805248305L) + ((int128)tmp_q[2] * 386485634293L) - ((((int128)tmp_q[3] * 964620256965L) - ((int128)tmp_q[4] * 952749757161L) + ((int128)tmp_q[5] * 3659116788996L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 3659116788996L) - ((int128)tmp_q[1] * 616702038507L) - ((int128)tmp_q[2] * 5472805248305L) + ((int128)tmp_q[3] * 386485634293L) - ((((int128)tmp_q[4] * 964620256965L) - ((int128)tmp_q[5] * 952749757161L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 952749757161L) + ((int128)tmp_q[1] * 3659116788996L) - ((int128)tmp_q[2] * 616702038507L) - ((int128)tmp_q[3] * 5472805248305L) + ((int128)tmp_q[4] * 386485634293L) - ((int128)tmp_q[5] * 1929240513930L);
	tmp_zero[5] = ((int128)tmp_q[0] * 964620256965L) - ((int128)tmp_q[1] * 952749757161L) + ((int128)tmp_q[2] * 3659116788996L) - ((int128)tmp_q[3] * 616702038507L) - ((int128)tmp_q[4] * 5472805248305L) + ((int128)tmp_q[5] * 386485634293L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

