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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12338405492130622081UL) + ((((uint64_t)op[1] * 477289389443688398UL) + ((uint64_t)op[2] * 9091045901372149943UL) + ((uint64_t)op[3] * 9045812057505338659UL) + ((uint64_t)op[4] * 17797899346729359247UL) + ((uint64_t)op[5] * 2169095786611653048UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 2169095786611653048UL) + ((uint64_t)op[1] * 12338405492130622081UL) + ((((uint64_t)op[2] * 477289389443688398UL) + ((uint64_t)op[3] * 9091045901372149943UL) + ((uint64_t)op[4] * 9045812057505338659UL) + ((uint64_t)op[5] * 17797899346729359247UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 17797899346729359247UL) + ((uint64_t)op[1] * 2169095786611653048UL) + ((uint64_t)op[2] * 12338405492130622081UL) + ((((uint64_t)op[3] * 477289389443688398UL) + ((uint64_t)op[4] * 9091045901372149943UL) + ((uint64_t)op[5] * 9045812057505338659UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 9045812057505338659UL) + ((uint64_t)op[1] * 17797899346729359247UL) + ((uint64_t)op[2] * 2169095786611653048UL) + ((uint64_t)op[3] * 12338405492130622081UL) + ((((uint64_t)op[4] * 477289389443688398UL) + ((uint64_t)op[5] * 9091045901372149943UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 9091045901372149943UL) + ((uint64_t)op[1] * 9045812057505338659UL) + ((uint64_t)op[2] * 17797899346729359247UL) + ((uint64_t)op[3] * 2169095786611653048UL) + ((uint64_t)op[4] * 12338405492130622081UL) + ((uint64_t)op[5] * 15583007737047421228UL);
	tmp_q[5] = ((uint64_t)op[0] * 477289389443688398UL) + ((uint64_t)op[1] * 9091045901372149943UL) + ((uint64_t)op[2] * 9045812057505338659UL) + ((uint64_t)op[3] * 17797899346729359247UL) + ((uint64_t)op[4] * 2169095786611653048UL) + ((uint64_t)op[5] * 12338405492130622081UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 140102732905L) - ((((int128)tmp_q[1] * 1006480866780L) + ((int128)tmp_q[2] * 2437763280696L) - ((int128)tmp_q[3] * 2778135051023L) + ((int128)tmp_q[4] * 4350279776019L) - ((int128)tmp_q[5] * 3676477289434L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 3676477289434L) - ((int128)tmp_q[1] * 140102732905L) - ((((int128)tmp_q[2] * 1006480866780L) + ((int128)tmp_q[3] * 2437763280696L) - ((int128)tmp_q[4] * 2778135051023L) + ((int128)tmp_q[5] * 4350279776019L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 4350279776019L) - ((int128)tmp_q[1] * 3676477289434L) - ((int128)tmp_q[2] * 140102732905L) - ((((int128)tmp_q[3] * 1006480866780L) + ((int128)tmp_q[4] * 2437763280696L) - ((int128)tmp_q[5] * 2778135051023L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 2778135051023L) + ((int128)tmp_q[1] * 4350279776019L) - ((int128)tmp_q[2] * 3676477289434L) - ((int128)tmp_q[3] * 140102732905L) - ((((int128)tmp_q[4] * 1006480866780L) + ((int128)tmp_q[5] * 2437763280696L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 2437763280696L) - ((int128)tmp_q[1] * 2778135051023L) + ((int128)tmp_q[2] * 4350279776019L) - ((int128)tmp_q[3] * 3676477289434L) - ((int128)tmp_q[4] * 140102732905L) - ((int128)tmp_q[5] * 6038885200680L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1006480866780L) + ((int128)tmp_q[1] * 2437763280696L) - ((int128)tmp_q[2] * 2778135051023L) + ((int128)tmp_q[3] * 4350279776019L) - ((int128)tmp_q[4] * 3676477289434L) - ((int128)tmp_q[5] * 140102732905L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

