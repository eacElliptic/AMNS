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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8878384086413837487UL) + ((((uint64_t)op[1] * 12277110042431094825UL) + ((uint64_t)op[2] * 7924220360739379701UL) + ((uint64_t)op[3] * 6869868569321039061UL) + ((uint64_t)op[4] * 9725446241631998340UL) + ((uint64_t)op[5] * 7309768032706497762UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 7309768032706497762UL) + ((uint64_t)op[1] * 8878384086413837487UL) + ((((uint64_t)op[2] * 12277110042431094825UL) + ((uint64_t)op[3] * 7924220360739379701UL) + ((uint64_t)op[4] * 6869868569321039061UL) + ((uint64_t)op[5] * 9725446241631998340UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 9725446241631998340UL) + ((uint64_t)op[1] * 7309768032706497762UL) + ((uint64_t)op[2] * 8878384086413837487UL) + ((((uint64_t)op[3] * 12277110042431094825UL) + ((uint64_t)op[4] * 7924220360739379701UL) + ((uint64_t)op[5] * 6869868569321039061UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 6869868569321039061UL) + ((uint64_t)op[1] * 9725446241631998340UL) + ((uint64_t)op[2] * 7309768032706497762UL) + ((uint64_t)op[3] * 8878384086413837487UL) + ((((uint64_t)op[4] * 12277110042431094825UL) + ((uint64_t)op[5] * 7924220360739379701UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 7924220360739379701UL) + ((uint64_t)op[1] * 6869868569321039061UL) + ((uint64_t)op[2] * 9725446241631998340UL) + ((uint64_t)op[3] * 7309768032706497762UL) + ((uint64_t)op[4] * 8878384086413837487UL) + ((uint64_t)op[5] * 12463584102808551096UL);
	tmp_q[5] = ((uint64_t)op[0] * 12277110042431094825UL) + ((uint64_t)op[1] * 7924220360739379701UL) + ((uint64_t)op[2] * 6869868569321039061UL) + ((uint64_t)op[3] * 9725446241631998340UL) + ((uint64_t)op[4] * 7309768032706497762UL) + ((uint64_t)op[5] * 8878384086413837487UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 30972965777577L) - ((-((int128)tmp_q[1] * 38839150371455L) + ((int128)tmp_q[2] * 14373606639929L) + ((int128)tmp_q[3] * 55408717178181L) + ((int128)tmp_q[4] * 69132734875120L) + ((int128)tmp_q[5] * 37414146181954L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 37414146181954L) + ((int128)tmp_q[1] * 30972965777577L) - ((-((int128)tmp_q[2] * 38839150371455L) + ((int128)tmp_q[3] * 14373606639929L) + ((int128)tmp_q[4] * 55408717178181L) + ((int128)tmp_q[5] * 69132734875120L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 69132734875120L) + ((int128)tmp_q[1] * 37414146181954L) + ((int128)tmp_q[2] * 30972965777577L) - ((-((int128)tmp_q[3] * 38839150371455L) + ((int128)tmp_q[4] * 14373606639929L) + ((int128)tmp_q[5] * 55408717178181L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 55408717178181L) + ((int128)tmp_q[1] * 69132734875120L) + ((int128)tmp_q[2] * 37414146181954L) + ((int128)tmp_q[3] * 30972965777577L) - ((-((int128)tmp_q[4] * 38839150371455L) + ((int128)tmp_q[5] * 14373606639929L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 14373606639929L) + ((int128)tmp_q[1] * 55408717178181L) + ((int128)tmp_q[2] * 69132734875120L) + ((int128)tmp_q[3] * 37414146181954L) + ((int128)tmp_q[4] * 30972965777577L) + ((int128)tmp_q[5] * 310713202971640L);
	tmp_zero[5] = -((int128)tmp_q[0] * 38839150371455L) + ((int128)tmp_q[1] * 14373606639929L) + ((int128)tmp_q[2] * 55408717178181L) + ((int128)tmp_q[3] * 69132734875120L) + ((int128)tmp_q[4] * 37414146181954L) + ((int128)tmp_q[5] * 30972965777577L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

