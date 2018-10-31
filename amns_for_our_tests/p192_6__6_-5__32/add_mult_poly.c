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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11615627686365777392UL) + ((((uint64_t)op[1] * 6586389065406720550UL) + ((uint64_t)op[2] * 13657089869526210791UL) + ((uint64_t)op[3] * 14556818883108040482UL) + ((uint64_t)op[4] * 17709708585967350907UL) + ((uint64_t)op[5] * 12142101963764755313UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 12142101963764755313UL) + ((uint64_t)op[1] * 11615627686365777392UL) + ((((uint64_t)op[2] * 6586389065406720550UL) + ((uint64_t)op[3] * 13657089869526210791UL) + ((uint64_t)op[4] * 14556818883108040482UL) + ((uint64_t)op[5] * 17709708585967350907UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 17709708585967350907UL) + ((uint64_t)op[1] * 12142101963764755313UL) + ((uint64_t)op[2] * 11615627686365777392UL) + ((((uint64_t)op[3] * 6586389065406720550UL) + ((uint64_t)op[4] * 13657089869526210791UL) + ((uint64_t)op[5] * 14556818883108040482UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 14556818883108040482UL) + ((uint64_t)op[1] * 17709708585967350907UL) + ((uint64_t)op[2] * 12142101963764755313UL) + ((uint64_t)op[3] * 11615627686365777392UL) + ((((uint64_t)op[4] * 6586389065406720550UL) + ((uint64_t)op[5] * 13657089869526210791UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 13657089869526210791UL) + ((uint64_t)op[1] * 14556818883108040482UL) + ((uint64_t)op[2] * 17709708585967350907UL) + ((uint64_t)op[3] * 12142101963764755313UL) + ((uint64_t)op[4] * 11615627686365777392UL) + ((uint64_t)op[5] * 3961542820385500482UL);
	tmp_q[5] = ((uint64_t)op[0] * 6586389065406720550UL) + ((uint64_t)op[1] * 13657089869526210791UL) + ((uint64_t)op[2] * 14556818883108040482UL) + ((uint64_t)op[3] * 17709708585967350907UL) + ((uint64_t)op[4] * 12142101963764755313UL) + ((uint64_t)op[5] * 11615627686365777392UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 202720139L) - ((-((int128)tmp_q[1] * 1408662120L) + ((int128)tmp_q[2] * 1089999141L) + ((int128)tmp_q[3] * 2238135323L) - ((int128)tmp_q[4] * 4046568L) + ((int128)tmp_q[5] * 435273376L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 435273376L) + ((int128)tmp_q[1] * 202720139L) - ((-((int128)tmp_q[2] * 1408662120L) + ((int128)tmp_q[3] * 1089999141L) + ((int128)tmp_q[4] * 2238135323L) - ((int128)tmp_q[5] * 4046568L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 4046568L) + ((int128)tmp_q[1] * 435273376L) + ((int128)tmp_q[2] * 202720139L) - ((-((int128)tmp_q[3] * 1408662120L) + ((int128)tmp_q[4] * 1089999141L) + ((int128)tmp_q[5] * 2238135323L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 2238135323L) - ((int128)tmp_q[1] * 4046568L) + ((int128)tmp_q[2] * 435273376L) + ((int128)tmp_q[3] * 202720139L) - ((-((int128)tmp_q[4] * 1408662120L) + ((int128)tmp_q[5] * 1089999141L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1089999141L) + ((int128)tmp_q[1] * 2238135323L) - ((int128)tmp_q[2] * 4046568L) + ((int128)tmp_q[3] * 435273376L) + ((int128)tmp_q[4] * 202720139L) + ((int128)tmp_q[5] * 7043310600L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1408662120L) + ((int128)tmp_q[1] * 1089999141L) + ((int128)tmp_q[2] * 2238135323L) - ((int128)tmp_q[3] * 4046568L) + ((int128)tmp_q[4] * 435273376L) + ((int128)tmp_q[5] * 202720139L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

