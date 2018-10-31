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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13554379172486715533UL) + ((((uint64_t)op[1] * 13863390777763217745UL) + ((uint64_t)op[2] * 1045480459873565985UL) + ((uint64_t)op[3] * 16841669715671922322UL) + ((uint64_t)op[4] * 6904954291449178243UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 6904954291449178243UL) + ((uint64_t)op[1] * 13554379172486715533UL) + ((((uint64_t)op[2] * 13863390777763217745UL) + ((uint64_t)op[3] * 1045480459873565985UL) + ((uint64_t)op[4] * 16841669715671922322UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 16841669715671922322UL) + ((uint64_t)op[1] * 6904954291449178243UL) + ((uint64_t)op[2] * 13554379172486715533UL) + ((((uint64_t)op[3] * 13863390777763217745UL) + ((uint64_t)op[4] * 1045480459873565985UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 1045480459873565985UL) + ((uint64_t)op[1] * 16841669715671922322UL) + ((uint64_t)op[2] * 6904954291449178243UL) + ((uint64_t)op[3] * 13554379172486715533UL) + ((uint64_t)op[4] * 9166706591892667742UL);
	tmp_q[4] = ((uint64_t)op[0] * 13863390777763217745UL) + ((uint64_t)op[1] * 1045480459873565985UL) + ((uint64_t)op[2] * 16841669715671922322UL) + ((uint64_t)op[3] * 6904954291449178243UL) + ((uint64_t)op[4] * 13554379172486715533UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1350934316191595L) - ((-((int128)tmp_q[1] * 125347164713036L) - ((int128)tmp_q[2] * 1071675124402294L) - ((int128)tmp_q[3] * 1462942292505637L) - ((int128)tmp_q[4] * 306559457019683L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 306559457019683L) + ((int128)tmp_q[1] * 1350934316191595L) - ((-((int128)tmp_q[2] * 125347164713036L) - ((int128)tmp_q[3] * 1071675124402294L) - ((int128)tmp_q[4] * 1462942292505637L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 1462942292505637L) - ((int128)tmp_q[1] * 306559457019683L) + ((int128)tmp_q[2] * 1350934316191595L) - ((-((int128)tmp_q[3] * 125347164713036L) - ((int128)tmp_q[4] * 1071675124402294L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 1071675124402294L) - ((int128)tmp_q[1] * 1462942292505637L) - ((int128)tmp_q[2] * 306559457019683L) + ((int128)tmp_q[3] * 1350934316191595L) + ((int128)tmp_q[4] * 250694329426072L);
	tmp_zero[4] = -((int128)tmp_q[0] * 125347164713036L) - ((int128)tmp_q[1] * 1071675124402294L) - ((int128)tmp_q[2] * 1462942292505637L) - ((int128)tmp_q[3] * 306559457019683L) + ((int128)tmp_q[4] * 1350934316191595L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

