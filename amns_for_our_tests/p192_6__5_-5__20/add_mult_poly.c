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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13781941190975465115UL) + ((((uint64_t)op[1] * 12494358785435343403UL) + ((uint64_t)op[2] * 15284305611394758266UL) + ((uint64_t)op[3] * 6879669878514351069UL) + ((uint64_t)op[4] * 17671257609711609014UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 17671257609711609014UL) + ((uint64_t)op[1] * 13781941190975465115UL) + ((((uint64_t)op[2] * 12494358785435343403UL) + ((uint64_t)op[3] * 15284305611394758266UL) + ((uint64_t)op[4] * 6879669878514351069UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 6879669878514351069UL) + ((uint64_t)op[1] * 17671257609711609014UL) + ((uint64_t)op[2] * 13781941190975465115UL) + ((((uint64_t)op[3] * 12494358785435343403UL) + ((uint64_t)op[4] * 15284305611394758266UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 15284305611394758266UL) + ((uint64_t)op[1] * 6879669878514351069UL) + ((uint64_t)op[2] * 17671257609711609014UL) + ((uint64_t)op[3] * 13781941190975465115UL) + ((uint64_t)op[4] * 11315182367661489449UL);
	tmp_q[4] = ((uint64_t)op[0] * 12494358785435343403UL) + ((uint64_t)op[1] * 15284305611394758266UL) + ((uint64_t)op[2] * 6879669878514351069UL) + ((uint64_t)op[3] * 17671257609711609014UL) + ((uint64_t)op[4] * 13781941190975465115UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 38013252294L) - ((((int128)tmp_q[1] * 140135843177L) + ((int128)tmp_q[2] * 101148663593L) + ((int128)tmp_q[3] * 118824269993L) + ((int128)tmp_q[4] * 110703237562L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 110703237562L) - ((int128)tmp_q[1] * 38013252294L) - ((((int128)tmp_q[2] * 140135843177L) + ((int128)tmp_q[3] * 101148663593L) + ((int128)tmp_q[4] * 118824269993L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 118824269993L) + ((int128)tmp_q[1] * 110703237562L) - ((int128)tmp_q[2] * 38013252294L) - ((((int128)tmp_q[3] * 140135843177L) + ((int128)tmp_q[4] * 101148663593L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 101148663593L) + ((int128)tmp_q[1] * 118824269993L) + ((int128)tmp_q[2] * 110703237562L) - ((int128)tmp_q[3] * 38013252294L) - ((int128)tmp_q[4] * 700679215885L);
	tmp_zero[4] = ((int128)tmp_q[0] * 140135843177L) + ((int128)tmp_q[1] * 101148663593L) + ((int128)tmp_q[2] * 118824269993L) + ((int128)tmp_q[3] * 110703237562L) - ((int128)tmp_q[4] * 38013252294L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

