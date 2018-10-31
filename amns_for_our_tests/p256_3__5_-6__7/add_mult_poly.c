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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2484407108893508567UL) + ((((uint64_t)op[1] * 8103400089933420426UL) + ((uint64_t)op[2] * 7989408392120257708UL) + ((uint64_t)op[3] * 1934375871653383501UL) + ((uint64_t)op[4] * 17454886645046703909UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 17454886645046703909UL) + ((uint64_t)op[1] * 2484407108893508567UL) + ((((uint64_t)op[2] * 8103400089933420426UL) + ((uint64_t)op[3] * 7989408392120257708UL) + ((uint64_t)op[4] * 1934375871653383501UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 1934375871653383501UL) + ((uint64_t)op[1] * 17454886645046703909UL) + ((uint64_t)op[2] * 2484407108893508567UL) + ((((uint64_t)op[3] * 8103400089933420426UL) + ((uint64_t)op[4] * 7989408392120257708UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 7989408392120257708UL) + ((uint64_t)op[1] * 1934375871653383501UL) + ((uint64_t)op[2] * 17454886645046703909UL) + ((uint64_t)op[3] * 2484407108893508567UL) + ((uint64_t)op[4] * 6719831681528132292UL);
	tmp_q[4] = ((uint64_t)op[0] * 8103400089933420426UL) + ((uint64_t)op[1] * 7989408392120257708UL) + ((uint64_t)op[2] * 1934375871653383501UL) + ((uint64_t)op[3] * 17454886645046703909UL) + ((uint64_t)op[4] * 2484407108893508567UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 841281320421299L) - ((((int128)tmp_q[1] * 1056592571912353L) - ((int128)tmp_q[2] * 195621063241061L) + ((int128)tmp_q[3] * 78198353308940L) + ((int128)tmp_q[4] * 817851721443315L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 817851721443315L) - ((int128)tmp_q[1] * 841281320421299L) - ((((int128)tmp_q[2] * 1056592571912353L) - ((int128)tmp_q[3] * 195621063241061L) + ((int128)tmp_q[4] * 78198353308940L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 78198353308940L) + ((int128)tmp_q[1] * 817851721443315L) - ((int128)tmp_q[2] * 841281320421299L) - ((((int128)tmp_q[3] * 1056592571912353L) - ((int128)tmp_q[4] * 195621063241061L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 195621063241061L) + ((int128)tmp_q[1] * 78198353308940L) + ((int128)tmp_q[2] * 817851721443315L) - ((int128)tmp_q[3] * 841281320421299L) - ((int128)tmp_q[4] * 6339555431474118L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1056592571912353L) - ((int128)tmp_q[1] * 195621063241061L) + ((int128)tmp_q[2] * 78198353308940L) + ((int128)tmp_q[3] * 817851721443315L) - ((int128)tmp_q[4] * 841281320421299L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

