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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4596368476408638760UL) + ((((uint64_t)op[1] * 14950397632264957497UL) + ((uint64_t)op[2] * 7567316585608783756UL) + ((uint64_t)op[3] * 1668235830785023579UL) + ((uint64_t)op[4] * 17831004001636385115UL) + ((uint64_t)op[5] * 1666677813883851975UL) + ((uint64_t)op[6] * 7352794297248518981UL) + ((uint64_t)op[7] * 14475035243038172198UL) + ((uint64_t)op[8] * 9438713371523988679UL) + ((uint64_t)op[9] * 17343221395997402425UL) + ((uint64_t)op[10] * 16747214901427369314UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16747214901427369314UL) + ((uint64_t)op[1] * 4596368476408638760UL) + ((((uint64_t)op[2] * 14950397632264957497UL) + ((uint64_t)op[3] * 7567316585608783756UL) + ((uint64_t)op[4] * 1668235830785023579UL) + ((uint64_t)op[5] * 17831004001636385115UL) + ((uint64_t)op[6] * 1666677813883851975UL) + ((uint64_t)op[7] * 7352794297248518981UL) + ((uint64_t)op[8] * 14475035243038172198UL) + ((uint64_t)op[9] * 9438713371523988679UL) + ((uint64_t)op[10] * 17343221395997402425UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 17343221395997402425UL) + ((uint64_t)op[1] * 16747214901427369314UL) + ((uint64_t)op[2] * 4596368476408638760UL) + ((((uint64_t)op[3] * 14950397632264957497UL) + ((uint64_t)op[4] * 7567316585608783756UL) + ((uint64_t)op[5] * 1668235830785023579UL) + ((uint64_t)op[6] * 17831004001636385115UL) + ((uint64_t)op[7] * 1666677813883851975UL) + ((uint64_t)op[8] * 7352794297248518981UL) + ((uint64_t)op[9] * 14475035243038172198UL) + ((uint64_t)op[10] * 9438713371523988679UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 9438713371523988679UL) + ((uint64_t)op[1] * 17343221395997402425UL) + ((uint64_t)op[2] * 16747214901427369314UL) + ((uint64_t)op[3] * 4596368476408638760UL) + ((((uint64_t)op[4] * 14950397632264957497UL) + ((uint64_t)op[5] * 7567316585608783756UL) + ((uint64_t)op[6] * 1668235830785023579UL) + ((uint64_t)op[7] * 17831004001636385115UL) + ((uint64_t)op[8] * 1666677813883851975UL) + ((uint64_t)op[9] * 7352794297248518981UL) + ((uint64_t)op[10] * 14475035243038172198UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 14475035243038172198UL) + ((uint64_t)op[1] * 9438713371523988679UL) + ((uint64_t)op[2] * 17343221395997402425UL) + ((uint64_t)op[3] * 16747214901427369314UL) + ((uint64_t)op[4] * 4596368476408638760UL) + ((((uint64_t)op[5] * 14950397632264957497UL) + ((uint64_t)op[6] * 7567316585608783756UL) + ((uint64_t)op[7] * 1668235830785023579UL) + ((uint64_t)op[8] * 17831004001636385115UL) + ((uint64_t)op[9] * 1666677813883851975UL) + ((uint64_t)op[10] * 7352794297248518981UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 7352794297248518981UL) + ((uint64_t)op[1] * 14475035243038172198UL) + ((uint64_t)op[2] * 9438713371523988679UL) + ((uint64_t)op[3] * 17343221395997402425UL) + ((uint64_t)op[4] * 16747214901427369314UL) + ((uint64_t)op[5] * 4596368476408638760UL) + ((((uint64_t)op[6] * 14950397632264957497UL) + ((uint64_t)op[7] * 7567316585608783756UL) + ((uint64_t)op[8] * 1668235830785023579UL) + ((uint64_t)op[9] * 17831004001636385115UL) + ((uint64_t)op[10] * 1666677813883851975UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 1666677813883851975UL) + ((uint64_t)op[1] * 7352794297248518981UL) + ((uint64_t)op[2] * 14475035243038172198UL) + ((uint64_t)op[3] * 9438713371523988679UL) + ((uint64_t)op[4] * 17343221395997402425UL) + ((uint64_t)op[5] * 16747214901427369314UL) + ((uint64_t)op[6] * 4596368476408638760UL) + ((((uint64_t)op[7] * 14950397632264957497UL) + ((uint64_t)op[8] * 7567316585608783756UL) + ((uint64_t)op[9] * 1668235830785023579UL) + ((uint64_t)op[10] * 17831004001636385115UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 17831004001636385115UL) + ((uint64_t)op[1] * 1666677813883851975UL) + ((uint64_t)op[2] * 7352794297248518981UL) + ((uint64_t)op[3] * 14475035243038172198UL) + ((uint64_t)op[4] * 9438713371523988679UL) + ((uint64_t)op[5] * 17343221395997402425UL) + ((uint64_t)op[6] * 16747214901427369314UL) + ((uint64_t)op[7] * 4596368476408638760UL) + ((((uint64_t)op[8] * 14950397632264957497UL) + ((uint64_t)op[9] * 7567316585608783756UL) + ((uint64_t)op[10] * 1668235830785023579UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 1668235830785023579UL) + ((uint64_t)op[1] * 17831004001636385115UL) + ((uint64_t)op[2] * 1666677813883851975UL) + ((uint64_t)op[3] * 7352794297248518981UL) + ((uint64_t)op[4] * 14475035243038172198UL) + ((uint64_t)op[5] * 9438713371523988679UL) + ((uint64_t)op[6] * 17343221395997402425UL) + ((uint64_t)op[7] * 16747214901427369314UL) + ((uint64_t)op[8] * 4596368476408638760UL) + ((((uint64_t)op[9] * 14950397632264957497UL) + ((uint64_t)op[10] * 7567316585608783756UL)) * 18446744073709551613);
	tmp_q[9] = ((uint64_t)op[0] * 7567316585608783756UL) + ((uint64_t)op[1] * 1668235830785023579UL) + ((uint64_t)op[2] * 17831004001636385115UL) + ((uint64_t)op[3] * 1666677813883851975UL) + ((uint64_t)op[4] * 7352794297248518981UL) + ((uint64_t)op[5] * 14475035243038172198UL) + ((uint64_t)op[6] * 9438713371523988679UL) + ((uint64_t)op[7] * 17343221395997402425UL) + ((uint64_t)op[8] * 16747214901427369314UL) + ((uint64_t)op[9] * 4596368476408638760UL) + ((uint64_t)op[10] * 10489039324333782357UL);
	tmp_q[10] = ((uint64_t)op[0] * 14950397632264957497UL) + ((uint64_t)op[1] * 7567316585608783756UL) + ((uint64_t)op[2] * 1668235830785023579UL) + ((uint64_t)op[3] * 17831004001636385115UL) + ((uint64_t)op[4] * 1666677813883851975UL) + ((uint64_t)op[5] * 7352794297248518981UL) + ((uint64_t)op[6] * 14475035243038172198UL) + ((uint64_t)op[7] * 9438713371523988679UL) + ((uint64_t)op[8] * 17343221395997402425UL) + ((uint64_t)op[9] * 16747214901427369314UL) + ((uint64_t)op[10] * 4596368476408638760UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 52001840031308L) - ((-((int128)tmp_q[1] * 30830197480814L) - ((int128)tmp_q[2] * 85527218674374L) + ((int128)tmp_q[3] * 85362351060129L) - ((int128)tmp_q[4] * 49750569740116L) + ((int128)tmp_q[5] * 31531415086882L) - ((int128)tmp_q[6] * 34107822775284L) - ((int128)tmp_q[7] * 129388969751772L) + ((int128)tmp_q[8] * 104390065713229L) - ((int128)tmp_q[9] * 13504404300340L) - ((int128)tmp_q[10] * 80144319299779L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 80144319299779L) + ((int128)tmp_q[1] * 52001840031308L) - ((-((int128)tmp_q[2] * 30830197480814L) - ((int128)tmp_q[3] * 85527218674374L) + ((int128)tmp_q[4] * 85362351060129L) - ((int128)tmp_q[5] * 49750569740116L) + ((int128)tmp_q[6] * 31531415086882L) - ((int128)tmp_q[7] * 34107822775284L) - ((int128)tmp_q[8] * 129388969751772L) + ((int128)tmp_q[9] * 104390065713229L) - ((int128)tmp_q[10] * 13504404300340L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 13504404300340L) - ((int128)tmp_q[1] * 80144319299779L) + ((int128)tmp_q[2] * 52001840031308L) - ((-((int128)tmp_q[3] * 30830197480814L) - ((int128)tmp_q[4] * 85527218674374L) + ((int128)tmp_q[5] * 85362351060129L) - ((int128)tmp_q[6] * 49750569740116L) + ((int128)tmp_q[7] * 31531415086882L) - ((int128)tmp_q[8] * 34107822775284L) - ((int128)tmp_q[9] * 129388969751772L) + ((int128)tmp_q[10] * 104390065713229L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 104390065713229L) - ((int128)tmp_q[1] * 13504404300340L) - ((int128)tmp_q[2] * 80144319299779L) + ((int128)tmp_q[3] * 52001840031308L) - ((-((int128)tmp_q[4] * 30830197480814L) - ((int128)tmp_q[5] * 85527218674374L) + ((int128)tmp_q[6] * 85362351060129L) - ((int128)tmp_q[7] * 49750569740116L) + ((int128)tmp_q[8] * 31531415086882L) - ((int128)tmp_q[9] * 34107822775284L) - ((int128)tmp_q[10] * 129388969751772L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 129388969751772L) + ((int128)tmp_q[1] * 104390065713229L) - ((int128)tmp_q[2] * 13504404300340L) - ((int128)tmp_q[3] * 80144319299779L) + ((int128)tmp_q[4] * 52001840031308L) - ((-((int128)tmp_q[5] * 30830197480814L) - ((int128)tmp_q[6] * 85527218674374L) + ((int128)tmp_q[7] * 85362351060129L) - ((int128)tmp_q[8] * 49750569740116L) + ((int128)tmp_q[9] * 31531415086882L) - ((int128)tmp_q[10] * 34107822775284L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 34107822775284L) - ((int128)tmp_q[1] * 129388969751772L) + ((int128)tmp_q[2] * 104390065713229L) - ((int128)tmp_q[3] * 13504404300340L) - ((int128)tmp_q[4] * 80144319299779L) + ((int128)tmp_q[5] * 52001840031308L) - ((-((int128)tmp_q[6] * 30830197480814L) - ((int128)tmp_q[7] * 85527218674374L) + ((int128)tmp_q[8] * 85362351060129L) - ((int128)tmp_q[9] * 49750569740116L) + ((int128)tmp_q[10] * 31531415086882L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 31531415086882L) - ((int128)tmp_q[1] * 34107822775284L) - ((int128)tmp_q[2] * 129388969751772L) + ((int128)tmp_q[3] * 104390065713229L) - ((int128)tmp_q[4] * 13504404300340L) - ((int128)tmp_q[5] * 80144319299779L) + ((int128)tmp_q[6] * 52001840031308L) - ((-((int128)tmp_q[7] * 30830197480814L) - ((int128)tmp_q[8] * 85527218674374L) + ((int128)tmp_q[9] * 85362351060129L) - ((int128)tmp_q[10] * 49750569740116L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 49750569740116L) + ((int128)tmp_q[1] * 31531415086882L) - ((int128)tmp_q[2] * 34107822775284L) - ((int128)tmp_q[3] * 129388969751772L) + ((int128)tmp_q[4] * 104390065713229L) - ((int128)tmp_q[5] * 13504404300340L) - ((int128)tmp_q[6] * 80144319299779L) + ((int128)tmp_q[7] * 52001840031308L) - ((-((int128)tmp_q[8] * 30830197480814L) - ((int128)tmp_q[9] * 85527218674374L) + ((int128)tmp_q[10] * 85362351060129L)) * 3);
	tmp_zero[8] = ((int128)tmp_q[0] * 85362351060129L) - ((int128)tmp_q[1] * 49750569740116L) + ((int128)tmp_q[2] * 31531415086882L) - ((int128)tmp_q[3] * 34107822775284L) - ((int128)tmp_q[4] * 129388969751772L) + ((int128)tmp_q[5] * 104390065713229L) - ((int128)tmp_q[6] * 13504404300340L) - ((int128)tmp_q[7] * 80144319299779L) + ((int128)tmp_q[8] * 52001840031308L) - ((-((int128)tmp_q[9] * 30830197480814L) - ((int128)tmp_q[10] * 85527218674374L)) * 3);
	tmp_zero[9] = -((int128)tmp_q[0] * 85527218674374L) + ((int128)tmp_q[1] * 85362351060129L) - ((int128)tmp_q[2] * 49750569740116L) + ((int128)tmp_q[3] * 31531415086882L) - ((int128)tmp_q[4] * 34107822775284L) - ((int128)tmp_q[5] * 129388969751772L) + ((int128)tmp_q[6] * 104390065713229L) - ((int128)tmp_q[7] * 13504404300340L) - ((int128)tmp_q[8] * 80144319299779L) + ((int128)tmp_q[9] * 52001840031308L) + ((int128)tmp_q[10] * 92490592442442L);
	tmp_zero[10] = -((int128)tmp_q[0] * 30830197480814L) - ((int128)tmp_q[1] * 85527218674374L) + ((int128)tmp_q[2] * 85362351060129L) - ((int128)tmp_q[3] * 49750569740116L) + ((int128)tmp_q[4] * 31531415086882L) - ((int128)tmp_q[5] * 34107822775284L) - ((int128)tmp_q[6] * 129388969751772L) + ((int128)tmp_q[7] * 104390065713229L) - ((int128)tmp_q[8] * 13504404300340L) - ((int128)tmp_q[9] * 80144319299779L) + ((int128)tmp_q[10] * 52001840031308L);

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

