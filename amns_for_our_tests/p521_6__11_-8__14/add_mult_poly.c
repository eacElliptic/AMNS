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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15488618013881649969UL) + ((((uint64_t)op[1] * 4828324716223588908UL) + ((uint64_t)op[2] * 8510083420278861684UL) + ((uint64_t)op[3] * 12515823256918824852UL) + ((uint64_t)op[4] * 16637438444982903207UL) + ((uint64_t)op[5] * 17170857391390421488UL) + ((uint64_t)op[6] * 7607083196848096439UL) + ((uint64_t)op[7] * 12778897894695782642UL) + ((uint64_t)op[8] * 370051892187834129UL) + ((uint64_t)op[9] * 6940957202246346977UL) + ((uint64_t)op[10] * 3398190459821698682UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 3398190459821698682UL) + ((uint64_t)op[1] * 15488618013881649969UL) + ((((uint64_t)op[2] * 4828324716223588908UL) + ((uint64_t)op[3] * 8510083420278861684UL) + ((uint64_t)op[4] * 12515823256918824852UL) + ((uint64_t)op[5] * 16637438444982903207UL) + ((uint64_t)op[6] * 17170857391390421488UL) + ((uint64_t)op[7] * 7607083196848096439UL) + ((uint64_t)op[8] * 12778897894695782642UL) + ((uint64_t)op[9] * 370051892187834129UL) + ((uint64_t)op[10] * 6940957202246346977UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 6940957202246346977UL) + ((uint64_t)op[1] * 3398190459821698682UL) + ((uint64_t)op[2] * 15488618013881649969UL) + ((((uint64_t)op[3] * 4828324716223588908UL) + ((uint64_t)op[4] * 8510083420278861684UL) + ((uint64_t)op[5] * 12515823256918824852UL) + ((uint64_t)op[6] * 16637438444982903207UL) + ((uint64_t)op[7] * 17170857391390421488UL) + ((uint64_t)op[8] * 7607083196848096439UL) + ((uint64_t)op[9] * 12778897894695782642UL) + ((uint64_t)op[10] * 370051892187834129UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 370051892187834129UL) + ((uint64_t)op[1] * 6940957202246346977UL) + ((uint64_t)op[2] * 3398190459821698682UL) + ((uint64_t)op[3] * 15488618013881649969UL) + ((((uint64_t)op[4] * 4828324716223588908UL) + ((uint64_t)op[5] * 8510083420278861684UL) + ((uint64_t)op[6] * 12515823256918824852UL) + ((uint64_t)op[7] * 16637438444982903207UL) + ((uint64_t)op[8] * 17170857391390421488UL) + ((uint64_t)op[9] * 7607083196848096439UL) + ((uint64_t)op[10] * 12778897894695782642UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 12778897894695782642UL) + ((uint64_t)op[1] * 370051892187834129UL) + ((uint64_t)op[2] * 6940957202246346977UL) + ((uint64_t)op[3] * 3398190459821698682UL) + ((uint64_t)op[4] * 15488618013881649969UL) + ((((uint64_t)op[5] * 4828324716223588908UL) + ((uint64_t)op[6] * 8510083420278861684UL) + ((uint64_t)op[7] * 12515823256918824852UL) + ((uint64_t)op[8] * 16637438444982903207UL) + ((uint64_t)op[9] * 17170857391390421488UL) + ((uint64_t)op[10] * 7607083196848096439UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 7607083196848096439UL) + ((uint64_t)op[1] * 12778897894695782642UL) + ((uint64_t)op[2] * 370051892187834129UL) + ((uint64_t)op[3] * 6940957202246346977UL) + ((uint64_t)op[4] * 3398190459821698682UL) + ((uint64_t)op[5] * 15488618013881649969UL) + ((((uint64_t)op[6] * 4828324716223588908UL) + ((uint64_t)op[7] * 8510083420278861684UL) + ((uint64_t)op[8] * 12515823256918824852UL) + ((uint64_t)op[9] * 16637438444982903207UL) + ((uint64_t)op[10] * 17170857391390421488UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 17170857391390421488UL) + ((uint64_t)op[1] * 7607083196848096439UL) + ((uint64_t)op[2] * 12778897894695782642UL) + ((uint64_t)op[3] * 370051892187834129UL) + ((uint64_t)op[4] * 6940957202246346977UL) + ((uint64_t)op[5] * 3398190459821698682UL) + ((uint64_t)op[6] * 15488618013881649969UL) + ((((uint64_t)op[7] * 4828324716223588908UL) + ((uint64_t)op[8] * 8510083420278861684UL) + ((uint64_t)op[9] * 12515823256918824852UL) + ((uint64_t)op[10] * 16637438444982903207UL)) * 18446744073709551608);
	tmp_q[7] = ((uint64_t)op[0] * 16637438444982903207UL) + ((uint64_t)op[1] * 17170857391390421488UL) + ((uint64_t)op[2] * 7607083196848096439UL) + ((uint64_t)op[3] * 12778897894695782642UL) + ((uint64_t)op[4] * 370051892187834129UL) + ((uint64_t)op[5] * 6940957202246346977UL) + ((uint64_t)op[6] * 3398190459821698682UL) + ((uint64_t)op[7] * 15488618013881649969UL) + ((((uint64_t)op[8] * 4828324716223588908UL) + ((uint64_t)op[9] * 8510083420278861684UL) + ((uint64_t)op[10] * 12515823256918824852UL)) * 18446744073709551608);
	tmp_q[8] = ((uint64_t)op[0] * 12515823256918824852UL) + ((uint64_t)op[1] * 16637438444982903207UL) + ((uint64_t)op[2] * 17170857391390421488UL) + ((uint64_t)op[3] * 7607083196848096439UL) + ((uint64_t)op[4] * 12778897894695782642UL) + ((uint64_t)op[5] * 370051892187834129UL) + ((uint64_t)op[6] * 6940957202246346977UL) + ((uint64_t)op[7] * 3398190459821698682UL) + ((uint64_t)op[8] * 15488618013881649969UL) + ((((uint64_t)op[9] * 4828324716223588908UL) + ((uint64_t)op[10] * 8510083420278861684UL)) * 18446744073709551608);
	tmp_q[9] = ((uint64_t)op[0] * 8510083420278861684UL) + ((uint64_t)op[1] * 12515823256918824852UL) + ((uint64_t)op[2] * 16637438444982903207UL) + ((uint64_t)op[3] * 17170857391390421488UL) + ((uint64_t)op[4] * 7607083196848096439UL) + ((uint64_t)op[5] * 12778897894695782642UL) + ((uint64_t)op[6] * 370051892187834129UL) + ((uint64_t)op[7] * 6940957202246346977UL) + ((uint64_t)op[8] * 3398190459821698682UL) + ((uint64_t)op[9] * 15488618013881649969UL) + ((uint64_t)op[10] * 16713634491339943584UL);
	tmp_q[10] = ((uint64_t)op[0] * 4828324716223588908UL) + ((uint64_t)op[1] * 8510083420278861684UL) + ((uint64_t)op[2] * 12515823256918824852UL) + ((uint64_t)op[3] * 16637438444982903207UL) + ((uint64_t)op[4] * 17170857391390421488UL) + ((uint64_t)op[5] * 7607083196848096439UL) + ((uint64_t)op[6] * 12778897894695782642UL) + ((uint64_t)op[7] * 370051892187834129UL) + ((uint64_t)op[8] * 6940957202246346977UL) + ((uint64_t)op[9] * 3398190459821698682UL) + ((uint64_t)op[10] * 15488618013881649969UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 52128857076153L) - ((-((int128)tmp_q[1] * 110824474675708L) - ((int128)tmp_q[2] * 134494585353046L) - ((int128)tmp_q[3] * 79263308277598L) - ((int128)tmp_q[4] * 50218530210166L) + ((int128)tmp_q[5] * 56145655103172L) + ((int128)tmp_q[6] * 54249251549495L) - ((int128)tmp_q[7] * 7085824329399L) + ((int128)tmp_q[8] * 61142234045285L) - ((int128)tmp_q[9] * 114841692503099L) + ((int128)tmp_q[10] * 9350659158322L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 9350659158322L) - ((int128)tmp_q[1] * 52128857076153L) - ((-((int128)tmp_q[2] * 110824474675708L) - ((int128)tmp_q[3] * 134494585353046L) - ((int128)tmp_q[4] * 79263308277598L) - ((int128)tmp_q[5] * 50218530210166L) + ((int128)tmp_q[6] * 56145655103172L) + ((int128)tmp_q[7] * 54249251549495L) - ((int128)tmp_q[8] * 7085824329399L) + ((int128)tmp_q[9] * 61142234045285L) - ((int128)tmp_q[10] * 114841692503099L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 114841692503099L) + ((int128)tmp_q[1] * 9350659158322L) - ((int128)tmp_q[2] * 52128857076153L) - ((-((int128)tmp_q[3] * 110824474675708L) - ((int128)tmp_q[4] * 134494585353046L) - ((int128)tmp_q[5] * 79263308277598L) - ((int128)tmp_q[6] * 50218530210166L) + ((int128)tmp_q[7] * 56145655103172L) + ((int128)tmp_q[8] * 54249251549495L) - ((int128)tmp_q[9] * 7085824329399L) + ((int128)tmp_q[10] * 61142234045285L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 61142234045285L) - ((int128)tmp_q[1] * 114841692503099L) + ((int128)tmp_q[2] * 9350659158322L) - ((int128)tmp_q[3] * 52128857076153L) - ((-((int128)tmp_q[4] * 110824474675708L) - ((int128)tmp_q[5] * 134494585353046L) - ((int128)tmp_q[6] * 79263308277598L) - ((int128)tmp_q[7] * 50218530210166L) + ((int128)tmp_q[8] * 56145655103172L) + ((int128)tmp_q[9] * 54249251549495L) - ((int128)tmp_q[10] * 7085824329399L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 7085824329399L) + ((int128)tmp_q[1] * 61142234045285L) - ((int128)tmp_q[2] * 114841692503099L) + ((int128)tmp_q[3] * 9350659158322L) - ((int128)tmp_q[4] * 52128857076153L) - ((-((int128)tmp_q[5] * 110824474675708L) - ((int128)tmp_q[6] * 134494585353046L) - ((int128)tmp_q[7] * 79263308277598L) - ((int128)tmp_q[8] * 50218530210166L) + ((int128)tmp_q[9] * 56145655103172L) + ((int128)tmp_q[10] * 54249251549495L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 54249251549495L) - ((int128)tmp_q[1] * 7085824329399L) + ((int128)tmp_q[2] * 61142234045285L) - ((int128)tmp_q[3] * 114841692503099L) + ((int128)tmp_q[4] * 9350659158322L) - ((int128)tmp_q[5] * 52128857076153L) - ((-((int128)tmp_q[6] * 110824474675708L) - ((int128)tmp_q[7] * 134494585353046L) - ((int128)tmp_q[8] * 79263308277598L) - ((int128)tmp_q[9] * 50218530210166L) + ((int128)tmp_q[10] * 56145655103172L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 56145655103172L) + ((int128)tmp_q[1] * 54249251549495L) - ((int128)tmp_q[2] * 7085824329399L) + ((int128)tmp_q[3] * 61142234045285L) - ((int128)tmp_q[4] * 114841692503099L) + ((int128)tmp_q[5] * 9350659158322L) - ((int128)tmp_q[6] * 52128857076153L) - ((-((int128)tmp_q[7] * 110824474675708L) - ((int128)tmp_q[8] * 134494585353046L) - ((int128)tmp_q[9] * 79263308277598L) - ((int128)tmp_q[10] * 50218530210166L)) * 8);
	tmp_zero[7] = -((int128)tmp_q[0] * 50218530210166L) + ((int128)tmp_q[1] * 56145655103172L) + ((int128)tmp_q[2] * 54249251549495L) - ((int128)tmp_q[3] * 7085824329399L) + ((int128)tmp_q[4] * 61142234045285L) - ((int128)tmp_q[5] * 114841692503099L) + ((int128)tmp_q[6] * 9350659158322L) - ((int128)tmp_q[7] * 52128857076153L) - ((-((int128)tmp_q[8] * 110824474675708L) - ((int128)tmp_q[9] * 134494585353046L) - ((int128)tmp_q[10] * 79263308277598L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 79263308277598L) - ((int128)tmp_q[1] * 50218530210166L) + ((int128)tmp_q[2] * 56145655103172L) + ((int128)tmp_q[3] * 54249251549495L) - ((int128)tmp_q[4] * 7085824329399L) + ((int128)tmp_q[5] * 61142234045285L) - ((int128)tmp_q[6] * 114841692503099L) + ((int128)tmp_q[7] * 9350659158322L) - ((int128)tmp_q[8] * 52128857076153L) - ((-((int128)tmp_q[9] * 110824474675708L) - ((int128)tmp_q[10] * 134494585353046L)) * 8);
	tmp_zero[9] = -((int128)tmp_q[0] * 134494585353046L) - ((int128)tmp_q[1] * 79263308277598L) - ((int128)tmp_q[2] * 50218530210166L) + ((int128)tmp_q[3] * 56145655103172L) + ((int128)tmp_q[4] * 54249251549495L) - ((int128)tmp_q[5] * 7085824329399L) + ((int128)tmp_q[6] * 61142234045285L) - ((int128)tmp_q[7] * 114841692503099L) + ((int128)tmp_q[8] * 9350659158322L) - ((int128)tmp_q[9] * 52128857076153L) + ((int128)tmp_q[10] * 886595797405664L);
	tmp_zero[10] = -((int128)tmp_q[0] * 110824474675708L) - ((int128)tmp_q[1] * 134494585353046L) - ((int128)tmp_q[2] * 79263308277598L) - ((int128)tmp_q[3] * 50218530210166L) + ((int128)tmp_q[4] * 56145655103172L) + ((int128)tmp_q[5] * 54249251549495L) - ((int128)tmp_q[6] * 7085824329399L) + ((int128)tmp_q[7] * 61142234045285L) - ((int128)tmp_q[8] * 114841692503099L) + ((int128)tmp_q[9] * 9350659158322L) - ((int128)tmp_q[10] * 52128857076153L);

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

