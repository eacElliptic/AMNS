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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[9] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14388014452064965306UL) + ((((uint64_t)op[1] * 14500254052397834256UL) + ((uint64_t)op[2] * 1522052409609893955UL) + ((uint64_t)op[3] * 9563123063442680770UL) + ((uint64_t)op[4] * 11395344622979772653UL) + ((uint64_t)op[5] * 9234355762929124086UL) + ((uint64_t)op[6] * 12615267642993295766UL) + ((uint64_t)op[7] * 1119498230262676579UL) + ((uint64_t)op[8] * 17388357529840499673UL) + ((uint64_t)op[9] * 1728183534309831885UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 1728183534309831885UL) + ((uint64_t)op[1] * 14388014452064965306UL) + ((((uint64_t)op[2] * 14500254052397834256UL) + ((uint64_t)op[3] * 1522052409609893955UL) + ((uint64_t)op[4] * 9563123063442680770UL) + ((uint64_t)op[5] * 11395344622979772653UL) + ((uint64_t)op[6] * 9234355762929124086UL) + ((uint64_t)op[7] * 12615267642993295766UL) + ((uint64_t)op[8] * 1119498230262676579UL) + ((uint64_t)op[9] * 17388357529840499673UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 17388357529840499673UL) + ((uint64_t)op[1] * 1728183534309831885UL) + ((uint64_t)op[2] * 14388014452064965306UL) + ((((uint64_t)op[3] * 14500254052397834256UL) + ((uint64_t)op[4] * 1522052409609893955UL) + ((uint64_t)op[5] * 9563123063442680770UL) + ((uint64_t)op[6] * 11395344622979772653UL) + ((uint64_t)op[7] * 9234355762929124086UL) + ((uint64_t)op[8] * 12615267642993295766UL) + ((uint64_t)op[9] * 1119498230262676579UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 1119498230262676579UL) + ((uint64_t)op[1] * 17388357529840499673UL) + ((uint64_t)op[2] * 1728183534309831885UL) + ((uint64_t)op[3] * 14388014452064965306UL) + ((((uint64_t)op[4] * 14500254052397834256UL) + ((uint64_t)op[5] * 1522052409609893955UL) + ((uint64_t)op[6] * 9563123063442680770UL) + ((uint64_t)op[7] * 11395344622979772653UL) + ((uint64_t)op[8] * 9234355762929124086UL) + ((uint64_t)op[9] * 12615267642993295766UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 12615267642993295766UL) + ((uint64_t)op[1] * 1119498230262676579UL) + ((uint64_t)op[2] * 17388357529840499673UL) + ((uint64_t)op[3] * 1728183534309831885UL) + ((uint64_t)op[4] * 14388014452064965306UL) + ((((uint64_t)op[5] * 14500254052397834256UL) + ((uint64_t)op[6] * 1522052409609893955UL) + ((uint64_t)op[7] * 9563123063442680770UL) + ((uint64_t)op[8] * 11395344622979772653UL) + ((uint64_t)op[9] * 9234355762929124086UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 9234355762929124086UL) + ((uint64_t)op[1] * 12615267642993295766UL) + ((uint64_t)op[2] * 1119498230262676579UL) + ((uint64_t)op[3] * 17388357529840499673UL) + ((uint64_t)op[4] * 1728183534309831885UL) + ((uint64_t)op[5] * 14388014452064965306UL) + ((((uint64_t)op[6] * 14500254052397834256UL) + ((uint64_t)op[7] * 1522052409609893955UL) + ((uint64_t)op[8] * 9563123063442680770UL) + ((uint64_t)op[9] * 11395344622979772653UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 11395344622979772653UL) + ((uint64_t)op[1] * 9234355762929124086UL) + ((uint64_t)op[2] * 12615267642993295766UL) + ((uint64_t)op[3] * 1119498230262676579UL) + ((uint64_t)op[4] * 17388357529840499673UL) + ((uint64_t)op[5] * 1728183534309831885UL) + ((uint64_t)op[6] * 14388014452064965306UL) + ((((uint64_t)op[7] * 14500254052397834256UL) + ((uint64_t)op[8] * 1522052409609893955UL) + ((uint64_t)op[9] * 9563123063442680770UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 9563123063442680770UL) + ((uint64_t)op[1] * 11395344622979772653UL) + ((uint64_t)op[2] * 9234355762929124086UL) + ((uint64_t)op[3] * 12615267642993295766UL) + ((uint64_t)op[4] * 1119498230262676579UL) + ((uint64_t)op[5] * 17388357529840499673UL) + ((uint64_t)op[6] * 1728183534309831885UL) + ((uint64_t)op[7] * 14388014452064965306UL) + ((((uint64_t)op[8] * 14500254052397834256UL) + ((uint64_t)op[9] * 1522052409609893955UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 1522052409609893955UL) + ((uint64_t)op[1] * 9563123063442680770UL) + ((uint64_t)op[2] * 11395344622979772653UL) + ((uint64_t)op[3] * 9234355762929124086UL) + ((uint64_t)op[4] * 12615267642993295766UL) + ((uint64_t)op[5] * 1119498230262676579UL) + ((uint64_t)op[6] * 17388357529840499673UL) + ((uint64_t)op[7] * 1728183534309831885UL) + ((uint64_t)op[8] * 14388014452064965306UL) + ((uint64_t)op[9] * 6607274009774399536UL);
	tmp_q[9] = ((uint64_t)op[0] * 14500254052397834256UL) + ((uint64_t)op[1] * 1522052409609893955UL) + ((uint64_t)op[2] * 9563123063442680770UL) + ((uint64_t)op[3] * 11395344622979772653UL) + ((uint64_t)op[4] * 9234355762929124086UL) + ((uint64_t)op[5] * 12615267642993295766UL) + ((uint64_t)op[6] * 1119498230262676579UL) + ((uint64_t)op[7] * 17388357529840499673UL) + ((uint64_t)op[8] * 1728183534309831885UL) + ((uint64_t)op[9] * 14388014452064965306UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1848367198612744L) + ((((int128)tmp_q[1] * 345338518145379L) - ((int128)tmp_q[2] * 658297748985759L) - ((int128)tmp_q[3] * 350310328176447L) + ((int128)tmp_q[4] * 1844402004255788L) - ((int128)tmp_q[5] * 1035536181941808L) + ((int128)tmp_q[6] * 1259746539004481L) - ((int128)tmp_q[7] * 2360625965467582L) + ((int128)tmp_q[8] * 75759054970051L) - ((int128)tmp_q[9] * 1489761720247418L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 1489761720247418L) + ((int128)tmp_q[1] * 1848367198612744L) + ((((int128)tmp_q[2] * 345338518145379L) - ((int128)tmp_q[3] * 658297748985759L) - ((int128)tmp_q[4] * 350310328176447L) + ((int128)tmp_q[5] * 1844402004255788L) - ((int128)tmp_q[6] * 1035536181941808L) + ((int128)tmp_q[7] * 1259746539004481L) - ((int128)tmp_q[8] * 2360625965467582L) + ((int128)tmp_q[9] * 75759054970051L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 75759054970051L) - ((int128)tmp_q[1] * 1489761720247418L) + ((int128)tmp_q[2] * 1848367198612744L) + ((((int128)tmp_q[3] * 345338518145379L) - ((int128)tmp_q[4] * 658297748985759L) - ((int128)tmp_q[5] * 350310328176447L) + ((int128)tmp_q[6] * 1844402004255788L) - ((int128)tmp_q[7] * 1035536181941808L) + ((int128)tmp_q[8] * 1259746539004481L) - ((int128)tmp_q[9] * 2360625965467582L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 2360625965467582L) + ((int128)tmp_q[1] * 75759054970051L) - ((int128)tmp_q[2] * 1489761720247418L) + ((int128)tmp_q[3] * 1848367198612744L) + ((((int128)tmp_q[4] * 345338518145379L) - ((int128)tmp_q[5] * 658297748985759L) - ((int128)tmp_q[6] * 350310328176447L) + ((int128)tmp_q[7] * 1844402004255788L) - ((int128)tmp_q[8] * 1035536181941808L) + ((int128)tmp_q[9] * 1259746539004481L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1259746539004481L) - ((int128)tmp_q[1] * 2360625965467582L) + ((int128)tmp_q[2] * 75759054970051L) - ((int128)tmp_q[3] * 1489761720247418L) + ((int128)tmp_q[4] * 1848367198612744L) + ((((int128)tmp_q[5] * 345338518145379L) - ((int128)tmp_q[6] * 658297748985759L) - ((int128)tmp_q[7] * 350310328176447L) + ((int128)tmp_q[8] * 1844402004255788L) - ((int128)tmp_q[9] * 1035536181941808L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 1035536181941808L) + ((int128)tmp_q[1] * 1259746539004481L) - ((int128)tmp_q[2] * 2360625965467582L) + ((int128)tmp_q[3] * 75759054970051L) - ((int128)tmp_q[4] * 1489761720247418L) + ((int128)tmp_q[5] * 1848367198612744L) + ((((int128)tmp_q[6] * 345338518145379L) - ((int128)tmp_q[7] * 658297748985759L) - ((int128)tmp_q[8] * 350310328176447L) + ((int128)tmp_q[9] * 1844402004255788L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 1844402004255788L) - ((int128)tmp_q[1] * 1035536181941808L) + ((int128)tmp_q[2] * 1259746539004481L) - ((int128)tmp_q[3] * 2360625965467582L) + ((int128)tmp_q[4] * 75759054970051L) - ((int128)tmp_q[5] * 1489761720247418L) + ((int128)tmp_q[6] * 1848367198612744L) + ((((int128)tmp_q[7] * 345338518145379L) - ((int128)tmp_q[8] * 658297748985759L) - ((int128)tmp_q[9] * 350310328176447L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 350310328176447L) + ((int128)tmp_q[1] * 1844402004255788L) - ((int128)tmp_q[2] * 1035536181941808L) + ((int128)tmp_q[3] * 1259746539004481L) - ((int128)tmp_q[4] * 2360625965467582L) + ((int128)tmp_q[5] * 75759054970051L) - ((int128)tmp_q[6] * 1489761720247418L) + ((int128)tmp_q[7] * 1848367198612744L) + ((((int128)tmp_q[8] * 345338518145379L) - ((int128)tmp_q[9] * 658297748985759L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 658297748985759L) - ((int128)tmp_q[1] * 350310328176447L) + ((int128)tmp_q[2] * 1844402004255788L) - ((int128)tmp_q[3] * 1035536181941808L) + ((int128)tmp_q[4] * 1259746539004481L) - ((int128)tmp_q[5] * 2360625965467582L) + ((int128)tmp_q[6] * 75759054970051L) - ((int128)tmp_q[7] * 1489761720247418L) + ((int128)tmp_q[8] * 1848367198612744L) + ((int128)tmp_q[9] * 1036015554436137L);
	tmp_zero[9] = ((int128)tmp_q[0] * 345338518145379L) - ((int128)tmp_q[1] * 658297748985759L) - ((int128)tmp_q[2] * 350310328176447L) + ((int128)tmp_q[3] * 1844402004255788L) - ((int128)tmp_q[4] * 1035536181941808L) + ((int128)tmp_q[5] * 1259746539004481L) - ((int128)tmp_q[6] * 2360625965467582L) + ((int128)tmp_q[7] * 75759054970051L) - ((int128)tmp_q[8] * 1489761720247418L) + ((int128)tmp_q[9] * 1848367198612744L);

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
}

