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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18061724865669022941UL) + ((((uint64_t)op[1] * 2621030190369430568UL) + ((uint64_t)op[2] * 14209222310999897774UL) + ((uint64_t)op[3] * 5876322622772162443UL) + ((uint64_t)op[4] * 13345560754365163148UL) + ((uint64_t)op[5] * 10584907608698382375UL) + ((uint64_t)op[6] * 14200400205823003712UL) + ((uint64_t)op[7] * 11411136238901499872UL) + ((uint64_t)op[8] * 3948862291773971777UL) + ((uint64_t)op[9] * 495443341604413720UL) + ((uint64_t)op[10] * 4956923348049395765UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 4956923348049395765UL) + ((uint64_t)op[1] * 18061724865669022941UL) + ((((uint64_t)op[2] * 2621030190369430568UL) + ((uint64_t)op[3] * 14209222310999897774UL) + ((uint64_t)op[4] * 5876322622772162443UL) + ((uint64_t)op[5] * 13345560754365163148UL) + ((uint64_t)op[6] * 10584907608698382375UL) + ((uint64_t)op[7] * 14200400205823003712UL) + ((uint64_t)op[8] * 11411136238901499872UL) + ((uint64_t)op[9] * 3948862291773971777UL) + ((uint64_t)op[10] * 495443341604413720UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 495443341604413720UL) + ((uint64_t)op[1] * 4956923348049395765UL) + ((uint64_t)op[2] * 18061724865669022941UL) + ((((uint64_t)op[3] * 2621030190369430568UL) + ((uint64_t)op[4] * 14209222310999897774UL) + ((uint64_t)op[5] * 5876322622772162443UL) + ((uint64_t)op[6] * 13345560754365163148UL) + ((uint64_t)op[7] * 10584907608698382375UL) + ((uint64_t)op[8] * 14200400205823003712UL) + ((uint64_t)op[9] * 11411136238901499872UL) + ((uint64_t)op[10] * 3948862291773971777UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 3948862291773971777UL) + ((uint64_t)op[1] * 495443341604413720UL) + ((uint64_t)op[2] * 4956923348049395765UL) + ((uint64_t)op[3] * 18061724865669022941UL) + ((((uint64_t)op[4] * 2621030190369430568UL) + ((uint64_t)op[5] * 14209222310999897774UL) + ((uint64_t)op[6] * 5876322622772162443UL) + ((uint64_t)op[7] * 13345560754365163148UL) + ((uint64_t)op[8] * 10584907608698382375UL) + ((uint64_t)op[9] * 14200400205823003712UL) + ((uint64_t)op[10] * 11411136238901499872UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 11411136238901499872UL) + ((uint64_t)op[1] * 3948862291773971777UL) + ((uint64_t)op[2] * 495443341604413720UL) + ((uint64_t)op[3] * 4956923348049395765UL) + ((uint64_t)op[4] * 18061724865669022941UL) + ((((uint64_t)op[5] * 2621030190369430568UL) + ((uint64_t)op[6] * 14209222310999897774UL) + ((uint64_t)op[7] * 5876322622772162443UL) + ((uint64_t)op[8] * 13345560754365163148UL) + ((uint64_t)op[9] * 10584907608698382375UL) + ((uint64_t)op[10] * 14200400205823003712UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 14200400205823003712UL) + ((uint64_t)op[1] * 11411136238901499872UL) + ((uint64_t)op[2] * 3948862291773971777UL) + ((uint64_t)op[3] * 495443341604413720UL) + ((uint64_t)op[4] * 4956923348049395765UL) + ((uint64_t)op[5] * 18061724865669022941UL) + ((((uint64_t)op[6] * 2621030190369430568UL) + ((uint64_t)op[7] * 14209222310999897774UL) + ((uint64_t)op[8] * 5876322622772162443UL) + ((uint64_t)op[9] * 13345560754365163148UL) + ((uint64_t)op[10] * 10584907608698382375UL)) * 18446744073709551612);
	tmp_q[6] = ((uint64_t)op[0] * 10584907608698382375UL) + ((uint64_t)op[1] * 14200400205823003712UL) + ((uint64_t)op[2] * 11411136238901499872UL) + ((uint64_t)op[3] * 3948862291773971777UL) + ((uint64_t)op[4] * 495443341604413720UL) + ((uint64_t)op[5] * 4956923348049395765UL) + ((uint64_t)op[6] * 18061724865669022941UL) + ((((uint64_t)op[7] * 2621030190369430568UL) + ((uint64_t)op[8] * 14209222310999897774UL) + ((uint64_t)op[9] * 5876322622772162443UL) + ((uint64_t)op[10] * 13345560754365163148UL)) * 18446744073709551612);
	tmp_q[7] = ((uint64_t)op[0] * 13345560754365163148UL) + ((uint64_t)op[1] * 10584907608698382375UL) + ((uint64_t)op[2] * 14200400205823003712UL) + ((uint64_t)op[3] * 11411136238901499872UL) + ((uint64_t)op[4] * 3948862291773971777UL) + ((uint64_t)op[5] * 495443341604413720UL) + ((uint64_t)op[6] * 4956923348049395765UL) + ((uint64_t)op[7] * 18061724865669022941UL) + ((((uint64_t)op[8] * 2621030190369430568UL) + ((uint64_t)op[9] * 14209222310999897774UL) + ((uint64_t)op[10] * 5876322622772162443UL)) * 18446744073709551612);
	tmp_q[8] = ((uint64_t)op[0] * 5876322622772162443UL) + ((uint64_t)op[1] * 13345560754365163148UL) + ((uint64_t)op[2] * 10584907608698382375UL) + ((uint64_t)op[3] * 14200400205823003712UL) + ((uint64_t)op[4] * 11411136238901499872UL) + ((uint64_t)op[5] * 3948862291773971777UL) + ((uint64_t)op[6] * 495443341604413720UL) + ((uint64_t)op[7] * 4956923348049395765UL) + ((uint64_t)op[8] * 18061724865669022941UL) + ((((uint64_t)op[9] * 2621030190369430568UL) + ((uint64_t)op[10] * 14209222310999897774UL)) * 18446744073709551612);
	tmp_q[9] = ((uint64_t)op[0] * 14209222310999897774UL) + ((uint64_t)op[1] * 5876322622772162443UL) + ((uint64_t)op[2] * 13345560754365163148UL) + ((uint64_t)op[3] * 10584907608698382375UL) + ((uint64_t)op[4] * 14200400205823003712UL) + ((uint64_t)op[5] * 11411136238901499872UL) + ((uint64_t)op[6] * 3948862291773971777UL) + ((uint64_t)op[7] * 495443341604413720UL) + ((uint64_t)op[8] * 4956923348049395765UL) + ((uint64_t)op[9] * 18061724865669022941UL) + ((uint64_t)op[10] * 7962623312231829344UL);
	tmp_q[10] = ((uint64_t)op[0] * 2621030190369430568UL) + ((uint64_t)op[1] * 14209222310999897774UL) + ((uint64_t)op[2] * 5876322622772162443UL) + ((uint64_t)op[3] * 13345560754365163148UL) + ((uint64_t)op[4] * 10584907608698382375UL) + ((uint64_t)op[5] * 14200400205823003712UL) + ((uint64_t)op[6] * 11411136238901499872UL) + ((uint64_t)op[7] * 3948862291773971777UL) + ((uint64_t)op[8] * 495443341604413720UL) + ((uint64_t)op[9] * 4956923348049395765UL) + ((uint64_t)op[10] * 18061724865669022941UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 23011257119489L) - ((-((int128)tmp_q[1] * 94856605580542L) - ((int128)tmp_q[2] * 96081295137959L) + ((int128)tmp_q[3] * 11304059522671L) - ((int128)tmp_q[4] * 24707607174577L) - ((int128)tmp_q[5] * 59764687414799L) - ((int128)tmp_q[6] * 37587564078240L) - ((int128)tmp_q[7] * 77342149871519L) + ((int128)tmp_q[8] * 117996134200086L) + ((int128)tmp_q[9] * 12788458602111L) + ((int128)tmp_q[10] * 42527414772037L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 42527414772037L) - ((int128)tmp_q[1] * 23011257119489L) - ((-((int128)tmp_q[2] * 94856605580542L) - ((int128)tmp_q[3] * 96081295137959L) + ((int128)tmp_q[4] * 11304059522671L) - ((int128)tmp_q[5] * 24707607174577L) - ((int128)tmp_q[6] * 59764687414799L) - ((int128)tmp_q[7] * 37587564078240L) - ((int128)tmp_q[8] * 77342149871519L) + ((int128)tmp_q[9] * 117996134200086L) + ((int128)tmp_q[10] * 12788458602111L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 12788458602111L) + ((int128)tmp_q[1] * 42527414772037L) - ((int128)tmp_q[2] * 23011257119489L) - ((-((int128)tmp_q[3] * 94856605580542L) - ((int128)tmp_q[4] * 96081295137959L) + ((int128)tmp_q[5] * 11304059522671L) - ((int128)tmp_q[6] * 24707607174577L) - ((int128)tmp_q[7] * 59764687414799L) - ((int128)tmp_q[8] * 37587564078240L) - ((int128)tmp_q[9] * 77342149871519L) + ((int128)tmp_q[10] * 117996134200086L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 117996134200086L) + ((int128)tmp_q[1] * 12788458602111L) + ((int128)tmp_q[2] * 42527414772037L) - ((int128)tmp_q[3] * 23011257119489L) - ((-((int128)tmp_q[4] * 94856605580542L) - ((int128)tmp_q[5] * 96081295137959L) + ((int128)tmp_q[6] * 11304059522671L) - ((int128)tmp_q[7] * 24707607174577L) - ((int128)tmp_q[8] * 59764687414799L) - ((int128)tmp_q[9] * 37587564078240L) - ((int128)tmp_q[10] * 77342149871519L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 77342149871519L) + ((int128)tmp_q[1] * 117996134200086L) + ((int128)tmp_q[2] * 12788458602111L) + ((int128)tmp_q[3] * 42527414772037L) - ((int128)tmp_q[4] * 23011257119489L) - ((-((int128)tmp_q[5] * 94856605580542L) - ((int128)tmp_q[6] * 96081295137959L) + ((int128)tmp_q[7] * 11304059522671L) - ((int128)tmp_q[8] * 24707607174577L) - ((int128)tmp_q[9] * 59764687414799L) - ((int128)tmp_q[10] * 37587564078240L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 37587564078240L) - ((int128)tmp_q[1] * 77342149871519L) + ((int128)tmp_q[2] * 117996134200086L) + ((int128)tmp_q[3] * 12788458602111L) + ((int128)tmp_q[4] * 42527414772037L) - ((int128)tmp_q[5] * 23011257119489L) - ((-((int128)tmp_q[6] * 94856605580542L) - ((int128)tmp_q[7] * 96081295137959L) + ((int128)tmp_q[8] * 11304059522671L) - ((int128)tmp_q[9] * 24707607174577L) - ((int128)tmp_q[10] * 59764687414799L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 59764687414799L) - ((int128)tmp_q[1] * 37587564078240L) - ((int128)tmp_q[2] * 77342149871519L) + ((int128)tmp_q[3] * 117996134200086L) + ((int128)tmp_q[4] * 12788458602111L) + ((int128)tmp_q[5] * 42527414772037L) - ((int128)tmp_q[6] * 23011257119489L) - ((-((int128)tmp_q[7] * 94856605580542L) - ((int128)tmp_q[8] * 96081295137959L) + ((int128)tmp_q[9] * 11304059522671L) - ((int128)tmp_q[10] * 24707607174577L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 24707607174577L) - ((int128)tmp_q[1] * 59764687414799L) - ((int128)tmp_q[2] * 37587564078240L) - ((int128)tmp_q[3] * 77342149871519L) + ((int128)tmp_q[4] * 117996134200086L) + ((int128)tmp_q[5] * 12788458602111L) + ((int128)tmp_q[6] * 42527414772037L) - ((int128)tmp_q[7] * 23011257119489L) - ((-((int128)tmp_q[8] * 94856605580542L) - ((int128)tmp_q[9] * 96081295137959L) + ((int128)tmp_q[10] * 11304059522671L)) * 4);
	tmp_zero[8] = ((int128)tmp_q[0] * 11304059522671L) - ((int128)tmp_q[1] * 24707607174577L) - ((int128)tmp_q[2] * 59764687414799L) - ((int128)tmp_q[3] * 37587564078240L) - ((int128)tmp_q[4] * 77342149871519L) + ((int128)tmp_q[5] * 117996134200086L) + ((int128)tmp_q[6] * 12788458602111L) + ((int128)tmp_q[7] * 42527414772037L) - ((int128)tmp_q[8] * 23011257119489L) - ((-((int128)tmp_q[9] * 94856605580542L) - ((int128)tmp_q[10] * 96081295137959L)) * 4);
	tmp_zero[9] = -((int128)tmp_q[0] * 96081295137959L) + ((int128)tmp_q[1] * 11304059522671L) - ((int128)tmp_q[2] * 24707607174577L) - ((int128)tmp_q[3] * 59764687414799L) - ((int128)tmp_q[4] * 37587564078240L) - ((int128)tmp_q[5] * 77342149871519L) + ((int128)tmp_q[6] * 117996134200086L) + ((int128)tmp_q[7] * 12788458602111L) + ((int128)tmp_q[8] * 42527414772037L) - ((int128)tmp_q[9] * 23011257119489L) + ((int128)tmp_q[10] * 379426422322168L);
	tmp_zero[10] = -((int128)tmp_q[0] * 94856605580542L) - ((int128)tmp_q[1] * 96081295137959L) + ((int128)tmp_q[2] * 11304059522671L) - ((int128)tmp_q[3] * 24707607174577L) - ((int128)tmp_q[4] * 59764687414799L) - ((int128)tmp_q[5] * 37587564078240L) - ((int128)tmp_q[6] * 77342149871519L) + ((int128)tmp_q[7] * 117996134200086L) + ((int128)tmp_q[8] * 12788458602111L) + ((int128)tmp_q[9] * 42527414772037L) - ((int128)tmp_q[10] * 23011257119489L);

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

