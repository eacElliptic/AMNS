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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9853283285484146145UL) + ((((uint64_t)op[1] * 9128513613314683837UL) + ((uint64_t)op[2] * 16980825150343992163UL) + ((uint64_t)op[3] * 4686180917643436506UL) + ((uint64_t)op[4] * 12083967061098105209UL) + ((uint64_t)op[5] * 6891257908598056655UL) + ((uint64_t)op[6] * 16965687210961733078UL) + ((uint64_t)op[7] * 3821938256574841625UL) + ((uint64_t)op[8] * 11757085864734567060UL) + ((uint64_t)op[9] * 2637623599782042593UL) + ((uint64_t)op[10] * 10834384576476133944UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 10834384576476133944UL) + ((uint64_t)op[1] * 9853283285484146145UL) + ((((uint64_t)op[2] * 9128513613314683837UL) + ((uint64_t)op[3] * 16980825150343992163UL) + ((uint64_t)op[4] * 4686180917643436506UL) + ((uint64_t)op[5] * 12083967061098105209UL) + ((uint64_t)op[6] * 6891257908598056655UL) + ((uint64_t)op[7] * 16965687210961733078UL) + ((uint64_t)op[8] * 3821938256574841625UL) + ((uint64_t)op[9] * 11757085864734567060UL) + ((uint64_t)op[10] * 2637623599782042593UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 2637623599782042593UL) + ((uint64_t)op[1] * 10834384576476133944UL) + ((uint64_t)op[2] * 9853283285484146145UL) + ((((uint64_t)op[3] * 9128513613314683837UL) + ((uint64_t)op[4] * 16980825150343992163UL) + ((uint64_t)op[5] * 4686180917643436506UL) + ((uint64_t)op[6] * 12083967061098105209UL) + ((uint64_t)op[7] * 6891257908598056655UL) + ((uint64_t)op[8] * 16965687210961733078UL) + ((uint64_t)op[9] * 3821938256574841625UL) + ((uint64_t)op[10] * 11757085864734567060UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 11757085864734567060UL) + ((uint64_t)op[1] * 2637623599782042593UL) + ((uint64_t)op[2] * 10834384576476133944UL) + ((uint64_t)op[3] * 9853283285484146145UL) + ((((uint64_t)op[4] * 9128513613314683837UL) + ((uint64_t)op[5] * 16980825150343992163UL) + ((uint64_t)op[6] * 4686180917643436506UL) + ((uint64_t)op[7] * 12083967061098105209UL) + ((uint64_t)op[8] * 6891257908598056655UL) + ((uint64_t)op[9] * 16965687210961733078UL) + ((uint64_t)op[10] * 3821938256574841625UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 3821938256574841625UL) + ((uint64_t)op[1] * 11757085864734567060UL) + ((uint64_t)op[2] * 2637623599782042593UL) + ((uint64_t)op[3] * 10834384576476133944UL) + ((uint64_t)op[4] * 9853283285484146145UL) + ((((uint64_t)op[5] * 9128513613314683837UL) + ((uint64_t)op[6] * 16980825150343992163UL) + ((uint64_t)op[7] * 4686180917643436506UL) + ((uint64_t)op[8] * 12083967061098105209UL) + ((uint64_t)op[9] * 6891257908598056655UL) + ((uint64_t)op[10] * 16965687210961733078UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 16965687210961733078UL) + ((uint64_t)op[1] * 3821938256574841625UL) + ((uint64_t)op[2] * 11757085864734567060UL) + ((uint64_t)op[3] * 2637623599782042593UL) + ((uint64_t)op[4] * 10834384576476133944UL) + ((uint64_t)op[5] * 9853283285484146145UL) + ((((uint64_t)op[6] * 9128513613314683837UL) + ((uint64_t)op[7] * 16980825150343992163UL) + ((uint64_t)op[8] * 4686180917643436506UL) + ((uint64_t)op[9] * 12083967061098105209UL) + ((uint64_t)op[10] * 6891257908598056655UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 6891257908598056655UL) + ((uint64_t)op[1] * 16965687210961733078UL) + ((uint64_t)op[2] * 3821938256574841625UL) + ((uint64_t)op[3] * 11757085864734567060UL) + ((uint64_t)op[4] * 2637623599782042593UL) + ((uint64_t)op[5] * 10834384576476133944UL) + ((uint64_t)op[6] * 9853283285484146145UL) + ((((uint64_t)op[7] * 9128513613314683837UL) + ((uint64_t)op[8] * 16980825150343992163UL) + ((uint64_t)op[9] * 4686180917643436506UL) + ((uint64_t)op[10] * 12083967061098105209UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 12083967061098105209UL) + ((uint64_t)op[1] * 6891257908598056655UL) + ((uint64_t)op[2] * 16965687210961733078UL) + ((uint64_t)op[3] * 3821938256574841625UL) + ((uint64_t)op[4] * 11757085864734567060UL) + ((uint64_t)op[5] * 2637623599782042593UL) + ((uint64_t)op[6] * 10834384576476133944UL) + ((uint64_t)op[7] * 9853283285484146145UL) + ((((uint64_t)op[8] * 9128513613314683837UL) + ((uint64_t)op[9] * 16980825150343992163UL) + ((uint64_t)op[10] * 4686180917643436506UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 4686180917643436506UL) + ((uint64_t)op[1] * 12083967061098105209UL) + ((uint64_t)op[2] * 6891257908598056655UL) + ((uint64_t)op[3] * 16965687210961733078UL) + ((uint64_t)op[4] * 3821938256574841625UL) + ((uint64_t)op[5] * 11757085864734567060UL) + ((uint64_t)op[6] * 2637623599782042593UL) + ((uint64_t)op[7] * 10834384576476133944UL) + ((uint64_t)op[8] * 9853283285484146145UL) + ((((uint64_t)op[9] * 9128513613314683837UL) + ((uint64_t)op[10] * 16980825150343992163UL)) * 3);
	tmp_q[9] = ((uint64_t)op[0] * 16980825150343992163UL) + ((uint64_t)op[1] * 4686180917643436506UL) + ((uint64_t)op[2] * 12083967061098105209UL) + ((uint64_t)op[3] * 6891257908598056655UL) + ((uint64_t)op[4] * 16965687210961733078UL) + ((uint64_t)op[5] * 3821938256574841625UL) + ((uint64_t)op[6] * 11757085864734567060UL) + ((uint64_t)op[7] * 2637623599782042593UL) + ((uint64_t)op[8] * 10834384576476133944UL) + ((uint64_t)op[9] * 9853283285484146145UL) + ((uint64_t)op[10] * 8938796766234499895UL);
	tmp_q[10] = ((uint64_t)op[0] * 9128513613314683837UL) + ((uint64_t)op[1] * 16980825150343992163UL) + ((uint64_t)op[2] * 4686180917643436506UL) + ((uint64_t)op[3] * 12083967061098105209UL) + ((uint64_t)op[4] * 6891257908598056655UL) + ((uint64_t)op[5] * 16965687210961733078UL) + ((uint64_t)op[6] * 3821938256574841625UL) + ((uint64_t)op[7] * 11757085864734567060UL) + ((uint64_t)op[8] * 2637623599782042593UL) + ((uint64_t)op[9] * 10834384576476133944UL) + ((uint64_t)op[10] * 9853283285484146145UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 63800393038753L) + ((((int128)tmp_q[1] * 100234367967271L) - ((int128)tmp_q[2] * 53281151049492L) + ((int128)tmp_q[3] * 14740502324586L) + ((int128)tmp_q[4] * 96801719312474L) - ((int128)tmp_q[5] * 95679293498708L) + ((int128)tmp_q[6] * 34277528353011L) + ((int128)tmp_q[7] * 48450626635358L) + ((int128)tmp_q[8] * 18157203451227L) + ((int128)tmp_q[9] * 71845021202114L) + ((int128)tmp_q[10] * 59664278664703L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 59664278664703L) - ((int128)tmp_q[1] * 63800393038753L) + ((((int128)tmp_q[2] * 100234367967271L) - ((int128)tmp_q[3] * 53281151049492L) + ((int128)tmp_q[4] * 14740502324586L) + ((int128)tmp_q[5] * 96801719312474L) - ((int128)tmp_q[6] * 95679293498708L) + ((int128)tmp_q[7] * 34277528353011L) + ((int128)tmp_q[8] * 48450626635358L) + ((int128)tmp_q[9] * 18157203451227L) + ((int128)tmp_q[10] * 71845021202114L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 71845021202114L) + ((int128)tmp_q[1] * 59664278664703L) - ((int128)tmp_q[2] * 63800393038753L) + ((((int128)tmp_q[3] * 100234367967271L) - ((int128)tmp_q[4] * 53281151049492L) + ((int128)tmp_q[5] * 14740502324586L) + ((int128)tmp_q[6] * 96801719312474L) - ((int128)tmp_q[7] * 95679293498708L) + ((int128)tmp_q[8] * 34277528353011L) + ((int128)tmp_q[9] * 48450626635358L) + ((int128)tmp_q[10] * 18157203451227L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 18157203451227L) + ((int128)tmp_q[1] * 71845021202114L) + ((int128)tmp_q[2] * 59664278664703L) - ((int128)tmp_q[3] * 63800393038753L) + ((((int128)tmp_q[4] * 100234367967271L) - ((int128)tmp_q[5] * 53281151049492L) + ((int128)tmp_q[6] * 14740502324586L) + ((int128)tmp_q[7] * 96801719312474L) - ((int128)tmp_q[8] * 95679293498708L) + ((int128)tmp_q[9] * 34277528353011L) + ((int128)tmp_q[10] * 48450626635358L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 48450626635358L) + ((int128)tmp_q[1] * 18157203451227L) + ((int128)tmp_q[2] * 71845021202114L) + ((int128)tmp_q[3] * 59664278664703L) - ((int128)tmp_q[4] * 63800393038753L) + ((((int128)tmp_q[5] * 100234367967271L) - ((int128)tmp_q[6] * 53281151049492L) + ((int128)tmp_q[7] * 14740502324586L) + ((int128)tmp_q[8] * 96801719312474L) - ((int128)tmp_q[9] * 95679293498708L) + ((int128)tmp_q[10] * 34277528353011L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 34277528353011L) + ((int128)tmp_q[1] * 48450626635358L) + ((int128)tmp_q[2] * 18157203451227L) + ((int128)tmp_q[3] * 71845021202114L) + ((int128)tmp_q[4] * 59664278664703L) - ((int128)tmp_q[5] * 63800393038753L) + ((((int128)tmp_q[6] * 100234367967271L) - ((int128)tmp_q[7] * 53281151049492L) + ((int128)tmp_q[8] * 14740502324586L) + ((int128)tmp_q[9] * 96801719312474L) - ((int128)tmp_q[10] * 95679293498708L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 95679293498708L) + ((int128)tmp_q[1] * 34277528353011L) + ((int128)tmp_q[2] * 48450626635358L) + ((int128)tmp_q[3] * 18157203451227L) + ((int128)tmp_q[4] * 71845021202114L) + ((int128)tmp_q[5] * 59664278664703L) - ((int128)tmp_q[6] * 63800393038753L) + ((((int128)tmp_q[7] * 100234367967271L) - ((int128)tmp_q[8] * 53281151049492L) + ((int128)tmp_q[9] * 14740502324586L) + ((int128)tmp_q[10] * 96801719312474L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 96801719312474L) - ((int128)tmp_q[1] * 95679293498708L) + ((int128)tmp_q[2] * 34277528353011L) + ((int128)tmp_q[3] * 48450626635358L) + ((int128)tmp_q[4] * 18157203451227L) + ((int128)tmp_q[5] * 71845021202114L) + ((int128)tmp_q[6] * 59664278664703L) - ((int128)tmp_q[7] * 63800393038753L) + ((((int128)tmp_q[8] * 100234367967271L) - ((int128)tmp_q[9] * 53281151049492L) + ((int128)tmp_q[10] * 14740502324586L)) * 3);
	tmp_zero[8] = ((int128)tmp_q[0] * 14740502324586L) + ((int128)tmp_q[1] * 96801719312474L) - ((int128)tmp_q[2] * 95679293498708L) + ((int128)tmp_q[3] * 34277528353011L) + ((int128)tmp_q[4] * 48450626635358L) + ((int128)tmp_q[5] * 18157203451227L) + ((int128)tmp_q[6] * 71845021202114L) + ((int128)tmp_q[7] * 59664278664703L) - ((int128)tmp_q[8] * 63800393038753L) + ((((int128)tmp_q[9] * 100234367967271L) - ((int128)tmp_q[10] * 53281151049492L)) * 3);
	tmp_zero[9] = -((int128)tmp_q[0] * 53281151049492L) + ((int128)tmp_q[1] * 14740502324586L) + ((int128)tmp_q[2] * 96801719312474L) - ((int128)tmp_q[3] * 95679293498708L) + ((int128)tmp_q[4] * 34277528353011L) + ((int128)tmp_q[5] * 48450626635358L) + ((int128)tmp_q[6] * 18157203451227L) + ((int128)tmp_q[7] * 71845021202114L) + ((int128)tmp_q[8] * 59664278664703L) - ((int128)tmp_q[9] * 63800393038753L) + ((int128)tmp_q[10] * 300703103901813L);
	tmp_zero[10] = ((int128)tmp_q[0] * 100234367967271L) - ((int128)tmp_q[1] * 53281151049492L) + ((int128)tmp_q[2] * 14740502324586L) + ((int128)tmp_q[3] * 96801719312474L) - ((int128)tmp_q[4] * 95679293498708L) + ((int128)tmp_q[5] * 34277528353011L) + ((int128)tmp_q[6] * 48450626635358L) + ((int128)tmp_q[7] * 18157203451227L) + ((int128)tmp_q[8] * 71845021202114L) + ((int128)tmp_q[9] * 59664278664703L) - ((int128)tmp_q[10] * 63800393038753L);

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

