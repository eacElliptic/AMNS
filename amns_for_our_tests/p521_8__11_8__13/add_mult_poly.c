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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 443491073024246939UL) + ((((uint64_t)op[1] * 12499470936714128278UL) + ((uint64_t)op[2] * 5039109474470487823UL) + ((uint64_t)op[3] * 18186642624921087499UL) + ((uint64_t)op[4] * 10149515908848611353UL) + ((uint64_t)op[5] * 10573049681814167774UL) + ((uint64_t)op[6] * 14523886747070143831UL) + ((uint64_t)op[7] * 16163226236915093516UL) + ((uint64_t)op[8] * 10940522967764167559UL) + ((uint64_t)op[9] * 12112436242021193608UL) + ((uint64_t)op[10] * 13279177399391500368UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 13279177399391500368UL) + ((uint64_t)op[1] * 443491073024246939UL) + ((((uint64_t)op[2] * 12499470936714128278UL) + ((uint64_t)op[3] * 5039109474470487823UL) + ((uint64_t)op[4] * 18186642624921087499UL) + ((uint64_t)op[5] * 10149515908848611353UL) + ((uint64_t)op[6] * 10573049681814167774UL) + ((uint64_t)op[7] * 14523886747070143831UL) + ((uint64_t)op[8] * 16163226236915093516UL) + ((uint64_t)op[9] * 10940522967764167559UL) + ((uint64_t)op[10] * 12112436242021193608UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 12112436242021193608UL) + ((uint64_t)op[1] * 13279177399391500368UL) + ((uint64_t)op[2] * 443491073024246939UL) + ((((uint64_t)op[3] * 12499470936714128278UL) + ((uint64_t)op[4] * 5039109474470487823UL) + ((uint64_t)op[5] * 18186642624921087499UL) + ((uint64_t)op[6] * 10149515908848611353UL) + ((uint64_t)op[7] * 10573049681814167774UL) + ((uint64_t)op[8] * 14523886747070143831UL) + ((uint64_t)op[9] * 16163226236915093516UL) + ((uint64_t)op[10] * 10940522967764167559UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 10940522967764167559UL) + ((uint64_t)op[1] * 12112436242021193608UL) + ((uint64_t)op[2] * 13279177399391500368UL) + ((uint64_t)op[3] * 443491073024246939UL) + ((((uint64_t)op[4] * 12499470936714128278UL) + ((uint64_t)op[5] * 5039109474470487823UL) + ((uint64_t)op[6] * 18186642624921087499UL) + ((uint64_t)op[7] * 10149515908848611353UL) + ((uint64_t)op[8] * 10573049681814167774UL) + ((uint64_t)op[9] * 14523886747070143831UL) + ((uint64_t)op[10] * 16163226236915093516UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 16163226236915093516UL) + ((uint64_t)op[1] * 10940522967764167559UL) + ((uint64_t)op[2] * 12112436242021193608UL) + ((uint64_t)op[3] * 13279177399391500368UL) + ((uint64_t)op[4] * 443491073024246939UL) + ((((uint64_t)op[5] * 12499470936714128278UL) + ((uint64_t)op[6] * 5039109474470487823UL) + ((uint64_t)op[7] * 18186642624921087499UL) + ((uint64_t)op[8] * 10149515908848611353UL) + ((uint64_t)op[9] * 10573049681814167774UL) + ((uint64_t)op[10] * 14523886747070143831UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 14523886747070143831UL) + ((uint64_t)op[1] * 16163226236915093516UL) + ((uint64_t)op[2] * 10940522967764167559UL) + ((uint64_t)op[3] * 12112436242021193608UL) + ((uint64_t)op[4] * 13279177399391500368UL) + ((uint64_t)op[5] * 443491073024246939UL) + ((((uint64_t)op[6] * 12499470936714128278UL) + ((uint64_t)op[7] * 5039109474470487823UL) + ((uint64_t)op[8] * 18186642624921087499UL) + ((uint64_t)op[9] * 10149515908848611353UL) + ((uint64_t)op[10] * 10573049681814167774UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 10573049681814167774UL) + ((uint64_t)op[1] * 14523886747070143831UL) + ((uint64_t)op[2] * 16163226236915093516UL) + ((uint64_t)op[3] * 10940522967764167559UL) + ((uint64_t)op[4] * 12112436242021193608UL) + ((uint64_t)op[5] * 13279177399391500368UL) + ((uint64_t)op[6] * 443491073024246939UL) + ((((uint64_t)op[7] * 12499470936714128278UL) + ((uint64_t)op[8] * 5039109474470487823UL) + ((uint64_t)op[9] * 18186642624921087499UL) + ((uint64_t)op[10] * 10149515908848611353UL)) * 8);
	tmp_q[7] = ((uint64_t)op[0] * 10149515908848611353UL) + ((uint64_t)op[1] * 10573049681814167774UL) + ((uint64_t)op[2] * 14523886747070143831UL) + ((uint64_t)op[3] * 16163226236915093516UL) + ((uint64_t)op[4] * 10940522967764167559UL) + ((uint64_t)op[5] * 12112436242021193608UL) + ((uint64_t)op[6] * 13279177399391500368UL) + ((uint64_t)op[7] * 443491073024246939UL) + ((((uint64_t)op[8] * 12499470936714128278UL) + ((uint64_t)op[9] * 5039109474470487823UL) + ((uint64_t)op[10] * 18186642624921087499UL)) * 8);
	tmp_q[8] = ((uint64_t)op[0] * 18186642624921087499UL) + ((uint64_t)op[1] * 10149515908848611353UL) + ((uint64_t)op[2] * 10573049681814167774UL) + ((uint64_t)op[3] * 14523886747070143831UL) + ((uint64_t)op[4] * 16163226236915093516UL) + ((uint64_t)op[5] * 10940522967764167559UL) + ((uint64_t)op[6] * 12112436242021193608UL) + ((uint64_t)op[7] * 13279177399391500368UL) + ((uint64_t)op[8] * 443491073024246939UL) + ((((uint64_t)op[9] * 12499470936714128278UL) + ((uint64_t)op[10] * 5039109474470487823UL)) * 8);
	tmp_q[9] = ((uint64_t)op[0] * 5039109474470487823UL) + ((uint64_t)op[1] * 18186642624921087499UL) + ((uint64_t)op[2] * 10149515908848611353UL) + ((uint64_t)op[3] * 10573049681814167774UL) + ((uint64_t)op[4] * 14523886747070143831UL) + ((uint64_t)op[5] * 16163226236915093516UL) + ((uint64_t)op[6] * 10940522967764167559UL) + ((uint64_t)op[7] * 12112436242021193608UL) + ((uint64_t)op[8] * 13279177399391500368UL) + ((uint64_t)op[9] * 443491073024246939UL) + ((uint64_t)op[10] * 7762047125165268144UL);
	tmp_q[10] = ((uint64_t)op[0] * 12499470936714128278UL) + ((uint64_t)op[1] * 5039109474470487823UL) + ((uint64_t)op[2] * 18186642624921087499UL) + ((uint64_t)op[3] * 10149515908848611353UL) + ((uint64_t)op[4] * 10573049681814167774UL) + ((uint64_t)op[5] * 14523886747070143831UL) + ((uint64_t)op[6] * 16163226236915093516UL) + ((uint64_t)op[7] * 10940522967764167559UL) + ((uint64_t)op[8] * 12112436242021193608UL) + ((uint64_t)op[9] * 13279177399391500368UL) + ((uint64_t)op[10] * 443491073024246939UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 48967629741051L) + ((-((int128)tmp_q[1] * 4805457113187L) + ((int128)tmp_q[2] * 36708160099130L) - ((int128)tmp_q[3] * 725490341995L) - ((int128)tmp_q[4] * 80710406256863L) + ((int128)tmp_q[5] * 19667594123531L) - ((int128)tmp_q[6] * 18833164932169L) - ((int128)tmp_q[7] * 93909315923548L) - ((int128)tmp_q[8] * 74695018882657L) - ((int128)tmp_q[9] * 11494169998696L) + ((int128)tmp_q[10] * 28542094856744L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 28542094856744L) - ((int128)tmp_q[1] * 48967629741051L) + ((-((int128)tmp_q[2] * 4805457113187L) + ((int128)tmp_q[3] * 36708160099130L) - ((int128)tmp_q[4] * 725490341995L) - ((int128)tmp_q[5] * 80710406256863L) + ((int128)tmp_q[6] * 19667594123531L) - ((int128)tmp_q[7] * 18833164932169L) - ((int128)tmp_q[8] * 93909315923548L) - ((int128)tmp_q[9] * 74695018882657L) - ((int128)tmp_q[10] * 11494169998696L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 11494169998696L) + ((int128)tmp_q[1] * 28542094856744L) - ((int128)tmp_q[2] * 48967629741051L) + ((-((int128)tmp_q[3] * 4805457113187L) + ((int128)tmp_q[4] * 36708160099130L) - ((int128)tmp_q[5] * 725490341995L) - ((int128)tmp_q[6] * 80710406256863L) + ((int128)tmp_q[7] * 19667594123531L) - ((int128)tmp_q[8] * 18833164932169L) - ((int128)tmp_q[9] * 93909315923548L) - ((int128)tmp_q[10] * 74695018882657L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 74695018882657L) - ((int128)tmp_q[1] * 11494169998696L) + ((int128)tmp_q[2] * 28542094856744L) - ((int128)tmp_q[3] * 48967629741051L) + ((-((int128)tmp_q[4] * 4805457113187L) + ((int128)tmp_q[5] * 36708160099130L) - ((int128)tmp_q[6] * 725490341995L) - ((int128)tmp_q[7] * 80710406256863L) + ((int128)tmp_q[8] * 19667594123531L) - ((int128)tmp_q[9] * 18833164932169L) - ((int128)tmp_q[10] * 93909315923548L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 93909315923548L) - ((int128)tmp_q[1] * 74695018882657L) - ((int128)tmp_q[2] * 11494169998696L) + ((int128)tmp_q[3] * 28542094856744L) - ((int128)tmp_q[4] * 48967629741051L) + ((-((int128)tmp_q[5] * 4805457113187L) + ((int128)tmp_q[6] * 36708160099130L) - ((int128)tmp_q[7] * 725490341995L) - ((int128)tmp_q[8] * 80710406256863L) + ((int128)tmp_q[9] * 19667594123531L) - ((int128)tmp_q[10] * 18833164932169L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 18833164932169L) - ((int128)tmp_q[1] * 93909315923548L) - ((int128)tmp_q[2] * 74695018882657L) - ((int128)tmp_q[3] * 11494169998696L) + ((int128)tmp_q[4] * 28542094856744L) - ((int128)tmp_q[5] * 48967629741051L) + ((-((int128)tmp_q[6] * 4805457113187L) + ((int128)tmp_q[7] * 36708160099130L) - ((int128)tmp_q[8] * 725490341995L) - ((int128)tmp_q[9] * 80710406256863L) + ((int128)tmp_q[10] * 19667594123531L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 19667594123531L) - ((int128)tmp_q[1] * 18833164932169L) - ((int128)tmp_q[2] * 93909315923548L) - ((int128)tmp_q[3] * 74695018882657L) - ((int128)tmp_q[4] * 11494169998696L) + ((int128)tmp_q[5] * 28542094856744L) - ((int128)tmp_q[6] * 48967629741051L) + ((-((int128)tmp_q[7] * 4805457113187L) + ((int128)tmp_q[8] * 36708160099130L) - ((int128)tmp_q[9] * 725490341995L) - ((int128)tmp_q[10] * 80710406256863L)) * 8);
	tmp_zero[7] = -((int128)tmp_q[0] * 80710406256863L) + ((int128)tmp_q[1] * 19667594123531L) - ((int128)tmp_q[2] * 18833164932169L) - ((int128)tmp_q[3] * 93909315923548L) - ((int128)tmp_q[4] * 74695018882657L) - ((int128)tmp_q[5] * 11494169998696L) + ((int128)tmp_q[6] * 28542094856744L) - ((int128)tmp_q[7] * 48967629741051L) + ((-((int128)tmp_q[8] * 4805457113187L) + ((int128)tmp_q[9] * 36708160099130L) - ((int128)tmp_q[10] * 725490341995L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 725490341995L) - ((int128)tmp_q[1] * 80710406256863L) + ((int128)tmp_q[2] * 19667594123531L) - ((int128)tmp_q[3] * 18833164932169L) - ((int128)tmp_q[4] * 93909315923548L) - ((int128)tmp_q[5] * 74695018882657L) - ((int128)tmp_q[6] * 11494169998696L) + ((int128)tmp_q[7] * 28542094856744L) - ((int128)tmp_q[8] * 48967629741051L) + ((-((int128)tmp_q[9] * 4805457113187L) + ((int128)tmp_q[10] * 36708160099130L)) * 8);
	tmp_zero[9] = ((int128)tmp_q[0] * 36708160099130L) - ((int128)tmp_q[1] * 725490341995L) - ((int128)tmp_q[2] * 80710406256863L) + ((int128)tmp_q[3] * 19667594123531L) - ((int128)tmp_q[4] * 18833164932169L) - ((int128)tmp_q[5] * 93909315923548L) - ((int128)tmp_q[6] * 74695018882657L) - ((int128)tmp_q[7] * 11494169998696L) + ((int128)tmp_q[8] * 28542094856744L) - ((int128)tmp_q[9] * 48967629741051L) - ((int128)tmp_q[10] * 38443656905496L);
	tmp_zero[10] = -((int128)tmp_q[0] * 4805457113187L) + ((int128)tmp_q[1] * 36708160099130L) - ((int128)tmp_q[2] * 725490341995L) - ((int128)tmp_q[3] * 80710406256863L) + ((int128)tmp_q[4] * 19667594123531L) - ((int128)tmp_q[5] * 18833164932169L) - ((int128)tmp_q[6] * 93909315923548L) - ((int128)tmp_q[7] * 74695018882657L) - ((int128)tmp_q[8] * 11494169998696L) + ((int128)tmp_q[9] * 28542094856744L) - ((int128)tmp_q[10] * 48967629741051L);

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

