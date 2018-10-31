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
	tmp_q[0] = ((uint64_t)op[0] * 9765468117783684351UL) + ((((uint64_t)op[1] * 1337262013104174518UL) + ((uint64_t)op[2] * 16345162209060407107UL) + ((uint64_t)op[3] * 4631442227946487937UL) + ((uint64_t)op[4] * 9734271055799886708UL) + ((uint64_t)op[5] * 4250727992788778437UL) + ((uint64_t)op[6] * 6982036278774956682UL) + ((uint64_t)op[7] * 6237789383306367664UL) + ((uint64_t)op[8] * 9044086548770055953UL) + ((uint64_t)op[9] * 13058553251053988183UL) + ((uint64_t)op[10] * 11898675230096574506UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 11898675230096574506UL) + ((uint64_t)op[1] * 9765468117783684351UL) + ((((uint64_t)op[2] * 1337262013104174518UL) + ((uint64_t)op[3] * 16345162209060407107UL) + ((uint64_t)op[4] * 4631442227946487937UL) + ((uint64_t)op[5] * 9734271055799886708UL) + ((uint64_t)op[6] * 4250727992788778437UL) + ((uint64_t)op[7] * 6982036278774956682UL) + ((uint64_t)op[8] * 6237789383306367664UL) + ((uint64_t)op[9] * 9044086548770055953UL) + ((uint64_t)op[10] * 13058553251053988183UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 13058553251053988183UL) + ((uint64_t)op[1] * 11898675230096574506UL) + ((uint64_t)op[2] * 9765468117783684351UL) + ((((uint64_t)op[3] * 1337262013104174518UL) + ((uint64_t)op[4] * 16345162209060407107UL) + ((uint64_t)op[5] * 4631442227946487937UL) + ((uint64_t)op[6] * 9734271055799886708UL) + ((uint64_t)op[7] * 4250727992788778437UL) + ((uint64_t)op[8] * 6982036278774956682UL) + ((uint64_t)op[9] * 6237789383306367664UL) + ((uint64_t)op[10] * 9044086548770055953UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 9044086548770055953UL) + ((uint64_t)op[1] * 13058553251053988183UL) + ((uint64_t)op[2] * 11898675230096574506UL) + ((uint64_t)op[3] * 9765468117783684351UL) + ((((uint64_t)op[4] * 1337262013104174518UL) + ((uint64_t)op[5] * 16345162209060407107UL) + ((uint64_t)op[6] * 4631442227946487937UL) + ((uint64_t)op[7] * 9734271055799886708UL) + ((uint64_t)op[8] * 4250727992788778437UL) + ((uint64_t)op[9] * 6982036278774956682UL) + ((uint64_t)op[10] * 6237789383306367664UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 6237789383306367664UL) + ((uint64_t)op[1] * 9044086548770055953UL) + ((uint64_t)op[2] * 13058553251053988183UL) + ((uint64_t)op[3] * 11898675230096574506UL) + ((uint64_t)op[4] * 9765468117783684351UL) + ((((uint64_t)op[5] * 1337262013104174518UL) + ((uint64_t)op[6] * 16345162209060407107UL) + ((uint64_t)op[7] * 4631442227946487937UL) + ((uint64_t)op[8] * 9734271055799886708UL) + ((uint64_t)op[9] * 4250727992788778437UL) + ((uint64_t)op[10] * 6982036278774956682UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 6982036278774956682UL) + ((uint64_t)op[1] * 6237789383306367664UL) + ((uint64_t)op[2] * 9044086548770055953UL) + ((uint64_t)op[3] * 13058553251053988183UL) + ((uint64_t)op[4] * 11898675230096574506UL) + ((uint64_t)op[5] * 9765468117783684351UL) + ((((uint64_t)op[6] * 1337262013104174518UL) + ((uint64_t)op[7] * 16345162209060407107UL) + ((uint64_t)op[8] * 4631442227946487937UL) + ((uint64_t)op[9] * 9734271055799886708UL) + ((uint64_t)op[10] * 4250727992788778437UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 4250727992788778437UL) + ((uint64_t)op[1] * 6982036278774956682UL) + ((uint64_t)op[2] * 6237789383306367664UL) + ((uint64_t)op[3] * 9044086548770055953UL) + ((uint64_t)op[4] * 13058553251053988183UL) + ((uint64_t)op[5] * 11898675230096574506UL) + ((uint64_t)op[6] * 9765468117783684351UL) + ((((uint64_t)op[7] * 1337262013104174518UL) + ((uint64_t)op[8] * 16345162209060407107UL) + ((uint64_t)op[9] * 4631442227946487937UL) + ((uint64_t)op[10] * 9734271055799886708UL)) * 8);
	tmp_q[7] = ((uint64_t)op[0] * 9734271055799886708UL) + ((uint64_t)op[1] * 4250727992788778437UL) + ((uint64_t)op[2] * 6982036278774956682UL) + ((uint64_t)op[3] * 6237789383306367664UL) + ((uint64_t)op[4] * 9044086548770055953UL) + ((uint64_t)op[5] * 13058553251053988183UL) + ((uint64_t)op[6] * 11898675230096574506UL) + ((uint64_t)op[7] * 9765468117783684351UL) + ((((uint64_t)op[8] * 1337262013104174518UL) + ((uint64_t)op[9] * 16345162209060407107UL) + ((uint64_t)op[10] * 4631442227946487937UL)) * 8);
	tmp_q[8] = ((uint64_t)op[0] * 4631442227946487937UL) + ((uint64_t)op[1] * 9734271055799886708UL) + ((uint64_t)op[2] * 4250727992788778437UL) + ((uint64_t)op[3] * 6982036278774956682UL) + ((uint64_t)op[4] * 6237789383306367664UL) + ((uint64_t)op[5] * 9044086548770055953UL) + ((uint64_t)op[6] * 13058553251053988183UL) + ((uint64_t)op[7] * 11898675230096574506UL) + ((uint64_t)op[8] * 9765468117783684351UL) + ((((uint64_t)op[9] * 1337262013104174518UL) + ((uint64_t)op[10] * 16345162209060407107UL)) * 8);
	tmp_q[9] = ((uint64_t)op[0] * 16345162209060407107UL) + ((uint64_t)op[1] * 4631442227946487937UL) + ((uint64_t)op[2] * 9734271055799886708UL) + ((uint64_t)op[3] * 4250727992788778437UL) + ((uint64_t)op[4] * 6982036278774956682UL) + ((uint64_t)op[5] * 6237789383306367664UL) + ((uint64_t)op[6] * 9044086548770055953UL) + ((uint64_t)op[7] * 13058553251053988183UL) + ((uint64_t)op[8] * 11898675230096574506UL) + ((uint64_t)op[9] * 9765468117783684351UL) + ((uint64_t)op[10] * 10698096104833396144UL);
	tmp_q[10] = ((uint64_t)op[0] * 1337262013104174518UL) + ((uint64_t)op[1] * 16345162209060407107UL) + ((uint64_t)op[2] * 4631442227946487937UL) + ((uint64_t)op[3] * 9734271055799886708UL) + ((uint64_t)op[4] * 4250727992788778437UL) + ((uint64_t)op[5] * 6982036278774956682UL) + ((uint64_t)op[6] * 6237789383306367664UL) + ((uint64_t)op[7] * 9044086548770055953UL) + ((uint64_t)op[8] * 13058553251053988183UL) + ((uint64_t)op[9] * 11898675230096574506UL) + ((uint64_t)op[10] * 9765468117783684351UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 89424227533801L) + ((-((int128)tmp_q[1] * 58237125464592L) + ((int128)tmp_q[2] * 89843554917138L) - ((int128)tmp_q[3] * 15102264754699L) + ((int128)tmp_q[4] * 82407524581861L) + ((int128)tmp_q[5] * 19276695739137L) - ((int128)tmp_q[6] * 37975213051750L) - ((int128)tmp_q[7] * 96465685441575L) + ((int128)tmp_q[8] * 33310968270757L) - ((int128)tmp_q[9] * 38586034257613L) - ((int128)tmp_q[10] * 57655347743214L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 57655347743214L) + ((int128)tmp_q[1] * 89424227533801L) + ((-((int128)tmp_q[2] * 58237125464592L) + ((int128)tmp_q[3] * 89843554917138L) - ((int128)tmp_q[4] * 15102264754699L) + ((int128)tmp_q[5] * 82407524581861L) + ((int128)tmp_q[6] * 19276695739137L) - ((int128)tmp_q[7] * 37975213051750L) - ((int128)tmp_q[8] * 96465685441575L) + ((int128)tmp_q[9] * 33310968270757L) - ((int128)tmp_q[10] * 38586034257613L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 38586034257613L) - ((int128)tmp_q[1] * 57655347743214L) + ((int128)tmp_q[2] * 89424227533801L) + ((-((int128)tmp_q[3] * 58237125464592L) + ((int128)tmp_q[4] * 89843554917138L) - ((int128)tmp_q[5] * 15102264754699L) + ((int128)tmp_q[6] * 82407524581861L) + ((int128)tmp_q[7] * 19276695739137L) - ((int128)tmp_q[8] * 37975213051750L) - ((int128)tmp_q[9] * 96465685441575L) + ((int128)tmp_q[10] * 33310968270757L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 33310968270757L) - ((int128)tmp_q[1] * 38586034257613L) - ((int128)tmp_q[2] * 57655347743214L) + ((int128)tmp_q[3] * 89424227533801L) + ((-((int128)tmp_q[4] * 58237125464592L) + ((int128)tmp_q[5] * 89843554917138L) - ((int128)tmp_q[6] * 15102264754699L) + ((int128)tmp_q[7] * 82407524581861L) + ((int128)tmp_q[8] * 19276695739137L) - ((int128)tmp_q[9] * 37975213051750L) - ((int128)tmp_q[10] * 96465685441575L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 96465685441575L) + ((int128)tmp_q[1] * 33310968270757L) - ((int128)tmp_q[2] * 38586034257613L) - ((int128)tmp_q[3] * 57655347743214L) + ((int128)tmp_q[4] * 89424227533801L) + ((-((int128)tmp_q[5] * 58237125464592L) + ((int128)tmp_q[6] * 89843554917138L) - ((int128)tmp_q[7] * 15102264754699L) + ((int128)tmp_q[8] * 82407524581861L) + ((int128)tmp_q[9] * 19276695739137L) - ((int128)tmp_q[10] * 37975213051750L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 37975213051750L) - ((int128)tmp_q[1] * 96465685441575L) + ((int128)tmp_q[2] * 33310968270757L) - ((int128)tmp_q[3] * 38586034257613L) - ((int128)tmp_q[4] * 57655347743214L) + ((int128)tmp_q[5] * 89424227533801L) + ((-((int128)tmp_q[6] * 58237125464592L) + ((int128)tmp_q[7] * 89843554917138L) - ((int128)tmp_q[8] * 15102264754699L) + ((int128)tmp_q[9] * 82407524581861L) + ((int128)tmp_q[10] * 19276695739137L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 19276695739137L) - ((int128)tmp_q[1] * 37975213051750L) - ((int128)tmp_q[2] * 96465685441575L) + ((int128)tmp_q[3] * 33310968270757L) - ((int128)tmp_q[4] * 38586034257613L) - ((int128)tmp_q[5] * 57655347743214L) + ((int128)tmp_q[6] * 89424227533801L) + ((-((int128)tmp_q[7] * 58237125464592L) + ((int128)tmp_q[8] * 89843554917138L) - ((int128)tmp_q[9] * 15102264754699L) + ((int128)tmp_q[10] * 82407524581861L)) * 8);
	tmp_zero[7] = ((int128)tmp_q[0] * 82407524581861L) + ((int128)tmp_q[1] * 19276695739137L) - ((int128)tmp_q[2] * 37975213051750L) - ((int128)tmp_q[3] * 96465685441575L) + ((int128)tmp_q[4] * 33310968270757L) - ((int128)tmp_q[5] * 38586034257613L) - ((int128)tmp_q[6] * 57655347743214L) + ((int128)tmp_q[7] * 89424227533801L) + ((-((int128)tmp_q[8] * 58237125464592L) + ((int128)tmp_q[9] * 89843554917138L) - ((int128)tmp_q[10] * 15102264754699L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 15102264754699L) + ((int128)tmp_q[1] * 82407524581861L) + ((int128)tmp_q[2] * 19276695739137L) - ((int128)tmp_q[3] * 37975213051750L) - ((int128)tmp_q[4] * 96465685441575L) + ((int128)tmp_q[5] * 33310968270757L) - ((int128)tmp_q[6] * 38586034257613L) - ((int128)tmp_q[7] * 57655347743214L) + ((int128)tmp_q[8] * 89424227533801L) + ((-((int128)tmp_q[9] * 58237125464592L) + ((int128)tmp_q[10] * 89843554917138L)) * 8);
	tmp_zero[9] = ((int128)tmp_q[0] * 89843554917138L) - ((int128)tmp_q[1] * 15102264754699L) + ((int128)tmp_q[2] * 82407524581861L) + ((int128)tmp_q[3] * 19276695739137L) - ((int128)tmp_q[4] * 37975213051750L) - ((int128)tmp_q[5] * 96465685441575L) + ((int128)tmp_q[6] * 33310968270757L) - ((int128)tmp_q[7] * 38586034257613L) - ((int128)tmp_q[8] * 57655347743214L) + ((int128)tmp_q[9] * 89424227533801L) - ((int128)tmp_q[10] * 465897003716736L);
	tmp_zero[10] = -((int128)tmp_q[0] * 58237125464592L) + ((int128)tmp_q[1] * 89843554917138L) - ((int128)tmp_q[2] * 15102264754699L) + ((int128)tmp_q[3] * 82407524581861L) + ((int128)tmp_q[4] * 19276695739137L) - ((int128)tmp_q[5] * 37975213051750L) - ((int128)tmp_q[6] * 96465685441575L) + ((int128)tmp_q[7] * 33310968270757L) - ((int128)tmp_q[8] * 38586034257613L) - ((int128)tmp_q[9] * 57655347743214L) + ((int128)tmp_q[10] * 89424227533801L);

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

