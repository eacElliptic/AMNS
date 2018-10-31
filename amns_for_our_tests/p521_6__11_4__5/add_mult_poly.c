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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10618670028315117875UL) + ((((uint64_t)op[1] * 2098793678053693488UL) + ((uint64_t)op[2] * 6132178058625961096UL) + ((uint64_t)op[3] * 18275314170071131833UL) + ((uint64_t)op[4] * 7015809469402054199UL) + ((uint64_t)op[5] * 2904525392357571422UL) + ((uint64_t)op[6] * 2007341490649990057UL) + ((uint64_t)op[7] * 14498876935781948838UL) + ((uint64_t)op[8] * 15508479383453377678UL) + ((uint64_t)op[9] * 2493304871107295038UL) + ((uint64_t)op[10] * 15284172028656315495UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 15284172028656315495UL) + ((uint64_t)op[1] * 10618670028315117875UL) + ((((uint64_t)op[2] * 2098793678053693488UL) + ((uint64_t)op[3] * 6132178058625961096UL) + ((uint64_t)op[4] * 18275314170071131833UL) + ((uint64_t)op[5] * 7015809469402054199UL) + ((uint64_t)op[6] * 2904525392357571422UL) + ((uint64_t)op[7] * 2007341490649990057UL) + ((uint64_t)op[8] * 14498876935781948838UL) + ((uint64_t)op[9] * 15508479383453377678UL) + ((uint64_t)op[10] * 2493304871107295038UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 2493304871107295038UL) + ((uint64_t)op[1] * 15284172028656315495UL) + ((uint64_t)op[2] * 10618670028315117875UL) + ((((uint64_t)op[3] * 2098793678053693488UL) + ((uint64_t)op[4] * 6132178058625961096UL) + ((uint64_t)op[5] * 18275314170071131833UL) + ((uint64_t)op[6] * 7015809469402054199UL) + ((uint64_t)op[7] * 2904525392357571422UL) + ((uint64_t)op[8] * 2007341490649990057UL) + ((uint64_t)op[9] * 14498876935781948838UL) + ((uint64_t)op[10] * 15508479383453377678UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 15508479383453377678UL) + ((uint64_t)op[1] * 2493304871107295038UL) + ((uint64_t)op[2] * 15284172028656315495UL) + ((uint64_t)op[3] * 10618670028315117875UL) + ((((uint64_t)op[4] * 2098793678053693488UL) + ((uint64_t)op[5] * 6132178058625961096UL) + ((uint64_t)op[6] * 18275314170071131833UL) + ((uint64_t)op[7] * 7015809469402054199UL) + ((uint64_t)op[8] * 2904525392357571422UL) + ((uint64_t)op[9] * 2007341490649990057UL) + ((uint64_t)op[10] * 14498876935781948838UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 14498876935781948838UL) + ((uint64_t)op[1] * 15508479383453377678UL) + ((uint64_t)op[2] * 2493304871107295038UL) + ((uint64_t)op[3] * 15284172028656315495UL) + ((uint64_t)op[4] * 10618670028315117875UL) + ((((uint64_t)op[5] * 2098793678053693488UL) + ((uint64_t)op[6] * 6132178058625961096UL) + ((uint64_t)op[7] * 18275314170071131833UL) + ((uint64_t)op[8] * 7015809469402054199UL) + ((uint64_t)op[9] * 2904525392357571422UL) + ((uint64_t)op[10] * 2007341490649990057UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 2007341490649990057UL) + ((uint64_t)op[1] * 14498876935781948838UL) + ((uint64_t)op[2] * 15508479383453377678UL) + ((uint64_t)op[3] * 2493304871107295038UL) + ((uint64_t)op[4] * 15284172028656315495UL) + ((uint64_t)op[5] * 10618670028315117875UL) + ((((uint64_t)op[6] * 2098793678053693488UL) + ((uint64_t)op[7] * 6132178058625961096UL) + ((uint64_t)op[8] * 18275314170071131833UL) + ((uint64_t)op[9] * 7015809469402054199UL) + ((uint64_t)op[10] * 2904525392357571422UL)) * 4);
	tmp_q[6] = ((uint64_t)op[0] * 2904525392357571422UL) + ((uint64_t)op[1] * 2007341490649990057UL) + ((uint64_t)op[2] * 14498876935781948838UL) + ((uint64_t)op[3] * 15508479383453377678UL) + ((uint64_t)op[4] * 2493304871107295038UL) + ((uint64_t)op[5] * 15284172028656315495UL) + ((uint64_t)op[6] * 10618670028315117875UL) + ((((uint64_t)op[7] * 2098793678053693488UL) + ((uint64_t)op[8] * 6132178058625961096UL) + ((uint64_t)op[9] * 18275314170071131833UL) + ((uint64_t)op[10] * 7015809469402054199UL)) * 4);
	tmp_q[7] = ((uint64_t)op[0] * 7015809469402054199UL) + ((uint64_t)op[1] * 2904525392357571422UL) + ((uint64_t)op[2] * 2007341490649990057UL) + ((uint64_t)op[3] * 14498876935781948838UL) + ((uint64_t)op[4] * 15508479383453377678UL) + ((uint64_t)op[5] * 2493304871107295038UL) + ((uint64_t)op[6] * 15284172028656315495UL) + ((uint64_t)op[7] * 10618670028315117875UL) + ((((uint64_t)op[8] * 2098793678053693488UL) + ((uint64_t)op[9] * 6132178058625961096UL) + ((uint64_t)op[10] * 18275314170071131833UL)) * 4);
	tmp_q[8] = ((uint64_t)op[0] * 18275314170071131833UL) + ((uint64_t)op[1] * 7015809469402054199UL) + ((uint64_t)op[2] * 2904525392357571422UL) + ((uint64_t)op[3] * 2007341490649990057UL) + ((uint64_t)op[4] * 14498876935781948838UL) + ((uint64_t)op[5] * 15508479383453377678UL) + ((uint64_t)op[6] * 2493304871107295038UL) + ((uint64_t)op[7] * 15284172028656315495UL) + ((uint64_t)op[8] * 10618670028315117875UL) + ((((uint64_t)op[9] * 2098793678053693488UL) + ((uint64_t)op[10] * 6132178058625961096UL)) * 4);
	tmp_q[9] = ((uint64_t)op[0] * 6132178058625961096UL) + ((uint64_t)op[1] * 18275314170071131833UL) + ((uint64_t)op[2] * 7015809469402054199UL) + ((uint64_t)op[3] * 2904525392357571422UL) + ((uint64_t)op[4] * 2007341490649990057UL) + ((uint64_t)op[5] * 14498876935781948838UL) + ((uint64_t)op[6] * 15508479383453377678UL) + ((uint64_t)op[7] * 2493304871107295038UL) + ((uint64_t)op[8] * 15284172028656315495UL) + ((uint64_t)op[9] * 10618670028315117875UL) + ((uint64_t)op[10] * 8395174712214773952UL);
	tmp_q[10] = ((uint64_t)op[0] * 2098793678053693488UL) + ((uint64_t)op[1] * 6132178058625961096UL) + ((uint64_t)op[2] * 18275314170071131833UL) + ((uint64_t)op[3] * 7015809469402054199UL) + ((uint64_t)op[4] * 2904525392357571422UL) + ((uint64_t)op[5] * 2007341490649990057UL) + ((uint64_t)op[6] * 14498876935781948838UL) + ((uint64_t)op[7] * 15508479383453377678UL) + ((uint64_t)op[8] * 2493304871107295038UL) + ((uint64_t)op[9] * 15284172028656315495UL) + ((uint64_t)op[10] * 10618670028315117875UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2416169928667L) + ((-((int128)tmp_q[1] * 77928104941131L) - ((int128)tmp_q[2] * 70796143036451L) + ((int128)tmp_q[3] * 40275477550758L) - ((int128)tmp_q[4] * 3371274605341L) - ((int128)tmp_q[5] * 12251470625231L) + ((int128)tmp_q[6] * 78886164007422L) - ((int128)tmp_q[7] * 56972919179867L) - ((int128)tmp_q[8] * 53420561126791L) + ((int128)tmp_q[9] * 9996521834287L) + ((int128)tmp_q[10] * 11331283053743L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 11331283053743L) - ((int128)tmp_q[1] * 2416169928667L) + ((-((int128)tmp_q[2] * 77928104941131L) - ((int128)tmp_q[3] * 70796143036451L) + ((int128)tmp_q[4] * 40275477550758L) - ((int128)tmp_q[5] * 3371274605341L) - ((int128)tmp_q[6] * 12251470625231L) + ((int128)tmp_q[7] * 78886164007422L) - ((int128)tmp_q[8] * 56972919179867L) - ((int128)tmp_q[9] * 53420561126791L) + ((int128)tmp_q[10] * 9996521834287L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 9996521834287L) + ((int128)tmp_q[1] * 11331283053743L) - ((int128)tmp_q[2] * 2416169928667L) + ((-((int128)tmp_q[3] * 77928104941131L) - ((int128)tmp_q[4] * 70796143036451L) + ((int128)tmp_q[5] * 40275477550758L) - ((int128)tmp_q[6] * 3371274605341L) - ((int128)tmp_q[7] * 12251470625231L) + ((int128)tmp_q[8] * 78886164007422L) - ((int128)tmp_q[9] * 56972919179867L) - ((int128)tmp_q[10] * 53420561126791L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 53420561126791L) + ((int128)tmp_q[1] * 9996521834287L) + ((int128)tmp_q[2] * 11331283053743L) - ((int128)tmp_q[3] * 2416169928667L) + ((-((int128)tmp_q[4] * 77928104941131L) - ((int128)tmp_q[5] * 70796143036451L) + ((int128)tmp_q[6] * 40275477550758L) - ((int128)tmp_q[7] * 3371274605341L) - ((int128)tmp_q[8] * 12251470625231L) + ((int128)tmp_q[9] * 78886164007422L) - ((int128)tmp_q[10] * 56972919179867L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 56972919179867L) - ((int128)tmp_q[1] * 53420561126791L) + ((int128)tmp_q[2] * 9996521834287L) + ((int128)tmp_q[3] * 11331283053743L) - ((int128)tmp_q[4] * 2416169928667L) + ((-((int128)tmp_q[5] * 77928104941131L) - ((int128)tmp_q[6] * 70796143036451L) + ((int128)tmp_q[7] * 40275477550758L) - ((int128)tmp_q[8] * 3371274605341L) - ((int128)tmp_q[9] * 12251470625231L) + ((int128)tmp_q[10] * 78886164007422L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 78886164007422L) - ((int128)tmp_q[1] * 56972919179867L) - ((int128)tmp_q[2] * 53420561126791L) + ((int128)tmp_q[3] * 9996521834287L) + ((int128)tmp_q[4] * 11331283053743L) - ((int128)tmp_q[5] * 2416169928667L) + ((-((int128)tmp_q[6] * 77928104941131L) - ((int128)tmp_q[7] * 70796143036451L) + ((int128)tmp_q[8] * 40275477550758L) - ((int128)tmp_q[9] * 3371274605341L) - ((int128)tmp_q[10] * 12251470625231L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 12251470625231L) + ((int128)tmp_q[1] * 78886164007422L) - ((int128)tmp_q[2] * 56972919179867L) - ((int128)tmp_q[3] * 53420561126791L) + ((int128)tmp_q[4] * 9996521834287L) + ((int128)tmp_q[5] * 11331283053743L) - ((int128)tmp_q[6] * 2416169928667L) + ((-((int128)tmp_q[7] * 77928104941131L) - ((int128)tmp_q[8] * 70796143036451L) + ((int128)tmp_q[9] * 40275477550758L) - ((int128)tmp_q[10] * 3371274605341L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 3371274605341L) - ((int128)tmp_q[1] * 12251470625231L) + ((int128)tmp_q[2] * 78886164007422L) - ((int128)tmp_q[3] * 56972919179867L) - ((int128)tmp_q[4] * 53420561126791L) + ((int128)tmp_q[5] * 9996521834287L) + ((int128)tmp_q[6] * 11331283053743L) - ((int128)tmp_q[7] * 2416169928667L) + ((-((int128)tmp_q[8] * 77928104941131L) - ((int128)tmp_q[9] * 70796143036451L) + ((int128)tmp_q[10] * 40275477550758L)) * 4);
	tmp_zero[8] = ((int128)tmp_q[0] * 40275477550758L) - ((int128)tmp_q[1] * 3371274605341L) - ((int128)tmp_q[2] * 12251470625231L) + ((int128)tmp_q[3] * 78886164007422L) - ((int128)tmp_q[4] * 56972919179867L) - ((int128)tmp_q[5] * 53420561126791L) + ((int128)tmp_q[6] * 9996521834287L) + ((int128)tmp_q[7] * 11331283053743L) - ((int128)tmp_q[8] * 2416169928667L) + ((-((int128)tmp_q[9] * 77928104941131L) - ((int128)tmp_q[10] * 70796143036451L)) * 4);
	tmp_zero[9] = -((int128)tmp_q[0] * 70796143036451L) + ((int128)tmp_q[1] * 40275477550758L) - ((int128)tmp_q[2] * 3371274605341L) - ((int128)tmp_q[3] * 12251470625231L) + ((int128)tmp_q[4] * 78886164007422L) - ((int128)tmp_q[5] * 56972919179867L) - ((int128)tmp_q[6] * 53420561126791L) + ((int128)tmp_q[7] * 9996521834287L) + ((int128)tmp_q[8] * 11331283053743L) - ((int128)tmp_q[9] * 2416169928667L) - ((int128)tmp_q[10] * 311712419764524L);
	tmp_zero[10] = -((int128)tmp_q[0] * 77928104941131L) - ((int128)tmp_q[1] * 70796143036451L) + ((int128)tmp_q[2] * 40275477550758L) - ((int128)tmp_q[3] * 3371274605341L) - ((int128)tmp_q[4] * 12251470625231L) + ((int128)tmp_q[5] * 78886164007422L) - ((int128)tmp_q[6] * 56972919179867L) - ((int128)tmp_q[7] * 53420561126791L) + ((int128)tmp_q[8] * 9996521834287L) + ((int128)tmp_q[9] * 11331283053743L) - ((int128)tmp_q[10] * 2416169928667L);

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

