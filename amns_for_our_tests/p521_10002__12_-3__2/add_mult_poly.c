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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 3);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 3);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1837462468107062542UL) + ((((uint64_t)op[1] * 7215755333360978146UL) + ((uint64_t)op[2] * 13157328056901798191UL) + ((uint64_t)op[3] * 16186874219411306488UL) + ((uint64_t)op[4] * 5488426770571784906UL) + ((uint64_t)op[5] * 3528317179278578294UL) + ((uint64_t)op[6] * 536350285443149889UL) + ((uint64_t)op[7] * 17438000240645488493UL) + ((uint64_t)op[8] * 8524667427531356093UL) + ((uint64_t)op[9] * 7513853617723120944UL) + ((uint64_t)op[10] * 15945754284570131111UL) + ((uint64_t)op[11] * 6306262455555777438UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 6306262455555777438UL) + ((uint64_t)op[1] * 1837462468107062542UL) + ((((uint64_t)op[2] * 7215755333360978146UL) + ((uint64_t)op[3] * 13157328056901798191UL) + ((uint64_t)op[4] * 16186874219411306488UL) + ((uint64_t)op[5] * 5488426770571784906UL) + ((uint64_t)op[6] * 3528317179278578294UL) + ((uint64_t)op[7] * 536350285443149889UL) + ((uint64_t)op[8] * 17438000240645488493UL) + ((uint64_t)op[9] * 8524667427531356093UL) + ((uint64_t)op[10] * 7513853617723120944UL) + ((uint64_t)op[11] * 15945754284570131111UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 15945754284570131111UL) + ((uint64_t)op[1] * 6306262455555777438UL) + ((uint64_t)op[2] * 1837462468107062542UL) + ((((uint64_t)op[3] * 7215755333360978146UL) + ((uint64_t)op[4] * 13157328056901798191UL) + ((uint64_t)op[5] * 16186874219411306488UL) + ((uint64_t)op[6] * 5488426770571784906UL) + ((uint64_t)op[7] * 3528317179278578294UL) + ((uint64_t)op[8] * 536350285443149889UL) + ((uint64_t)op[9] * 17438000240645488493UL) + ((uint64_t)op[10] * 8524667427531356093UL) + ((uint64_t)op[11] * 7513853617723120944UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7513853617723120944UL) + ((uint64_t)op[1] * 15945754284570131111UL) + ((uint64_t)op[2] * 6306262455555777438UL) + ((uint64_t)op[3] * 1837462468107062542UL) + ((((uint64_t)op[4] * 7215755333360978146UL) + ((uint64_t)op[5] * 13157328056901798191UL) + ((uint64_t)op[6] * 16186874219411306488UL) + ((uint64_t)op[7] * 5488426770571784906UL) + ((uint64_t)op[8] * 3528317179278578294UL) + ((uint64_t)op[9] * 536350285443149889UL) + ((uint64_t)op[10] * 17438000240645488493UL) + ((uint64_t)op[11] * 8524667427531356093UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 8524667427531356093UL) + ((uint64_t)op[1] * 7513853617723120944UL) + ((uint64_t)op[2] * 15945754284570131111UL) + ((uint64_t)op[3] * 6306262455555777438UL) + ((uint64_t)op[4] * 1837462468107062542UL) + ((((uint64_t)op[5] * 7215755333360978146UL) + ((uint64_t)op[6] * 13157328056901798191UL) + ((uint64_t)op[7] * 16186874219411306488UL) + ((uint64_t)op[8] * 5488426770571784906UL) + ((uint64_t)op[9] * 3528317179278578294UL) + ((uint64_t)op[10] * 536350285443149889UL) + ((uint64_t)op[11] * 17438000240645488493UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 17438000240645488493UL) + ((uint64_t)op[1] * 8524667427531356093UL) + ((uint64_t)op[2] * 7513853617723120944UL) + ((uint64_t)op[3] * 15945754284570131111UL) + ((uint64_t)op[4] * 6306262455555777438UL) + ((uint64_t)op[5] * 1837462468107062542UL) + ((((uint64_t)op[6] * 7215755333360978146UL) + ((uint64_t)op[7] * 13157328056901798191UL) + ((uint64_t)op[8] * 16186874219411306488UL) + ((uint64_t)op[9] * 5488426770571784906UL) + ((uint64_t)op[10] * 3528317179278578294UL) + ((uint64_t)op[11] * 536350285443149889UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 536350285443149889UL) + ((uint64_t)op[1] * 17438000240645488493UL) + ((uint64_t)op[2] * 8524667427531356093UL) + ((uint64_t)op[3] * 7513853617723120944UL) + ((uint64_t)op[4] * 15945754284570131111UL) + ((uint64_t)op[5] * 6306262455555777438UL) + ((uint64_t)op[6] * 1837462468107062542UL) + ((((uint64_t)op[7] * 7215755333360978146UL) + ((uint64_t)op[8] * 13157328056901798191UL) + ((uint64_t)op[9] * 16186874219411306488UL) + ((uint64_t)op[10] * 5488426770571784906UL) + ((uint64_t)op[11] * 3528317179278578294UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 3528317179278578294UL) + ((uint64_t)op[1] * 536350285443149889UL) + ((uint64_t)op[2] * 17438000240645488493UL) + ((uint64_t)op[3] * 8524667427531356093UL) + ((uint64_t)op[4] * 7513853617723120944UL) + ((uint64_t)op[5] * 15945754284570131111UL) + ((uint64_t)op[6] * 6306262455555777438UL) + ((uint64_t)op[7] * 1837462468107062542UL) + ((((uint64_t)op[8] * 7215755333360978146UL) + ((uint64_t)op[9] * 13157328056901798191UL) + ((uint64_t)op[10] * 16186874219411306488UL) + ((uint64_t)op[11] * 5488426770571784906UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 5488426770571784906UL) + ((uint64_t)op[1] * 3528317179278578294UL) + ((uint64_t)op[2] * 536350285443149889UL) + ((uint64_t)op[3] * 17438000240645488493UL) + ((uint64_t)op[4] * 8524667427531356093UL) + ((uint64_t)op[5] * 7513853617723120944UL) + ((uint64_t)op[6] * 15945754284570131111UL) + ((uint64_t)op[7] * 6306262455555777438UL) + ((uint64_t)op[8] * 1837462468107062542UL) + ((((uint64_t)op[9] * 7215755333360978146UL) + ((uint64_t)op[10] * 13157328056901798191UL) + ((uint64_t)op[11] * 16186874219411306488UL)) * 18446744073709551613);
	tmp_q[9] = ((uint64_t)op[0] * 16186874219411306488UL) + ((uint64_t)op[1] * 5488426770571784906UL) + ((uint64_t)op[2] * 3528317179278578294UL) + ((uint64_t)op[3] * 536350285443149889UL) + ((uint64_t)op[4] * 17438000240645488493UL) + ((uint64_t)op[5] * 8524667427531356093UL) + ((uint64_t)op[6] * 7513853617723120944UL) + ((uint64_t)op[7] * 15945754284570131111UL) + ((uint64_t)op[8] * 6306262455555777438UL) + ((uint64_t)op[9] * 1837462468107062542UL) + ((((uint64_t)op[10] * 7215755333360978146UL) + ((uint64_t)op[11] * 13157328056901798191UL)) * 18446744073709551613);
	tmp_q[10] = ((uint64_t)op[0] * 13157328056901798191UL) + ((uint64_t)op[1] * 16186874219411306488UL) + ((uint64_t)op[2] * 5488426770571784906UL) + ((uint64_t)op[3] * 3528317179278578294UL) + ((uint64_t)op[4] * 536350285443149889UL) + ((uint64_t)op[5] * 17438000240645488493UL) + ((uint64_t)op[6] * 8524667427531356093UL) + ((uint64_t)op[7] * 7513853617723120944UL) + ((uint64_t)op[8] * 15945754284570131111UL) + ((uint64_t)op[9] * 6306262455555777438UL) + ((uint64_t)op[10] * 1837462468107062542UL) + ((uint64_t)op[11] * 15246222147336168794UL);
	tmp_q[11] = ((uint64_t)op[0] * 7215755333360978146UL) + ((uint64_t)op[1] * 13157328056901798191UL) + ((uint64_t)op[2] * 16186874219411306488UL) + ((uint64_t)op[3] * 5488426770571784906UL) + ((uint64_t)op[4] * 3528317179278578294UL) + ((uint64_t)op[5] * 536350285443149889UL) + ((uint64_t)op[6] * 17438000240645488493UL) + ((uint64_t)op[7] * 8524667427531356093UL) + ((uint64_t)op[8] * 7513853617723120944UL) + ((uint64_t)op[9] * 15945754284570131111UL) + ((uint64_t)op[10] * 6306262455555777438UL) + ((uint64_t)op[11] * 1837462468107062542UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 305763315507L) - ((((int128)tmp_q[1] * 251994405312L) - ((int128)tmp_q[2] * 4415422271200L) + ((int128)tmp_q[3] * 2988802878691L) - ((int128)tmp_q[4] * 3780952059490L) + ((int128)tmp_q[5] * 4081351147662L) + ((int128)tmp_q[6] * 6644387537682L) + ((int128)tmp_q[7] * 1436766058395L) - ((int128)tmp_q[8] * 161419853902L) + ((int128)tmp_q[9] * 6729544726213L) + ((int128)tmp_q[10] * 5084286511307L) + ((int128)tmp_q[11] * 1277249840956L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 1277249840956L) + ((int128)tmp_q[1] * 305763315507L) - ((((int128)tmp_q[2] * 251994405312L) - ((int128)tmp_q[3] * 4415422271200L) + ((int128)tmp_q[4] * 2988802878691L) - ((int128)tmp_q[5] * 3780952059490L) + ((int128)tmp_q[6] * 4081351147662L) + ((int128)tmp_q[7] * 6644387537682L) + ((int128)tmp_q[8] * 1436766058395L) - ((int128)tmp_q[9] * 161419853902L) + ((int128)tmp_q[10] * 6729544726213L) + ((int128)tmp_q[11] * 5084286511307L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 5084286511307L) + ((int128)tmp_q[1] * 1277249840956L) + ((int128)tmp_q[2] * 305763315507L) - ((((int128)tmp_q[3] * 251994405312L) - ((int128)tmp_q[4] * 4415422271200L) + ((int128)tmp_q[5] * 2988802878691L) - ((int128)tmp_q[6] * 3780952059490L) + ((int128)tmp_q[7] * 4081351147662L) + ((int128)tmp_q[8] * 6644387537682L) + ((int128)tmp_q[9] * 1436766058395L) - ((int128)tmp_q[10] * 161419853902L) + ((int128)tmp_q[11] * 6729544726213L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 6729544726213L) + ((int128)tmp_q[1] * 5084286511307L) + ((int128)tmp_q[2] * 1277249840956L) + ((int128)tmp_q[3] * 305763315507L) - ((((int128)tmp_q[4] * 251994405312L) - ((int128)tmp_q[5] * 4415422271200L) + ((int128)tmp_q[6] * 2988802878691L) - ((int128)tmp_q[7] * 3780952059490L) + ((int128)tmp_q[8] * 4081351147662L) + ((int128)tmp_q[9] * 6644387537682L) + ((int128)tmp_q[10] * 1436766058395L) - ((int128)tmp_q[11] * 161419853902L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 161419853902L) + ((int128)tmp_q[1] * 6729544726213L) + ((int128)tmp_q[2] * 5084286511307L) + ((int128)tmp_q[3] * 1277249840956L) + ((int128)tmp_q[4] * 305763315507L) - ((((int128)tmp_q[5] * 251994405312L) - ((int128)tmp_q[6] * 4415422271200L) + ((int128)tmp_q[7] * 2988802878691L) - ((int128)tmp_q[8] * 3780952059490L) + ((int128)tmp_q[9] * 4081351147662L) + ((int128)tmp_q[10] * 6644387537682L) + ((int128)tmp_q[11] * 1436766058395L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 1436766058395L) - ((int128)tmp_q[1] * 161419853902L) + ((int128)tmp_q[2] * 6729544726213L) + ((int128)tmp_q[3] * 5084286511307L) + ((int128)tmp_q[4] * 1277249840956L) + ((int128)tmp_q[5] * 305763315507L) - ((((int128)tmp_q[6] * 251994405312L) - ((int128)tmp_q[7] * 4415422271200L) + ((int128)tmp_q[8] * 2988802878691L) - ((int128)tmp_q[9] * 3780952059490L) + ((int128)tmp_q[10] * 4081351147662L) + ((int128)tmp_q[11] * 6644387537682L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 6644387537682L) + ((int128)tmp_q[1] * 1436766058395L) - ((int128)tmp_q[2] * 161419853902L) + ((int128)tmp_q[3] * 6729544726213L) + ((int128)tmp_q[4] * 5084286511307L) + ((int128)tmp_q[5] * 1277249840956L) + ((int128)tmp_q[6] * 305763315507L) - ((((int128)tmp_q[7] * 251994405312L) - ((int128)tmp_q[8] * 4415422271200L) + ((int128)tmp_q[9] * 2988802878691L) - ((int128)tmp_q[10] * 3780952059490L) + ((int128)tmp_q[11] * 4081351147662L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 4081351147662L) + ((int128)tmp_q[1] * 6644387537682L) + ((int128)tmp_q[2] * 1436766058395L) - ((int128)tmp_q[3] * 161419853902L) + ((int128)tmp_q[4] * 6729544726213L) + ((int128)tmp_q[5] * 5084286511307L) + ((int128)tmp_q[6] * 1277249840956L) + ((int128)tmp_q[7] * 305763315507L) - ((((int128)tmp_q[8] * 251994405312L) - ((int128)tmp_q[9] * 4415422271200L) + ((int128)tmp_q[10] * 2988802878691L) - ((int128)tmp_q[11] * 3780952059490L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 3780952059490L) + ((int128)tmp_q[1] * 4081351147662L) + ((int128)tmp_q[2] * 6644387537682L) + ((int128)tmp_q[3] * 1436766058395L) - ((int128)tmp_q[4] * 161419853902L) + ((int128)tmp_q[5] * 6729544726213L) + ((int128)tmp_q[6] * 5084286511307L) + ((int128)tmp_q[7] * 1277249840956L) + ((int128)tmp_q[8] * 305763315507L) - ((((int128)tmp_q[9] * 251994405312L) - ((int128)tmp_q[10] * 4415422271200L) + ((int128)tmp_q[11] * 2988802878691L)) * 3);
	tmp_zero[9] = ((int128)tmp_q[0] * 2988802878691L) - ((int128)tmp_q[1] * 3780952059490L) + ((int128)tmp_q[2] * 4081351147662L) + ((int128)tmp_q[3] * 6644387537682L) + ((int128)tmp_q[4] * 1436766058395L) - ((int128)tmp_q[5] * 161419853902L) + ((int128)tmp_q[6] * 6729544726213L) + ((int128)tmp_q[7] * 5084286511307L) + ((int128)tmp_q[8] * 1277249840956L) + ((int128)tmp_q[9] * 305763315507L) - ((((int128)tmp_q[10] * 251994405312L) - ((int128)tmp_q[11] * 4415422271200L)) * 3);
	tmp_zero[10] = -((int128)tmp_q[0] * 4415422271200L) + ((int128)tmp_q[1] * 2988802878691L) - ((int128)tmp_q[2] * 3780952059490L) + ((int128)tmp_q[3] * 4081351147662L) + ((int128)tmp_q[4] * 6644387537682L) + ((int128)tmp_q[5] * 1436766058395L) - ((int128)tmp_q[6] * 161419853902L) + ((int128)tmp_q[7] * 6729544726213L) + ((int128)tmp_q[8] * 5084286511307L) + ((int128)tmp_q[9] * 1277249840956L) + ((int128)tmp_q[10] * 305763315507L) - ((int128)tmp_q[11] * 755983215936L);
	tmp_zero[11] = ((int128)tmp_q[0] * 251994405312L) - ((int128)tmp_q[1] * 4415422271200L) + ((int128)tmp_q[2] * 2988802878691L) - ((int128)tmp_q[3] * 3780952059490L) + ((int128)tmp_q[4] * 4081351147662L) + ((int128)tmp_q[5] * 6644387537682L) + ((int128)tmp_q[6] * 1436766058395L) - ((int128)tmp_q[7] * 161419853902L) + ((int128)tmp_q[8] * 6729544726213L) + ((int128)tmp_q[9] * 5084286511307L) + ((int128)tmp_q[10] * 1277249840956L) + ((int128)tmp_q[11] * 305763315507L);

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
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

