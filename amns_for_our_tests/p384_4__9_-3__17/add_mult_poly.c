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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12113326314446553913UL) + ((((uint64_t)op[1] * 9445215545442246895UL) + ((uint64_t)op[2] * 15302614963493398261UL) + ((uint64_t)op[3] * 12010479723982394374UL) + ((uint64_t)op[4] * 16769658320858369504UL) + ((uint64_t)op[5] * 12940187175818087714UL) + ((uint64_t)op[6] * 9932275039712982233UL) + ((uint64_t)op[7] * 16218010324248581234UL) + ((uint64_t)op[8] * 17067600261465572109UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 17067600261465572109UL) + ((uint64_t)op[1] * 12113326314446553913UL) + ((((uint64_t)op[2] * 9445215545442246895UL) + ((uint64_t)op[3] * 15302614963493398261UL) + ((uint64_t)op[4] * 12010479723982394374UL) + ((uint64_t)op[5] * 16769658320858369504UL) + ((uint64_t)op[6] * 12940187175818087714UL) + ((uint64_t)op[7] * 9932275039712982233UL) + ((uint64_t)op[8] * 16218010324248581234UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 16218010324248581234UL) + ((uint64_t)op[1] * 17067600261465572109UL) + ((uint64_t)op[2] * 12113326314446553913UL) + ((((uint64_t)op[3] * 9445215545442246895UL) + ((uint64_t)op[4] * 15302614963493398261UL) + ((uint64_t)op[5] * 12010479723982394374UL) + ((uint64_t)op[6] * 16769658320858369504UL) + ((uint64_t)op[7] * 12940187175818087714UL) + ((uint64_t)op[8] * 9932275039712982233UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 9932275039712982233UL) + ((uint64_t)op[1] * 16218010324248581234UL) + ((uint64_t)op[2] * 17067600261465572109UL) + ((uint64_t)op[3] * 12113326314446553913UL) + ((((uint64_t)op[4] * 9445215545442246895UL) + ((uint64_t)op[5] * 15302614963493398261UL) + ((uint64_t)op[6] * 12010479723982394374UL) + ((uint64_t)op[7] * 16769658320858369504UL) + ((uint64_t)op[8] * 12940187175818087714UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 12940187175818087714UL) + ((uint64_t)op[1] * 9932275039712982233UL) + ((uint64_t)op[2] * 16218010324248581234UL) + ((uint64_t)op[3] * 17067600261465572109UL) + ((uint64_t)op[4] * 12113326314446553913UL) + ((((uint64_t)op[5] * 9445215545442246895UL) + ((uint64_t)op[6] * 15302614963493398261UL) + ((uint64_t)op[7] * 12010479723982394374UL) + ((uint64_t)op[8] * 16769658320858369504UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 16769658320858369504UL) + ((uint64_t)op[1] * 12940187175818087714UL) + ((uint64_t)op[2] * 9932275039712982233UL) + ((uint64_t)op[3] * 16218010324248581234UL) + ((uint64_t)op[4] * 17067600261465572109UL) + ((uint64_t)op[5] * 12113326314446553913UL) + ((((uint64_t)op[6] * 9445215545442246895UL) + ((uint64_t)op[7] * 15302614963493398261UL) + ((uint64_t)op[8] * 12010479723982394374UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 12010479723982394374UL) + ((uint64_t)op[1] * 16769658320858369504UL) + ((uint64_t)op[2] * 12940187175818087714UL) + ((uint64_t)op[3] * 9932275039712982233UL) + ((uint64_t)op[4] * 16218010324248581234UL) + ((uint64_t)op[5] * 17067600261465572109UL) + ((uint64_t)op[6] * 12113326314446553913UL) + ((((uint64_t)op[7] * 9445215545442246895UL) + ((uint64_t)op[8] * 15302614963493398261UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 15302614963493398261UL) + ((uint64_t)op[1] * 12010479723982394374UL) + ((uint64_t)op[2] * 16769658320858369504UL) + ((uint64_t)op[3] * 12940187175818087714UL) + ((uint64_t)op[4] * 9932275039712982233UL) + ((uint64_t)op[5] * 16218010324248581234UL) + ((uint64_t)op[6] * 17067600261465572109UL) + ((uint64_t)op[7] * 12113326314446553913UL) + ((uint64_t)op[8] * 8557841511092362547UL);
	tmp_q[8] = ((uint64_t)op[0] * 9445215545442246895UL) + ((uint64_t)op[1] * 15302614963493398261UL) + ((uint64_t)op[2] * 12010479723982394374UL) + ((uint64_t)op[3] * 16769658320858369504UL) + ((uint64_t)op[4] * 12940187175818087714UL) + ((uint64_t)op[5] * 9932275039712982233UL) + ((uint64_t)op[6] * 16218010324248581234UL) + ((uint64_t)op[7] * 17067600261465572109UL) + ((uint64_t)op[8] * 12113326314446553913UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 950337045020L) - ((((int128)tmp_q[1] * 2740369006951L) + ((int128)tmp_q[2] * 3364752943782L) + ((int128)tmp_q[3] * 2080054203194L) + ((int128)tmp_q[4] * 875457292305L) - ((int128)tmp_q[5] * 3423011889999L) + ((int128)tmp_q[6] * 4805200429158L) + ((int128)tmp_q[7] * 4182985665400L) + ((int128)tmp_q[8] * 746534907632L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 746534907632L) + ((int128)tmp_q[1] * 950337045020L) - ((((int128)tmp_q[2] * 2740369006951L) + ((int128)tmp_q[3] * 3364752943782L) + ((int128)tmp_q[4] * 2080054203194L) + ((int128)tmp_q[5] * 875457292305L) - ((int128)tmp_q[6] * 3423011889999L) + ((int128)tmp_q[7] * 4805200429158L) + ((int128)tmp_q[8] * 4182985665400L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 4182985665400L) + ((int128)tmp_q[1] * 746534907632L) + ((int128)tmp_q[2] * 950337045020L) - ((((int128)tmp_q[3] * 2740369006951L) + ((int128)tmp_q[4] * 3364752943782L) + ((int128)tmp_q[5] * 2080054203194L) + ((int128)tmp_q[6] * 875457292305L) - ((int128)tmp_q[7] * 3423011889999L) + ((int128)tmp_q[8] * 4805200429158L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 4805200429158L) + ((int128)tmp_q[1] * 4182985665400L) + ((int128)tmp_q[2] * 746534907632L) + ((int128)tmp_q[3] * 950337045020L) - ((((int128)tmp_q[4] * 2740369006951L) + ((int128)tmp_q[5] * 3364752943782L) + ((int128)tmp_q[6] * 2080054203194L) + ((int128)tmp_q[7] * 875457292305L) - ((int128)tmp_q[8] * 3423011889999L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 3423011889999L) + ((int128)tmp_q[1] * 4805200429158L) + ((int128)tmp_q[2] * 4182985665400L) + ((int128)tmp_q[3] * 746534907632L) + ((int128)tmp_q[4] * 950337045020L) - ((((int128)tmp_q[5] * 2740369006951L) + ((int128)tmp_q[6] * 3364752943782L) + ((int128)tmp_q[7] * 2080054203194L) + ((int128)tmp_q[8] * 875457292305L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 875457292305L) - ((int128)tmp_q[1] * 3423011889999L) + ((int128)tmp_q[2] * 4805200429158L) + ((int128)tmp_q[3] * 4182985665400L) + ((int128)tmp_q[4] * 746534907632L) + ((int128)tmp_q[5] * 950337045020L) - ((((int128)tmp_q[6] * 2740369006951L) + ((int128)tmp_q[7] * 3364752943782L) + ((int128)tmp_q[8] * 2080054203194L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 2080054203194L) + ((int128)tmp_q[1] * 875457292305L) - ((int128)tmp_q[2] * 3423011889999L) + ((int128)tmp_q[3] * 4805200429158L) + ((int128)tmp_q[4] * 4182985665400L) + ((int128)tmp_q[5] * 746534907632L) + ((int128)tmp_q[6] * 950337045020L) - ((((int128)tmp_q[7] * 2740369006951L) + ((int128)tmp_q[8] * 3364752943782L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 3364752943782L) + ((int128)tmp_q[1] * 2080054203194L) + ((int128)tmp_q[2] * 875457292305L) - ((int128)tmp_q[3] * 3423011889999L) + ((int128)tmp_q[4] * 4805200429158L) + ((int128)tmp_q[5] * 4182985665400L) + ((int128)tmp_q[6] * 746534907632L) + ((int128)tmp_q[7] * 950337045020L) - ((int128)tmp_q[8] * 8221107020853L);
	tmp_zero[8] = ((int128)tmp_q[0] * 2740369006951L) + ((int128)tmp_q[1] * 3364752943782L) + ((int128)tmp_q[2] * 2080054203194L) + ((int128)tmp_q[3] * 875457292305L) - ((int128)tmp_q[4] * 3423011889999L) + ((int128)tmp_q[5] * 4805200429158L) + ((int128)tmp_q[6] * 4182985665400L) + ((int128)tmp_q[7] * 746534907632L) + ((int128)tmp_q[8] * 950337045020L);

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
}

