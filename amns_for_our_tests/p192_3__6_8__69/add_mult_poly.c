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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17043481124025746701UL) + ((((uint64_t)op[1] * 2544038061232150316UL) + ((uint64_t)op[2] * 6858279546203315670UL) + ((uint64_t)op[3] * 5296542555972285826UL) + ((uint64_t)op[4] * 6151498464916166641UL) + ((uint64_t)op[5] * 15417880498150267724UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 15417880498150267724UL) + ((uint64_t)op[1] * 17043481124025746701UL) + ((((uint64_t)op[2] * 2544038061232150316UL) + ((uint64_t)op[3] * 6858279546203315670UL) + ((uint64_t)op[4] * 5296542555972285826UL) + ((uint64_t)op[5] * 6151498464916166641UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 6151498464916166641UL) + ((uint64_t)op[1] * 15417880498150267724UL) + ((uint64_t)op[2] * 17043481124025746701UL) + ((((uint64_t)op[3] * 2544038061232150316UL) + ((uint64_t)op[4] * 6858279546203315670UL) + ((uint64_t)op[5] * 5296542555972285826UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 5296542555972285826UL) + ((uint64_t)op[1] * 6151498464916166641UL) + ((uint64_t)op[2] * 15417880498150267724UL) + ((uint64_t)op[3] * 17043481124025746701UL) + ((((uint64_t)op[4] * 2544038061232150316UL) + ((uint64_t)op[5] * 6858279546203315670UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 6858279546203315670UL) + ((uint64_t)op[1] * 5296542555972285826UL) + ((uint64_t)op[2] * 6151498464916166641UL) + ((uint64_t)op[3] * 15417880498150267724UL) + ((uint64_t)op[4] * 17043481124025746701UL) + ((uint64_t)op[5] * 1905560416147650912UL);
	tmp_q[5] = ((uint64_t)op[0] * 2544038061232150316UL) + ((uint64_t)op[1] * 6858279546203315670UL) + ((uint64_t)op[2] * 5296542555972285826UL) + ((uint64_t)op[3] * 6151498464916166641UL) + ((uint64_t)op[4] * 15417880498150267724UL) + ((uint64_t)op[5] * 17043481124025746701UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 35474025527165L) + ((-((int128)tmp_q[1] * 59902972877156L) - ((int128)tmp_q[2] * 16122398761295L) - ((int128)tmp_q[3] * 18844336668294L) + ((int128)tmp_q[4] * 79009250434977L) + ((int128)tmp_q[5] * 24397911553436L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 24397911553436L) - ((int128)tmp_q[1] * 35474025527165L) + ((-((int128)tmp_q[2] * 59902972877156L) - ((int128)tmp_q[3] * 16122398761295L) - ((int128)tmp_q[4] * 18844336668294L) + ((int128)tmp_q[5] * 79009250434977L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 79009250434977L) + ((int128)tmp_q[1] * 24397911553436L) - ((int128)tmp_q[2] * 35474025527165L) + ((-((int128)tmp_q[3] * 59902972877156L) - ((int128)tmp_q[4] * 16122398761295L) - ((int128)tmp_q[5] * 18844336668294L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 18844336668294L) + ((int128)tmp_q[1] * 79009250434977L) + ((int128)tmp_q[2] * 24397911553436L) - ((int128)tmp_q[3] * 35474025527165L) + ((-((int128)tmp_q[4] * 59902972877156L) - ((int128)tmp_q[5] * 16122398761295L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 16122398761295L) - ((int128)tmp_q[1] * 18844336668294L) + ((int128)tmp_q[2] * 79009250434977L) + ((int128)tmp_q[3] * 24397911553436L) - ((int128)tmp_q[4] * 35474025527165L) - ((int128)tmp_q[5] * 479223783017248L);
	tmp_zero[5] = -((int128)tmp_q[0] * 59902972877156L) - ((int128)tmp_q[1] * 16122398761295L) - ((int128)tmp_q[2] * 18844336668294L) + ((int128)tmp_q[3] * 79009250434977L) + ((int128)tmp_q[4] * 24397911553436L) - ((int128)tmp_q[5] * 35474025527165L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

