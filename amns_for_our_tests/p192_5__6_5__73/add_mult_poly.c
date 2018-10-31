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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16768392806365583444UL) + ((((uint64_t)op[1] * 3925909087157594210UL) + ((uint64_t)op[2] * 1407527245456972115UL) + ((uint64_t)op[3] * 13995746538069368282UL) + ((uint64_t)op[4] * 13838609441069565483UL) + ((uint64_t)op[5] * 8649093784893311429UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 8649093784893311429UL) + ((uint64_t)op[1] * 16768392806365583444UL) + ((((uint64_t)op[2] * 3925909087157594210UL) + ((uint64_t)op[3] * 1407527245456972115UL) + ((uint64_t)op[4] * 13995746538069368282UL) + ((uint64_t)op[5] * 13838609441069565483UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13838609441069565483UL) + ((uint64_t)op[1] * 8649093784893311429UL) + ((uint64_t)op[2] * 16768392806365583444UL) + ((((uint64_t)op[3] * 3925909087157594210UL) + ((uint64_t)op[4] * 1407527245456972115UL) + ((uint64_t)op[5] * 13995746538069368282UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 13995746538069368282UL) + ((uint64_t)op[1] * 13838609441069565483UL) + ((uint64_t)op[2] * 8649093784893311429UL) + ((uint64_t)op[3] * 16768392806365583444UL) + ((((uint64_t)op[4] * 3925909087157594210UL) + ((uint64_t)op[5] * 1407527245456972115UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 1407527245456972115UL) + ((uint64_t)op[1] * 13995746538069368282UL) + ((uint64_t)op[2] * 13838609441069565483UL) + ((uint64_t)op[3] * 8649093784893311429UL) + ((uint64_t)op[4] * 16768392806365583444UL) + ((uint64_t)op[5] * 1182801362078419434UL);
	tmp_q[5] = ((uint64_t)op[0] * 3925909087157594210UL) + ((uint64_t)op[1] * 1407527245456972115UL) + ((uint64_t)op[2] * 13995746538069368282UL) + ((uint64_t)op[3] * 13838609441069565483UL) + ((uint64_t)op[4] * 8649093784893311429UL) + ((uint64_t)op[5] * 16768392806365583444UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 859252843L) + ((((int128)tmp_q[1] * 549304040L) + ((int128)tmp_q[2] * 1600923481L) + ((int128)tmp_q[3] * 1354240125L) - ((int128)tmp_q[4] * 76085494L) - ((int128)tmp_q[5] * 1347668434L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 1347668434L) + ((int128)tmp_q[1] * 859252843L) + ((((int128)tmp_q[2] * 549304040L) + ((int128)tmp_q[3] * 1600923481L) + ((int128)tmp_q[4] * 1354240125L) - ((int128)tmp_q[5] * 76085494L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 76085494L) - ((int128)tmp_q[1] * 1347668434L) + ((int128)tmp_q[2] * 859252843L) + ((((int128)tmp_q[3] * 549304040L) + ((int128)tmp_q[4] * 1600923481L) + ((int128)tmp_q[5] * 1354240125L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1354240125L) - ((int128)tmp_q[1] * 76085494L) - ((int128)tmp_q[2] * 1347668434L) + ((int128)tmp_q[3] * 859252843L) + ((((int128)tmp_q[4] * 549304040L) + ((int128)tmp_q[5] * 1600923481L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1600923481L) + ((int128)tmp_q[1] * 1354240125L) - ((int128)tmp_q[2] * 76085494L) - ((int128)tmp_q[3] * 1347668434L) + ((int128)tmp_q[4] * 859252843L) + ((int128)tmp_q[5] * 2746520200L);
	tmp_zero[5] = ((int128)tmp_q[0] * 549304040L) + ((int128)tmp_q[1] * 1600923481L) + ((int128)tmp_q[2] * 1354240125L) - ((int128)tmp_q[3] * 76085494L) - ((int128)tmp_q[4] * 1347668434L) + ((int128)tmp_q[5] * 859252843L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

