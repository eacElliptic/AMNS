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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 7);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 14);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 14);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 7);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6776747070860136416UL) + ((((uint64_t)op[1] * 6051276546230758437UL) + ((uint64_t)op[2] * 10575070291882238638UL) + ((uint64_t)op[3] * 5855663889769609941UL) + ((uint64_t)op[4] * 18113040722864828597UL) + ((uint64_t)op[5] * 7981546177527876445UL) + ((uint64_t)op[6] * 10663111378470515808UL) + ((uint64_t)op[7] * 9054448456414183531UL) + ((uint64_t)op[8] * 10826458217989004336UL) + ((uint64_t)op[9] * 10662818522160244511UL) + ((uint64_t)op[10] * 11425847754594596236UL) + ((uint64_t)op[11] * 3205519207618457197UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 3205519207618457197UL) + ((uint64_t)op[1] * 6776747070860136416UL) + ((((uint64_t)op[2] * 6051276546230758437UL) + ((uint64_t)op[3] * 10575070291882238638UL) + ((uint64_t)op[4] * 5855663889769609941UL) + ((uint64_t)op[5] * 18113040722864828597UL) + ((uint64_t)op[6] * 7981546177527876445UL) + ((uint64_t)op[7] * 10663111378470515808UL) + ((uint64_t)op[8] * 9054448456414183531UL) + ((uint64_t)op[9] * 10826458217989004336UL) + ((uint64_t)op[10] * 10662818522160244511UL) + ((uint64_t)op[11] * 11425847754594596236UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 11425847754594596236UL) + ((uint64_t)op[1] * 3205519207618457197UL) + ((uint64_t)op[2] * 6776747070860136416UL) + ((((uint64_t)op[3] * 6051276546230758437UL) + ((uint64_t)op[4] * 10575070291882238638UL) + ((uint64_t)op[5] * 5855663889769609941UL) + ((uint64_t)op[6] * 18113040722864828597UL) + ((uint64_t)op[7] * 7981546177527876445UL) + ((uint64_t)op[8] * 10663111378470515808UL) + ((uint64_t)op[9] * 9054448456414183531UL) + ((uint64_t)op[10] * 10826458217989004336UL) + ((uint64_t)op[11] * 10662818522160244511UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 10662818522160244511UL) + ((uint64_t)op[1] * 11425847754594596236UL) + ((uint64_t)op[2] * 3205519207618457197UL) + ((uint64_t)op[3] * 6776747070860136416UL) + ((((uint64_t)op[4] * 6051276546230758437UL) + ((uint64_t)op[5] * 10575070291882238638UL) + ((uint64_t)op[6] * 5855663889769609941UL) + ((uint64_t)op[7] * 18113040722864828597UL) + ((uint64_t)op[8] * 7981546177527876445UL) + ((uint64_t)op[9] * 10663111378470515808UL) + ((uint64_t)op[10] * 9054448456414183531UL) + ((uint64_t)op[11] * 10826458217989004336UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 10826458217989004336UL) + ((uint64_t)op[1] * 10662818522160244511UL) + ((uint64_t)op[2] * 11425847754594596236UL) + ((uint64_t)op[3] * 3205519207618457197UL) + ((uint64_t)op[4] * 6776747070860136416UL) + ((((uint64_t)op[5] * 6051276546230758437UL) + ((uint64_t)op[6] * 10575070291882238638UL) + ((uint64_t)op[7] * 5855663889769609941UL) + ((uint64_t)op[8] * 18113040722864828597UL) + ((uint64_t)op[9] * 7981546177527876445UL) + ((uint64_t)op[10] * 10663111378470515808UL) + ((uint64_t)op[11] * 9054448456414183531UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 9054448456414183531UL) + ((uint64_t)op[1] * 10826458217989004336UL) + ((uint64_t)op[2] * 10662818522160244511UL) + ((uint64_t)op[3] * 11425847754594596236UL) + ((uint64_t)op[4] * 3205519207618457197UL) + ((uint64_t)op[5] * 6776747070860136416UL) + ((((uint64_t)op[6] * 6051276546230758437UL) + ((uint64_t)op[7] * 10575070291882238638UL) + ((uint64_t)op[8] * 5855663889769609941UL) + ((uint64_t)op[9] * 18113040722864828597UL) + ((uint64_t)op[10] * 7981546177527876445UL) + ((uint64_t)op[11] * 10663111378470515808UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 10663111378470515808UL) + ((uint64_t)op[1] * 9054448456414183531UL) + ((uint64_t)op[2] * 10826458217989004336UL) + ((uint64_t)op[3] * 10662818522160244511UL) + ((uint64_t)op[4] * 11425847754594596236UL) + ((uint64_t)op[5] * 3205519207618457197UL) + ((uint64_t)op[6] * 6776747070860136416UL) + ((((uint64_t)op[7] * 6051276546230758437UL) + ((uint64_t)op[8] * 10575070291882238638UL) + ((uint64_t)op[9] * 5855663889769609941UL) + ((uint64_t)op[10] * 18113040722864828597UL) + ((uint64_t)op[11] * 7981546177527876445UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 7981546177527876445UL) + ((uint64_t)op[1] * 10663111378470515808UL) + ((uint64_t)op[2] * 9054448456414183531UL) + ((uint64_t)op[3] * 10826458217989004336UL) + ((uint64_t)op[4] * 10662818522160244511UL) + ((uint64_t)op[5] * 11425847754594596236UL) + ((uint64_t)op[6] * 3205519207618457197UL) + ((uint64_t)op[7] * 6776747070860136416UL) + ((((uint64_t)op[8] * 6051276546230758437UL) + ((uint64_t)op[9] * 10575070291882238638UL) + ((uint64_t)op[10] * 5855663889769609941UL) + ((uint64_t)op[11] * 18113040722864828597UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 18113040722864828597UL) + ((uint64_t)op[1] * 7981546177527876445UL) + ((uint64_t)op[2] * 10663111378470515808UL) + ((uint64_t)op[3] * 9054448456414183531UL) + ((uint64_t)op[4] * 10826458217989004336UL) + ((uint64_t)op[5] * 10662818522160244511UL) + ((uint64_t)op[6] * 11425847754594596236UL) + ((uint64_t)op[7] * 3205519207618457197UL) + ((uint64_t)op[8] * 6776747070860136416UL) + ((((uint64_t)op[9] * 6051276546230758437UL) + ((uint64_t)op[10] * 10575070291882238638UL) + ((uint64_t)op[11] * 5855663889769609941UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 5855663889769609941UL) + ((uint64_t)op[1] * 18113040722864828597UL) + ((uint64_t)op[2] * 7981546177527876445UL) + ((uint64_t)op[3] * 10663111378470515808UL) + ((uint64_t)op[4] * 9054448456414183531UL) + ((uint64_t)op[5] * 10826458217989004336UL) + ((uint64_t)op[6] * 10662818522160244511UL) + ((uint64_t)op[7] * 11425847754594596236UL) + ((uint64_t)op[8] * 3205519207618457197UL) + ((uint64_t)op[9] * 6776747070860136416UL) + ((((uint64_t)op[10] * 6051276546230758437UL) + ((uint64_t)op[11] * 10575070291882238638UL)) * 7);
	tmp_q[10] = ((uint64_t)op[0] * 10575070291882238638UL) + ((uint64_t)op[1] * 5855663889769609941UL) + ((uint64_t)op[2] * 18113040722864828597UL) + ((uint64_t)op[3] * 7981546177527876445UL) + ((uint64_t)op[4] * 10663111378470515808UL) + ((uint64_t)op[5] * 9054448456414183531UL) + ((uint64_t)op[6] * 10826458217989004336UL) + ((uint64_t)op[7] * 10662818522160244511UL) + ((uint64_t)op[8] * 11425847754594596236UL) + ((uint64_t)op[9] * 3205519207618457197UL) + ((uint64_t)op[10] * 6776747070860136416UL) + ((uint64_t)op[11] * 5465447676196205827UL);
	tmp_q[11] = ((uint64_t)op[0] * 6051276546230758437UL) + ((uint64_t)op[1] * 10575070291882238638UL) + ((uint64_t)op[2] * 5855663889769609941UL) + ((uint64_t)op[3] * 18113040722864828597UL) + ((uint64_t)op[4] * 7981546177527876445UL) + ((uint64_t)op[5] * 10663111378470515808UL) + ((uint64_t)op[6] * 9054448456414183531UL) + ((uint64_t)op[7] * 10826458217989004336UL) + ((uint64_t)op[8] * 10662818522160244511UL) + ((uint64_t)op[9] * 11425847754594596236UL) + ((uint64_t)op[10] * 3205519207618457197UL) + ((uint64_t)op[11] * 6776747070860136416UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4754275579706L) + ((((int128)tmp_q[1] * 4175316162817L) - ((int128)tmp_q[2] * 6233655042184L) - ((int128)tmp_q[3] * 6317285983905L) + ((int128)tmp_q[4] * 388705012470L) - ((int128)tmp_q[5] * 984669767031L) - ((int128)tmp_q[6] * 3602056533710L) - ((int128)tmp_q[7] * 5607212224967L) + ((int128)tmp_q[8] * 2052178233219L) - ((int128)tmp_q[9] * 4120494057131L) - ((int128)tmp_q[10] * 543002214576L) + ((int128)tmp_q[11] * 4814121546703L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 4814121546703L) + ((int128)tmp_q[1] * 4754275579706L) + ((((int128)tmp_q[2] * 4175316162817L) - ((int128)tmp_q[3] * 6233655042184L) - ((int128)tmp_q[4] * 6317285983905L) + ((int128)tmp_q[5] * 388705012470L) - ((int128)tmp_q[6] * 984669767031L) - ((int128)tmp_q[7] * 3602056533710L) - ((int128)tmp_q[8] * 5607212224967L) + ((int128)tmp_q[9] * 2052178233219L) - ((int128)tmp_q[10] * 4120494057131L) - ((int128)tmp_q[11] * 543002214576L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 543002214576L) + ((int128)tmp_q[1] * 4814121546703L) + ((int128)tmp_q[2] * 4754275579706L) + ((((int128)tmp_q[3] * 4175316162817L) - ((int128)tmp_q[4] * 6233655042184L) - ((int128)tmp_q[5] * 6317285983905L) + ((int128)tmp_q[6] * 388705012470L) - ((int128)tmp_q[7] * 984669767031L) - ((int128)tmp_q[8] * 3602056533710L) - ((int128)tmp_q[9] * 5607212224967L) + ((int128)tmp_q[10] * 2052178233219L) - ((int128)tmp_q[11] * 4120494057131L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 4120494057131L) - ((int128)tmp_q[1] * 543002214576L) + ((int128)tmp_q[2] * 4814121546703L) + ((int128)tmp_q[3] * 4754275579706L) + ((((int128)tmp_q[4] * 4175316162817L) - ((int128)tmp_q[5] * 6233655042184L) - ((int128)tmp_q[6] * 6317285983905L) + ((int128)tmp_q[7] * 388705012470L) - ((int128)tmp_q[8] * 984669767031L) - ((int128)tmp_q[9] * 3602056533710L) - ((int128)tmp_q[10] * 5607212224967L) + ((int128)tmp_q[11] * 2052178233219L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 2052178233219L) - ((int128)tmp_q[1] * 4120494057131L) - ((int128)tmp_q[2] * 543002214576L) + ((int128)tmp_q[3] * 4814121546703L) + ((int128)tmp_q[4] * 4754275579706L) + ((((int128)tmp_q[5] * 4175316162817L) - ((int128)tmp_q[6] * 6233655042184L) - ((int128)tmp_q[7] * 6317285983905L) + ((int128)tmp_q[8] * 388705012470L) - ((int128)tmp_q[9] * 984669767031L) - ((int128)tmp_q[10] * 3602056533710L) - ((int128)tmp_q[11] * 5607212224967L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 5607212224967L) + ((int128)tmp_q[1] * 2052178233219L) - ((int128)tmp_q[2] * 4120494057131L) - ((int128)tmp_q[3] * 543002214576L) + ((int128)tmp_q[4] * 4814121546703L) + ((int128)tmp_q[5] * 4754275579706L) + ((((int128)tmp_q[6] * 4175316162817L) - ((int128)tmp_q[7] * 6233655042184L) - ((int128)tmp_q[8] * 6317285983905L) + ((int128)tmp_q[9] * 388705012470L) - ((int128)tmp_q[10] * 984669767031L) - ((int128)tmp_q[11] * 3602056533710L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 3602056533710L) - ((int128)tmp_q[1] * 5607212224967L) + ((int128)tmp_q[2] * 2052178233219L) - ((int128)tmp_q[3] * 4120494057131L) - ((int128)tmp_q[4] * 543002214576L) + ((int128)tmp_q[5] * 4814121546703L) + ((int128)tmp_q[6] * 4754275579706L) + ((((int128)tmp_q[7] * 4175316162817L) - ((int128)tmp_q[8] * 6233655042184L) - ((int128)tmp_q[9] * 6317285983905L) + ((int128)tmp_q[10] * 388705012470L) - ((int128)tmp_q[11] * 984669767031L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 984669767031L) - ((int128)tmp_q[1] * 3602056533710L) - ((int128)tmp_q[2] * 5607212224967L) + ((int128)tmp_q[3] * 2052178233219L) - ((int128)tmp_q[4] * 4120494057131L) - ((int128)tmp_q[5] * 543002214576L) + ((int128)tmp_q[6] * 4814121546703L) + ((int128)tmp_q[7] * 4754275579706L) + ((((int128)tmp_q[8] * 4175316162817L) - ((int128)tmp_q[9] * 6233655042184L) - ((int128)tmp_q[10] * 6317285983905L) + ((int128)tmp_q[11] * 388705012470L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 388705012470L) - ((int128)tmp_q[1] * 984669767031L) - ((int128)tmp_q[2] * 3602056533710L) - ((int128)tmp_q[3] * 5607212224967L) + ((int128)tmp_q[4] * 2052178233219L) - ((int128)tmp_q[5] * 4120494057131L) - ((int128)tmp_q[6] * 543002214576L) + ((int128)tmp_q[7] * 4814121546703L) + ((int128)tmp_q[8] * 4754275579706L) + ((((int128)tmp_q[9] * 4175316162817L) - ((int128)tmp_q[10] * 6233655042184L) - ((int128)tmp_q[11] * 6317285983905L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 6317285983905L) + ((int128)tmp_q[1] * 388705012470L) - ((int128)tmp_q[2] * 984669767031L) - ((int128)tmp_q[3] * 3602056533710L) - ((int128)tmp_q[4] * 5607212224967L) + ((int128)tmp_q[5] * 2052178233219L) - ((int128)tmp_q[6] * 4120494057131L) - ((int128)tmp_q[7] * 543002214576L) + ((int128)tmp_q[8] * 4814121546703L) + ((int128)tmp_q[9] * 4754275579706L) + ((((int128)tmp_q[10] * 4175316162817L) - ((int128)tmp_q[11] * 6233655042184L)) * 7);
	tmp_zero[10] = -((int128)tmp_q[0] * 6233655042184L) - ((int128)tmp_q[1] * 6317285983905L) + ((int128)tmp_q[2] * 388705012470L) - ((int128)tmp_q[3] * 984669767031L) - ((int128)tmp_q[4] * 3602056533710L) - ((int128)tmp_q[5] * 5607212224967L) + ((int128)tmp_q[6] * 2052178233219L) - ((int128)tmp_q[7] * 4120494057131L) - ((int128)tmp_q[8] * 543002214576L) + ((int128)tmp_q[9] * 4814121546703L) + ((int128)tmp_q[10] * 4754275579706L) + ((int128)tmp_q[11] * 29227213139719L);
	tmp_zero[11] = ((int128)tmp_q[0] * 4175316162817L) - ((int128)tmp_q[1] * 6233655042184L) - ((int128)tmp_q[2] * 6317285983905L) + ((int128)tmp_q[3] * 388705012470L) - ((int128)tmp_q[4] * 984669767031L) - ((int128)tmp_q[5] * 3602056533710L) - ((int128)tmp_q[6] * 5607212224967L) + ((int128)tmp_q[7] * 2052178233219L) - ((int128)tmp_q[8] * 4120494057131L) - ((int128)tmp_q[9] * 543002214576L) + ((int128)tmp_q[10] * 4814121546703L) + ((int128)tmp_q[11] * 4754275579706L);

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

