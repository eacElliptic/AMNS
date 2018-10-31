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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15507890400131003375UL) + ((((uint64_t)op[1] * 13696673563691068792UL) + ((uint64_t)op[2] * 4763369230216787317UL) + ((uint64_t)op[3] * 12389084315016052841UL) + ((uint64_t)op[4] * 8909343689634873775UL) + ((uint64_t)op[5] * 13995749126431180033UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 13995749126431180033UL) + ((uint64_t)op[1] * 15507890400131003375UL) + ((((uint64_t)op[2] * 13696673563691068792UL) + ((uint64_t)op[3] * 4763369230216787317UL) + ((uint64_t)op[4] * 12389084315016052841UL) + ((uint64_t)op[5] * 8909343689634873775UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 8909343689634873775UL) + ((uint64_t)op[1] * 13995749126431180033UL) + ((uint64_t)op[2] * 15507890400131003375UL) + ((((uint64_t)op[3] * 13696673563691068792UL) + ((uint64_t)op[4] * 4763369230216787317UL) + ((uint64_t)op[5] * 12389084315016052841UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 12389084315016052841UL) + ((uint64_t)op[1] * 8909343689634873775UL) + ((uint64_t)op[2] * 13995749126431180033UL) + ((uint64_t)op[3] * 15507890400131003375UL) + ((((uint64_t)op[4] * 13696673563691068792UL) + ((uint64_t)op[5] * 4763369230216787317UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 4763369230216787317UL) + ((uint64_t)op[1] * 12389084315016052841UL) + ((uint64_t)op[2] * 8909343689634873775UL) + ((uint64_t)op[3] * 13995749126431180033UL) + ((uint64_t)op[4] * 15507890400131003375UL) + ((uint64_t)op[5] * 3642994577289723464UL);
	tmp_q[5] = ((uint64_t)op[0] * 13696673563691068792UL) + ((uint64_t)op[1] * 4763369230216787317UL) + ((uint64_t)op[2] * 12389084315016052841UL) + ((uint64_t)op[3] * 8909343689634873775UL) + ((uint64_t)op[4] * 13995749126431180033UL) + ((uint64_t)op[5] * 15507890400131003375UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1360151532467L) + ((((int128)tmp_q[1] * 5171292234511L) + ((int128)tmp_q[2] * 895174965559L) - ((int128)tmp_q[3] * 3328815639813L) + ((int128)tmp_q[4] * 6256076767527L) - ((int128)tmp_q[5] * 678902175588L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 678902175588L) + ((int128)tmp_q[1] * 1360151532467L) + ((((int128)tmp_q[2] * 5171292234511L) + ((int128)tmp_q[3] * 895174965559L) - ((int128)tmp_q[4] * 3328815639813L) + ((int128)tmp_q[5] * 6256076767527L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 6256076767527L) - ((int128)tmp_q[1] * 678902175588L) + ((int128)tmp_q[2] * 1360151532467L) + ((((int128)tmp_q[3] * 5171292234511L) + ((int128)tmp_q[4] * 895174965559L) - ((int128)tmp_q[5] * 3328815639813L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 3328815639813L) + ((int128)tmp_q[1] * 6256076767527L) - ((int128)tmp_q[2] * 678902175588L) + ((int128)tmp_q[3] * 1360151532467L) + ((((int128)tmp_q[4] * 5171292234511L) + ((int128)tmp_q[5] * 895174965559L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 895174965559L) - ((int128)tmp_q[1] * 3328815639813L) + ((int128)tmp_q[2] * 6256076767527L) - ((int128)tmp_q[3] * 678902175588L) + ((int128)tmp_q[4] * 1360151532467L) + ((int128)tmp_q[5] * 36199045641577L);
	tmp_zero[5] = ((int128)tmp_q[0] * 5171292234511L) + ((int128)tmp_q[1] * 895174965559L) - ((int128)tmp_q[2] * 3328815639813L) + ((int128)tmp_q[3] * 6256076767527L) - ((int128)tmp_q[4] * 678902175588L) + ((int128)tmp_q[5] * 1360151532467L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

