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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9972679492486699989UL) + ((((uint64_t)op[1] * 4389512573700126710UL) + ((uint64_t)op[2] * 13502309057433309589UL) + ((uint64_t)op[3] * 4914298761371297028UL) + ((uint64_t)op[4] * 15811913075197461524UL) + ((uint64_t)op[5] * 10518359118283236949UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 10518359118283236949UL) + ((uint64_t)op[1] * 9972679492486699989UL) + ((((uint64_t)op[2] * 4389512573700126710UL) + ((uint64_t)op[3] * 13502309057433309589UL) + ((uint64_t)op[4] * 4914298761371297028UL) + ((uint64_t)op[5] * 15811913075197461524UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 15811913075197461524UL) + ((uint64_t)op[1] * 10518359118283236949UL) + ((uint64_t)op[2] * 9972679492486699989UL) + ((((uint64_t)op[3] * 4389512573700126710UL) + ((uint64_t)op[4] * 13502309057433309589UL) + ((uint64_t)op[5] * 4914298761371297028UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 4914298761371297028UL) + ((uint64_t)op[1] * 15811913075197461524UL) + ((uint64_t)op[2] * 10518359118283236949UL) + ((uint64_t)op[3] * 9972679492486699989UL) + ((((uint64_t)op[4] * 4389512573700126710UL) + ((uint64_t)op[5] * 13502309057433309589UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 13502309057433309589UL) + ((uint64_t)op[1] * 4914298761371297028UL) + ((uint64_t)op[2] * 15811913075197461524UL) + ((uint64_t)op[3] * 10518359118283236949UL) + ((uint64_t)op[4] * 9972679492486699989UL) + ((uint64_t)op[5] * 7890331368491208644UL);
	tmp_q[5] = ((uint64_t)op[0] * 4389512573700126710UL) + ((uint64_t)op[1] * 13502309057433309589UL) + ((uint64_t)op[2] * 4914298761371297028UL) + ((uint64_t)op[3] * 15811913075197461524UL) + ((uint64_t)op[4] * 10518359118283236949UL) + ((uint64_t)op[5] * 9972679492486699989UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 109048193615L) + ((((int128)tmp_q[1] * 81153224039L) - ((int128)tmp_q[2] * 94218192600L) + ((int128)tmp_q[3] * 9906161261L) - ((int128)tmp_q[4] * 51438102687L) + ((int128)tmp_q[5] * 73046274823L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 73046274823L) + ((int128)tmp_q[1] * 109048193615L) + ((((int128)tmp_q[2] * 81153224039L) - ((int128)tmp_q[3] * 94218192600L) + ((int128)tmp_q[4] * 9906161261L) - ((int128)tmp_q[5] * 51438102687L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 51438102687L) + ((int128)tmp_q[1] * 73046274823L) + ((int128)tmp_q[2] * 109048193615L) + ((((int128)tmp_q[3] * 81153224039L) - ((int128)tmp_q[4] * 94218192600L) + ((int128)tmp_q[5] * 9906161261L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 9906161261L) - ((int128)tmp_q[1] * 51438102687L) + ((int128)tmp_q[2] * 73046274823L) + ((int128)tmp_q[3] * 109048193615L) + ((((int128)tmp_q[4] * 81153224039L) - ((int128)tmp_q[5] * 94218192600L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 94218192600L) + ((int128)tmp_q[1] * 9906161261L) - ((int128)tmp_q[2] * 51438102687L) + ((int128)tmp_q[3] * 73046274823L) + ((int128)tmp_q[4] * 109048193615L) + ((int128)tmp_q[5] * 486919344234L);
	tmp_zero[5] = ((int128)tmp_q[0] * 81153224039L) - ((int128)tmp_q[1] * 94218192600L) + ((int128)tmp_q[2] * 9906161261L) - ((int128)tmp_q[3] * 51438102687L) + ((int128)tmp_q[4] * 73046274823L) + ((int128)tmp_q[5] * 109048193615L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

