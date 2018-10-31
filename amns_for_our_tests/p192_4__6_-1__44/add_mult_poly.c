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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17799548802849042539UL) + ((((uint64_t)op[1] * 10415349516312193518UL) + ((uint64_t)op[2] * 7564940420058751337UL) + ((uint64_t)op[3] * 8514600525396670734UL) + ((uint64_t)op[4] * 6917745149198242260UL) + ((uint64_t)op[5] * 16545995082794028833UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 16545995082794028833UL) + ((uint64_t)op[1] * 17799548802849042539UL) + ((((uint64_t)op[2] * 10415349516312193518UL) + ((uint64_t)op[3] * 7564940420058751337UL) + ((uint64_t)op[4] * 8514600525396670734UL) + ((uint64_t)op[5] * 6917745149198242260UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 6917745149198242260UL) + ((uint64_t)op[1] * 16545995082794028833UL) + ((uint64_t)op[2] * 17799548802849042539UL) + ((((uint64_t)op[3] * 10415349516312193518UL) + ((uint64_t)op[4] * 7564940420058751337UL) + ((uint64_t)op[5] * 8514600525396670734UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 8514600525396670734UL) + ((uint64_t)op[1] * 6917745149198242260UL) + ((uint64_t)op[2] * 16545995082794028833UL) + ((uint64_t)op[3] * 17799548802849042539UL) + ((((uint64_t)op[4] * 10415349516312193518UL) + ((uint64_t)op[5] * 7564940420058751337UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 7564940420058751337UL) + ((uint64_t)op[1] * 8514600525396670734UL) + ((uint64_t)op[2] * 6917745149198242260UL) + ((uint64_t)op[3] * 16545995082794028833UL) + ((uint64_t)op[4] * 17799548802849042539UL) + ((uint64_t)op[5] * 8031394557397358098UL);
	tmp_q[5] = ((uint64_t)op[0] * 10415349516312193518UL) + ((uint64_t)op[1] * 7564940420058751337UL) + ((uint64_t)op[2] * 8514600525396670734UL) + ((uint64_t)op[3] * 6917745149198242260UL) + ((uint64_t)op[4] * 16545995082794028833UL) + ((uint64_t)op[5] * 17799548802849042539UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 161360114276075L) - (-((int128)tmp_q[1] * 79338020744496L) - ((int128)tmp_q[2] * 113162679049875L) + ((int128)tmp_q[3] * 51234141398304L) + ((int128)tmp_q[4] * 48197435226200L) + ((int128)tmp_q[5] * 130572162142801L));
	tmp_zero[1] = ((int128)tmp_q[0] * 130572162142801L) + ((int128)tmp_q[1] * 161360114276075L) - (-((int128)tmp_q[2] * 79338020744496L) - ((int128)tmp_q[3] * 113162679049875L) + ((int128)tmp_q[4] * 51234141398304L) + ((int128)tmp_q[5] * 48197435226200L));
	tmp_zero[2] = ((int128)tmp_q[0] * 48197435226200L) + ((int128)tmp_q[1] * 130572162142801L) + ((int128)tmp_q[2] * 161360114276075L) - (-((int128)tmp_q[3] * 79338020744496L) - ((int128)tmp_q[4] * 113162679049875L) + ((int128)tmp_q[5] * 51234141398304L));
	tmp_zero[3] = ((int128)tmp_q[0] * 51234141398304L) + ((int128)tmp_q[1] * 48197435226200L) + ((int128)tmp_q[2] * 130572162142801L) + ((int128)tmp_q[3] * 161360114276075L) - (-((int128)tmp_q[4] * 79338020744496L) - ((int128)tmp_q[5] * 113162679049875L));
	tmp_zero[4] = -((int128)tmp_q[0] * 113162679049875L) + ((int128)tmp_q[1] * 51234141398304L) + ((int128)tmp_q[2] * 48197435226200L) + ((int128)tmp_q[3] * 130572162142801L) + ((int128)tmp_q[4] * 161360114276075L) + ((int128)tmp_q[5] * 79338020744496L);
	tmp_zero[5] = -((int128)tmp_q[0] * 79338020744496L) - ((int128)tmp_q[1] * 113162679049875L) + ((int128)tmp_q[2] * 51234141398304L) + ((int128)tmp_q[3] * 48197435226200L) + ((int128)tmp_q[4] * 130572162142801L) + ((int128)tmp_q[5] * 161360114276075L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

