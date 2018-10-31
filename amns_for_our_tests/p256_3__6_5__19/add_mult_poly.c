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
	tmp_q[0] = ((uint64_t)op[0] * 7725962757684489642UL) + ((((uint64_t)op[1] * 8338158098366143336UL) + ((uint64_t)op[2] * 2211100251135257706UL) + ((uint64_t)op[3] * 12290438885079254020UL) + ((uint64_t)op[4] * 3165279922311718743UL) + ((uint64_t)op[5] * 16799679511817836280UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 16799679511817836280UL) + ((uint64_t)op[1] * 7725962757684489642UL) + ((((uint64_t)op[2] * 8338158098366143336UL) + ((uint64_t)op[3] * 2211100251135257706UL) + ((uint64_t)op[4] * 12290438885079254020UL) + ((uint64_t)op[5] * 3165279922311718743UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 3165279922311718743UL) + ((uint64_t)op[1] * 16799679511817836280UL) + ((uint64_t)op[2] * 7725962757684489642UL) + ((((uint64_t)op[3] * 8338158098366143336UL) + ((uint64_t)op[4] * 2211100251135257706UL) + ((uint64_t)op[5] * 12290438885079254020UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 12290438885079254020UL) + ((uint64_t)op[1] * 3165279922311718743UL) + ((uint64_t)op[2] * 16799679511817836280UL) + ((uint64_t)op[3] * 7725962757684489642UL) + ((((uint64_t)op[4] * 8338158098366143336UL) + ((uint64_t)op[5] * 2211100251135257706UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 2211100251135257706UL) + ((uint64_t)op[1] * 12290438885079254020UL) + ((uint64_t)op[2] * 3165279922311718743UL) + ((uint64_t)op[3] * 16799679511817836280UL) + ((uint64_t)op[4] * 7725962757684489642UL) + ((uint64_t)op[5] * 4797302344411613448UL);
	tmp_q[5] = ((uint64_t)op[0] * 8338158098366143336UL) + ((uint64_t)op[1] * 2211100251135257706UL) + ((uint64_t)op[2] * 12290438885079254020UL) + ((uint64_t)op[3] * 3165279922311718743UL) + ((uint64_t)op[4] * 16799679511817836280UL) + ((uint64_t)op[5] * 7725962757684489642UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4891053490698L) + ((((int128)tmp_q[1] * 3217366122708L) - ((int128)tmp_q[2] * 1071404126339L) + ((int128)tmp_q[3] * 5656753521112L) + ((int128)tmp_q[4] * 2011588412718L) + ((int128)tmp_q[5] * 2800137368936L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 2800137368936L) - ((int128)tmp_q[1] * 4891053490698L) + ((((int128)tmp_q[2] * 3217366122708L) - ((int128)tmp_q[3] * 1071404126339L) + ((int128)tmp_q[4] * 5656753521112L) + ((int128)tmp_q[5] * 2011588412718L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 2011588412718L) + ((int128)tmp_q[1] * 2800137368936L) - ((int128)tmp_q[2] * 4891053490698L) + ((((int128)tmp_q[3] * 3217366122708L) - ((int128)tmp_q[4] * 1071404126339L) + ((int128)tmp_q[5] * 5656753521112L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 5656753521112L) + ((int128)tmp_q[1] * 2011588412718L) + ((int128)tmp_q[2] * 2800137368936L) - ((int128)tmp_q[3] * 4891053490698L) + ((((int128)tmp_q[4] * 3217366122708L) - ((int128)tmp_q[5] * 1071404126339L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 1071404126339L) + ((int128)tmp_q[1] * 5656753521112L) + ((int128)tmp_q[2] * 2011588412718L) + ((int128)tmp_q[3] * 2800137368936L) - ((int128)tmp_q[4] * 4891053490698L) + ((int128)tmp_q[5] * 16086830613540L);
	tmp_zero[5] = ((int128)tmp_q[0] * 3217366122708L) - ((int128)tmp_q[1] * 1071404126339L) + ((int128)tmp_q[2] * 5656753521112L) + ((int128)tmp_q[3] * 2011588412718L) + ((int128)tmp_q[4] * 2800137368936L) - ((int128)tmp_q[5] * 4891053490698L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

