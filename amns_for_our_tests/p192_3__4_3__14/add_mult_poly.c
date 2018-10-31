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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[3] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17452808223785974636UL) + ((((uint64_t)op[1] * 3320793618174945616UL) + ((uint64_t)op[2] * 10872905950601895525UL) + ((uint64_t)op[3] * 15277842632240448424UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 15277842632240448424UL) + ((uint64_t)op[1] * 17452808223785974636UL) + ((((uint64_t)op[2] * 3320793618174945616UL) + ((uint64_t)op[3] * 10872905950601895525UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 10872905950601895525UL) + ((uint64_t)op[1] * 15277842632240448424UL) + ((uint64_t)op[2] * 17452808223785974636UL) + ((uint64_t)op[3] * 9962380854524836848UL);
	tmp_q[3] = ((uint64_t)op[0] * 3320793618174945616UL) + ((uint64_t)op[1] * 10872905950601895525UL) + ((uint64_t)op[2] * 15277842632240448424UL) + ((uint64_t)op[3] * 17452808223785974636UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 91282236674620L) + ((((int128)tmp_q[1] * 93112523539504L) + ((int128)tmp_q[2] * 191785893977825L) - ((int128)tmp_q[3] * 17280853361928L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 17280853361928L) - ((int128)tmp_q[1] * 91282236674620L) + ((((int128)tmp_q[2] * 93112523539504L) + ((int128)tmp_q[3] * 191785893977825L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 191785893977825L) - ((int128)tmp_q[1] * 17280853361928L) - ((int128)tmp_q[2] * 91282236674620L) + ((int128)tmp_q[3] * 279337570618512L);
	tmp_zero[3] = ((int128)tmp_q[0] * 93112523539504L) + ((int128)tmp_q[1] * 191785893977825L) - ((int128)tmp_q[2] * 17280853361928L) - ((int128)tmp_q[3] * 91282236674620L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

