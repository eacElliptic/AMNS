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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5179064351418191735UL) + ((((uint64_t)op[1] * 10242265550306402080UL) + ((uint64_t)op[2] * 16882840375278795730UL) + ((uint64_t)op[3] * 16992647286176853823UL) + ((uint64_t)op[4] * 1673958603546851974UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 1673958603546851974UL) + ((uint64_t)op[1] * 5179064351418191735UL) + ((((uint64_t)op[2] * 10242265550306402080UL) + ((uint64_t)op[3] * 16882840375278795730UL) + ((uint64_t)op[4] * 16992647286176853823UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 16992647286176853823UL) + ((uint64_t)op[1] * 1673958603546851974UL) + ((uint64_t)op[2] * 5179064351418191735UL) + ((((uint64_t)op[3] * 10242265550306402080UL) + ((uint64_t)op[4] * 16882840375278795730UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 16882840375278795730UL) + ((uint64_t)op[1] * 16992647286176853823UL) + ((uint64_t)op[2] * 1673958603546851974UL) + ((uint64_t)op[3] * 5179064351418191735UL) + ((uint64_t)op[4] * 8151148107613010176UL);
	tmp_q[4] = ((uint64_t)op[0] * 10242265550306402080UL) + ((uint64_t)op[1] * 16882840375278795730UL) + ((uint64_t)op[2] * 16992647286176853823UL) + ((uint64_t)op[3] * 1673958603546851974UL) + ((uint64_t)op[4] * 5179064351418191735UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4522886751319L) + ((-((int128)tmp_q[1] * 15795518281227L) - ((int128)tmp_q[2] * 3845535361594L) + ((int128)tmp_q[3] * 14012179680323L) - ((int128)tmp_q[4] * 9577264729026L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 9577264729026L) - ((int128)tmp_q[1] * 4522886751319L) + ((-((int128)tmp_q[2] * 15795518281227L) - ((int128)tmp_q[3] * 3845535361594L) + ((int128)tmp_q[4] * 14012179680323L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 14012179680323L) - ((int128)tmp_q[1] * 9577264729026L) - ((int128)tmp_q[2] * 4522886751319L) + ((-((int128)tmp_q[3] * 15795518281227L) - ((int128)tmp_q[4] * 3845535361594L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 3845535361594L) + ((int128)tmp_q[1] * 14012179680323L) - ((int128)tmp_q[2] * 9577264729026L) - ((int128)tmp_q[3] * 4522886751319L) - ((int128)tmp_q[4] * 126364146249816L);
	tmp_zero[4] = -((int128)tmp_q[0] * 15795518281227L) - ((int128)tmp_q[1] * 3845535361594L) + ((int128)tmp_q[2] * 14012179680323L) - ((int128)tmp_q[3] * 9577264729026L) - ((int128)tmp_q[4] * 4522886751319L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

