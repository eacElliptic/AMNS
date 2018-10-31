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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10243021138657805364UL) + ((((uint64_t)op[1] * 13711656053896892399UL) + ((uint64_t)op[2] * 11353963354395397207UL) + ((uint64_t)op[3] * 5823662007008569409UL) + ((uint64_t)op[4] * 486673151884018679UL) + ((uint64_t)op[5] * 7535830812669180051UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 7535830812669180051UL) + ((uint64_t)op[1] * 10243021138657805364UL) + ((((uint64_t)op[2] * 13711656053896892399UL) + ((uint64_t)op[3] * 11353963354395397207UL) + ((uint64_t)op[4] * 5823662007008569409UL) + ((uint64_t)op[5] * 486673151884018679UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 486673151884018679UL) + ((uint64_t)op[1] * 7535830812669180051UL) + ((uint64_t)op[2] * 10243021138657805364UL) + ((((uint64_t)op[3] * 13711656053896892399UL) + ((uint64_t)op[4] * 11353963354395397207UL) + ((uint64_t)op[5] * 5823662007008569409UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 5823662007008569409UL) + ((uint64_t)op[1] * 486673151884018679UL) + ((uint64_t)op[2] * 7535830812669180051UL) + ((uint64_t)op[3] * 10243021138657805364UL) + ((((uint64_t)op[4] * 13711656053896892399UL) + ((uint64_t)op[5] * 11353963354395397207UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 11353963354395397207UL) + ((uint64_t)op[1] * 5823662007008569409UL) + ((uint64_t)op[2] * 486673151884018679UL) + ((uint64_t)op[3] * 7535830812669180051UL) + ((uint64_t)op[4] * 10243021138657805364UL) + ((uint64_t)op[5] * 4241480014271573965UL);
	tmp_q[5] = ((uint64_t)op[0] * 13711656053896892399UL) + ((uint64_t)op[1] * 11353963354395397207UL) + ((uint64_t)op[2] * 5823662007008569409UL) + ((uint64_t)op[3] * 486673151884018679UL) + ((uint64_t)op[4] * 7535830812669180051UL) + ((uint64_t)op[5] * 10243021138657805364UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4813532L) + ((-((int128)tmp_q[1] * 549465669L) + ((int128)tmp_q[2] * 1434019469L) + ((int128)tmp_q[3] * 3367901973L) - ((int128)tmp_q[4] * 340884107L) + ((int128)tmp_q[5] * 2171632815L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 2171632815L) - ((int128)tmp_q[1] * 4813532L) + ((-((int128)tmp_q[2] * 549465669L) + ((int128)tmp_q[3] * 1434019469L) + ((int128)tmp_q[4] * 3367901973L) - ((int128)tmp_q[5] * 340884107L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 340884107L) + ((int128)tmp_q[1] * 2171632815L) - ((int128)tmp_q[2] * 4813532L) + ((-((int128)tmp_q[3] * 549465669L) + ((int128)tmp_q[4] * 1434019469L) + ((int128)tmp_q[5] * 3367901973L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 3367901973L) - ((int128)tmp_q[1] * 340884107L) + ((int128)tmp_q[2] * 2171632815L) - ((int128)tmp_q[3] * 4813532L) + ((-((int128)tmp_q[4] * 549465669L) + ((int128)tmp_q[5] * 1434019469L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1434019469L) + ((int128)tmp_q[1] * 3367901973L) - ((int128)tmp_q[2] * 340884107L) + ((int128)tmp_q[3] * 2171632815L) - ((int128)tmp_q[4] * 4813532L) - ((int128)tmp_q[5] * 1648397007L);
	tmp_zero[5] = -((int128)tmp_q[0] * 549465669L) + ((int128)tmp_q[1] * 1434019469L) + ((int128)tmp_q[2] * 3367901973L) - ((int128)tmp_q[3] * 340884107L) + ((int128)tmp_q[4] * 2171632815L) - ((int128)tmp_q[5] * 4813532L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

