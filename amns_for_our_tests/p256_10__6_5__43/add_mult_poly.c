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
	tmp_q[0] = ((uint64_t)op[0] * 9361872834377544726UL) + ((((uint64_t)op[1] * 1543989275114648508UL) + ((uint64_t)op[2] * 1770306974485827033UL) + ((uint64_t)op[3] * 17019743932062298293UL) + ((uint64_t)op[4] * 11704915996319895280UL) + ((uint64_t)op[5] * 11056718165057623401UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 11056718165057623401UL) + ((uint64_t)op[1] * 9361872834377544726UL) + ((((uint64_t)op[2] * 1543989275114648508UL) + ((uint64_t)op[3] * 1770306974485827033UL) + ((uint64_t)op[4] * 17019743932062298293UL) + ((uint64_t)op[5] * 11704915996319895280UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 11704915996319895280UL) + ((uint64_t)op[1] * 11056718165057623401UL) + ((uint64_t)op[2] * 9361872834377544726UL) + ((((uint64_t)op[3] * 1543989275114648508UL) + ((uint64_t)op[4] * 1770306974485827033UL) + ((uint64_t)op[5] * 17019743932062298293UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 17019743932062298293UL) + ((uint64_t)op[1] * 11704915996319895280UL) + ((uint64_t)op[2] * 11056718165057623401UL) + ((uint64_t)op[3] * 9361872834377544726UL) + ((((uint64_t)op[4] * 1543989275114648508UL) + ((uint64_t)op[5] * 1770306974485827033UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 1770306974485827033UL) + ((uint64_t)op[1] * 17019743932062298293UL) + ((uint64_t)op[2] * 11704915996319895280UL) + ((uint64_t)op[3] * 11056718165057623401UL) + ((uint64_t)op[4] * 9361872834377544726UL) + ((uint64_t)op[5] * 7719946375573242540UL);
	tmp_q[5] = ((uint64_t)op[0] * 1543989275114648508UL) + ((uint64_t)op[1] * 1770306974485827033UL) + ((uint64_t)op[2] * 17019743932062298293UL) + ((uint64_t)op[3] * 11704915996319895280UL) + ((uint64_t)op[4] * 11056718165057623401UL) + ((uint64_t)op[5] * 9361872834377544726UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4075568064346L) + ((-((int128)tmp_q[1] * 4588332848798L) - ((int128)tmp_q[2] * 1183424282707L) - ((int128)tmp_q[3] * 529610864053L) - ((int128)tmp_q[4] * 3734176554278L) + ((int128)tmp_q[5] * 1853386629197L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1853386629197L) - ((int128)tmp_q[1] * 4075568064346L) + ((-((int128)tmp_q[2] * 4588332848798L) - ((int128)tmp_q[3] * 1183424282707L) - ((int128)tmp_q[4] * 529610864053L) - ((int128)tmp_q[5] * 3734176554278L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 3734176554278L) + ((int128)tmp_q[1] * 1853386629197L) - ((int128)tmp_q[2] * 4075568064346L) + ((-((int128)tmp_q[3] * 4588332848798L) - ((int128)tmp_q[4] * 1183424282707L) - ((int128)tmp_q[5] * 529610864053L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 529610864053L) - ((int128)tmp_q[1] * 3734176554278L) + ((int128)tmp_q[2] * 1853386629197L) - ((int128)tmp_q[3] * 4075568064346L) + ((-((int128)tmp_q[4] * 4588332848798L) - ((int128)tmp_q[5] * 1183424282707L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 1183424282707L) - ((int128)tmp_q[1] * 529610864053L) - ((int128)tmp_q[2] * 3734176554278L) + ((int128)tmp_q[3] * 1853386629197L) - ((int128)tmp_q[4] * 4075568064346L) - ((int128)tmp_q[5] * 22941664243990L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4588332848798L) - ((int128)tmp_q[1] * 1183424282707L) - ((int128)tmp_q[2] * 529610864053L) - ((int128)tmp_q[3] * 3734176554278L) + ((int128)tmp_q[4] * 1853386629197L) - ((int128)tmp_q[5] * 4075568064346L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

