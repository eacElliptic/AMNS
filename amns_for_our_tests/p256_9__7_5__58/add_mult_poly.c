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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 909595925563797539UL) + ((((uint64_t)op[1] * 891908557457746762UL) + ((uint64_t)op[2] * 3342459565277150559UL) + ((uint64_t)op[3] * 12310709338569831509UL) + ((uint64_t)op[4] * 13981288858818843935UL) + ((uint64_t)op[5] * 4109232940764438050UL) + ((uint64_t)op[6] * 10020048482822618971UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 10020048482822618971UL) + ((uint64_t)op[1] * 909595925563797539UL) + ((((uint64_t)op[2] * 891908557457746762UL) + ((uint64_t)op[3] * 3342459565277150559UL) + ((uint64_t)op[4] * 12310709338569831509UL) + ((uint64_t)op[5] * 13981288858818843935UL) + ((uint64_t)op[6] * 4109232940764438050UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 4109232940764438050UL) + ((uint64_t)op[1] * 10020048482822618971UL) + ((uint64_t)op[2] * 909595925563797539UL) + ((((uint64_t)op[3] * 891908557457746762UL) + ((uint64_t)op[4] * 3342459565277150559UL) + ((uint64_t)op[5] * 12310709338569831509UL) + ((uint64_t)op[6] * 13981288858818843935UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 13981288858818843935UL) + ((uint64_t)op[1] * 4109232940764438050UL) + ((uint64_t)op[2] * 10020048482822618971UL) + ((uint64_t)op[3] * 909595925563797539UL) + ((((uint64_t)op[4] * 891908557457746762UL) + ((uint64_t)op[5] * 3342459565277150559UL) + ((uint64_t)op[6] * 12310709338569831509UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 12310709338569831509UL) + ((uint64_t)op[1] * 13981288858818843935UL) + ((uint64_t)op[2] * 4109232940764438050UL) + ((uint64_t)op[3] * 10020048482822618971UL) + ((uint64_t)op[4] * 909595925563797539UL) + ((((uint64_t)op[5] * 891908557457746762UL) + ((uint64_t)op[6] * 3342459565277150559UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 3342459565277150559UL) + ((uint64_t)op[1] * 12310709338569831509UL) + ((uint64_t)op[2] * 13981288858818843935UL) + ((uint64_t)op[3] * 4109232940764438050UL) + ((uint64_t)op[4] * 10020048482822618971UL) + ((uint64_t)op[5] * 909595925563797539UL) + ((uint64_t)op[6] * 4459542787288733810UL);
	tmp_q[6] = ((uint64_t)op[0] * 891908557457746762UL) + ((uint64_t)op[1] * 3342459565277150559UL) + ((uint64_t)op[2] * 12310709338569831509UL) + ((uint64_t)op[3] * 13981288858818843935UL) + ((uint64_t)op[4] * 4109232940764438050UL) + ((uint64_t)op[5] * 10020048482822618971UL) + ((uint64_t)op[6] * 909595925563797539UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 43683561728L) + ((-((int128)tmp_q[1] * 32931941244L) - ((int128)tmp_q[2] * 28841211062L) + ((int128)tmp_q[3] * 38485075691L) + ((int128)tmp_q[4] * 51328367901L) + ((int128)tmp_q[5] * 16208001615L) - ((int128)tmp_q[6] * 47203684978L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 47203684978L) - ((int128)tmp_q[1] * 43683561728L) + ((-((int128)tmp_q[2] * 32931941244L) - ((int128)tmp_q[3] * 28841211062L) + ((int128)tmp_q[4] * 38485075691L) + ((int128)tmp_q[5] * 51328367901L) + ((int128)tmp_q[6] * 16208001615L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 16208001615L) - ((int128)tmp_q[1] * 47203684978L) - ((int128)tmp_q[2] * 43683561728L) + ((-((int128)tmp_q[3] * 32931941244L) - ((int128)tmp_q[4] * 28841211062L) + ((int128)tmp_q[5] * 38485075691L) + ((int128)tmp_q[6] * 51328367901L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 51328367901L) + ((int128)tmp_q[1] * 16208001615L) - ((int128)tmp_q[2] * 47203684978L) - ((int128)tmp_q[3] * 43683561728L) + ((-((int128)tmp_q[4] * 32931941244L) - ((int128)tmp_q[5] * 28841211062L) + ((int128)tmp_q[6] * 38485075691L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 38485075691L) + ((int128)tmp_q[1] * 51328367901L) + ((int128)tmp_q[2] * 16208001615L) - ((int128)tmp_q[3] * 47203684978L) - ((int128)tmp_q[4] * 43683561728L) + ((-((int128)tmp_q[5] * 32931941244L) - ((int128)tmp_q[6] * 28841211062L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 28841211062L) + ((int128)tmp_q[1] * 38485075691L) + ((int128)tmp_q[2] * 51328367901L) + ((int128)tmp_q[3] * 16208001615L) - ((int128)tmp_q[4] * 47203684978L) - ((int128)tmp_q[5] * 43683561728L) - ((int128)tmp_q[6] * 164659706220L);
	tmp_zero[6] = -((int128)tmp_q[0] * 32931941244L) - ((int128)tmp_q[1] * 28841211062L) + ((int128)tmp_q[2] * 38485075691L) + ((int128)tmp_q[3] * 51328367901L) + ((int128)tmp_q[4] * 16208001615L) - ((int128)tmp_q[5] * 47203684978L) - ((int128)tmp_q[6] * 43683561728L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

