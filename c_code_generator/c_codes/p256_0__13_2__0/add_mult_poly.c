#include "add_mult_poly.h"


void add_poly(int *rop, int *pa, int *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int *rop, int *pa, int *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

void neg_poly(int *rop, int *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int *rop, int *op, int scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ Computes pa(X)*pb(X) mod(X^n - c)
void mult_mod_poly(int *rop, int *pa, int *pb){

	llong tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (llong)pa[0] * pb[0] + (((llong)pa[1] * pb[12] + (llong)pa[2] * pb[11] + (llong)pa[3] * pb[10] + (llong)pa[4] * pb[9] + (llong)pa[5] * pb[8] + (llong)pa[6] * pb[7] + (llong)pa[7] * pb[6] + (llong)pa[8] * pb[5] + (llong)pa[9] * pb[4] + (llong)pa[10] * pb[3] + (llong)pa[11] * pb[2] + (llong)pa[12] * pb[1]) << 1);
	tmp_prod_result[1] = (llong)pa[0] * pb[1] + (llong)pa[1] * pb[0] + (((llong)pa[2] * pb[12] + (llong)pa[3] * pb[11] + (llong)pa[4] * pb[10] + (llong)pa[5] * pb[9] + (llong)pa[6] * pb[8] + (llong)pa[7] * pb[7] + (llong)pa[8] * pb[6] + (llong)pa[9] * pb[5] + (llong)pa[10] * pb[4] + (llong)pa[11] * pb[3] + (llong)pa[12] * pb[2]) << 1);
	tmp_prod_result[2] = (llong)pa[0] * pb[2] + (llong)pa[1] * pb[1] + (llong)pa[2] * pb[0] + (((llong)pa[3] * pb[12] + (llong)pa[4] * pb[11] + (llong)pa[5] * pb[10] + (llong)pa[6] * pb[9] + (llong)pa[7] * pb[8] + (llong)pa[8] * pb[7] + (llong)pa[9] * pb[6] + (llong)pa[10] * pb[5] + (llong)pa[11] * pb[4] + (llong)pa[12] * pb[3]) << 1);
	tmp_prod_result[3] = (llong)pa[0] * pb[3] + (llong)pa[1] * pb[2] + (llong)pa[2] * pb[1] + (llong)pa[3] * pb[0] + (((llong)pa[4] * pb[12] + (llong)pa[5] * pb[11] + (llong)pa[6] * pb[10] + (llong)pa[7] * pb[9] + (llong)pa[8] * pb[8] + (llong)pa[9] * pb[7] + (llong)pa[10] * pb[6] + (llong)pa[11] * pb[5] + (llong)pa[12] * pb[4]) << 1);
	tmp_prod_result[4] = (llong)pa[0] * pb[4] + (llong)pa[1] * pb[3] + (llong)pa[2] * pb[2] + (llong)pa[3] * pb[1] + (llong)pa[4] * pb[0] + (((llong)pa[5] * pb[12] + (llong)pa[6] * pb[11] + (llong)pa[7] * pb[10] + (llong)pa[8] * pb[9] + (llong)pa[9] * pb[8] + (llong)pa[10] * pb[7] + (llong)pa[11] * pb[6] + (llong)pa[12] * pb[5]) << 1);
	tmp_prod_result[5] = (llong)pa[0] * pb[5] + (llong)pa[1] * pb[4] + (llong)pa[2] * pb[3] + (llong)pa[3] * pb[2] + (llong)pa[4] * pb[1] + (llong)pa[5] * pb[0] + (((llong)pa[6] * pb[12] + (llong)pa[7] * pb[11] + (llong)pa[8] * pb[10] + (llong)pa[9] * pb[9] + (llong)pa[10] * pb[8] + (llong)pa[11] * pb[7] + (llong)pa[12] * pb[6]) << 1);
	tmp_prod_result[6] = (llong)pa[0] * pb[6] + (llong)pa[1] * pb[5] + (llong)pa[2] * pb[4] + (llong)pa[3] * pb[3] + (llong)pa[4] * pb[2] + (llong)pa[5] * pb[1] + (llong)pa[6] * pb[0] + (((llong)pa[7] * pb[12] + (llong)pa[8] * pb[11] + (llong)pa[9] * pb[10] + (llong)pa[10] * pb[9] + (llong)pa[11] * pb[8] + (llong)pa[12] * pb[7]) << 1);
	tmp_prod_result[7] = (llong)pa[0] * pb[7] + (llong)pa[1] * pb[6] + (llong)pa[2] * pb[5] + (llong)pa[3] * pb[4] + (llong)pa[4] * pb[3] + (llong)pa[5] * pb[2] + (llong)pa[6] * pb[1] + (llong)pa[7] * pb[0] + (((llong)pa[8] * pb[12] + (llong)pa[9] * pb[11] + (llong)pa[10] * pb[10] + (llong)pa[11] * pb[9] + (llong)pa[12] * pb[8]) << 1);
	tmp_prod_result[8] = (llong)pa[0] * pb[8] + (llong)pa[1] * pb[7] + (llong)pa[2] * pb[6] + (llong)pa[3] * pb[5] + (llong)pa[4] * pb[4] + (llong)pa[5] * pb[3] + (llong)pa[6] * pb[2] + (llong)pa[7] * pb[1] + (llong)pa[8] * pb[0] + (((llong)pa[9] * pb[12] + (llong)pa[10] * pb[11] + (llong)pa[11] * pb[10] + (llong)pa[12] * pb[9]) << 1);
	tmp_prod_result[9] = (llong)pa[0] * pb[9] + (llong)pa[1] * pb[8] + (llong)pa[2] * pb[7] + (llong)pa[3] * pb[6] + (llong)pa[4] * pb[5] + (llong)pa[5] * pb[4] + (llong)pa[6] * pb[3] + (llong)pa[7] * pb[2] + (llong)pa[8] * pb[1] + (llong)pa[9] * pb[0] + (((llong)pa[10] * pb[12] + (llong)pa[11] * pb[11] + (llong)pa[12] * pb[10]) << 1);
	tmp_prod_result[10] = (llong)pa[0] * pb[10] + (llong)pa[1] * pb[9] + (llong)pa[2] * pb[8] + (llong)pa[3] * pb[7] + (llong)pa[4] * pb[6] + (llong)pa[5] * pb[5] + (llong)pa[6] * pb[4] + (llong)pa[7] * pb[3] + (llong)pa[8] * pb[2] + (llong)pa[9] * pb[1] + (llong)pa[10] * pb[0] + (((llong)pa[11] * pb[12] + (llong)pa[12] * pb[11]) << 1);
	tmp_prod_result[11] = (llong)pa[0] * pb[11] + (llong)pa[1] * pb[10] + (llong)pa[2] * pb[9] + (llong)pa[3] * pb[8] + (llong)pa[4] * pb[7] + (llong)pa[5] * pb[6] + (llong)pa[6] * pb[5] + (llong)pa[7] * pb[4] + (llong)pa[8] * pb[3] + (llong)pa[9] * pb[2] + (llong)pa[10] * pb[1] + (llong)pa[11] * pb[0] + (((llong)pa[12] * pb[12]) << 1);
	tmp_prod_result[12] = (llong)pa[0] * pb[12] + (llong)pa[1] * pb[11] + (llong)pa[2] * pb[10] + (llong)pa[3] * pb[9] + (llong)pa[4] * pb[8] + (llong)pa[5] * pb[7] + (llong)pa[6] * pb[6] + (llong)pa[7] * pb[5] + (llong)pa[8] * pb[4] + (llong)pa[9] * pb[3] + (llong)pa[10] * pb[2] + (llong)pa[11] * pb[1] + (llong)pa[12] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int *rop, int *pa){

	llong tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (llong)pa[0] * pa[0] + (((llong)pa[7] * pa[6] + (llong)pa[8] * pa[5] + (llong)pa[9] * pa[4] + (llong)pa[10] * pa[3] + (llong)pa[11] * pa[2] + (llong)pa[12] * pa[1]) << 2);
	tmp_prod_result[1] = (((llong)pa[1] * pa[0]) << 1) + (((((llong)pa[8] * pa[6] + (llong)pa[9] * pa[5] + (llong)pa[10] * pa[4] + (llong)pa[11] * pa[3] + (llong)pa[12] * pa[2]) << 1) + (llong)pa[7] * pa[7]) << 1);
	tmp_prod_result[2] = (((llong)pa[2] * pa[0]) << 1) + (llong)pa[1] * pa[1] + (((llong)pa[8] * pa[7] + (llong)pa[9] * pa[6] + (llong)pa[10] * pa[5] + (llong)pa[11] * pa[4] + (llong)pa[12] * pa[3]) << 2);
	tmp_prod_result[3] = (((llong)pa[2] * pa[1] + (llong)pa[3] * pa[0]) << 1) + (((((llong)pa[9] * pa[7] + (llong)pa[10] * pa[6] + (llong)pa[11] * pa[5] + (llong)pa[12] * pa[4]) << 1) + (llong)pa[8] * pa[8]) << 1);
	tmp_prod_result[4] = (((llong)pa[3] * pa[1] + (llong)pa[4] * pa[0]) << 1) + (llong)pa[2] * pa[2] + (((llong)pa[9] * pa[8] + (llong)pa[10] * pa[7] + (llong)pa[11] * pa[6] + (llong)pa[12] * pa[5]) << 2);
	tmp_prod_result[5] = (((llong)pa[3] * pa[2] + (llong)pa[4] * pa[1] + (llong)pa[5] * pa[0]) << 1) + (((((llong)pa[10] * pa[8] + (llong)pa[11] * pa[7] + (llong)pa[12] * pa[6]) << 1) + (llong)pa[9] * pa[9]) << 1);
	tmp_prod_result[6] = (((llong)pa[4] * pa[2] + (llong)pa[5] * pa[1] + (llong)pa[6] * pa[0]) << 1) + (llong)pa[3] * pa[3] + (((llong)pa[10] * pa[9] + (llong)pa[11] * pa[8] + (llong)pa[12] * pa[7]) << 2);
	tmp_prod_result[7] = (((llong)pa[4] * pa[3] + (llong)pa[5] * pa[2] + (llong)pa[6] * pa[1] + (llong)pa[7] * pa[0]) << 1) + (((((llong)pa[11] * pa[9] + (llong)pa[12] * pa[8]) << 1) + (llong)pa[10] * pa[10]) << 1);
	tmp_prod_result[8] = (((llong)pa[5] * pa[3] + (llong)pa[6] * pa[2] + (llong)pa[7] * pa[1] + (llong)pa[8] * pa[0]) << 1) + (llong)pa[4] * pa[4] + (((llong)pa[11] * pa[10] + (llong)pa[12] * pa[9]) << 2);
	tmp_prod_result[9] = (((llong)pa[5] * pa[4] + (llong)pa[6] * pa[3] + (llong)pa[7] * pa[2] + (llong)pa[8] * pa[1] + (llong)pa[9] * pa[0]) << 1) + (((((llong)pa[12] * pa[10]) << 1) + (llong)pa[11] * pa[11]) << 1);
	tmp_prod_result[10] = (((llong)pa[6] * pa[4] + (llong)pa[7] * pa[3] + (llong)pa[8] * pa[2] + (llong)pa[9] * pa[1] + (llong)pa[10] * pa[0]) << 1) + (llong)pa[5] * pa[5] + (((llong)pa[12] * pa[11]) << 2);
	tmp_prod_result[11] = (((llong)pa[6] * pa[5] + (llong)pa[7] * pa[4] + (llong)pa[8] * pa[3] + (llong)pa[9] * pa[2] + (llong)pa[10] * pa[1] + (llong)pa[11] * pa[0]) << 1) + (((llong)pa[12] * pa[12]) << 1);
	tmp_prod_result[12] = (((llong)pa[7] * pa[5] + (llong)pa[8] * pa[4] + (llong)pa[9] * pa[3] + (llong)pa[10] * pa[2] + (llong)pa[11] * pa[1] + (llong)pa[12] * pa[0]) << 1) + (llong)pa[6] * pa[6];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int *rop, llong *op){

	uint tmp_q[NB_COEFF];
	llong tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint)op[0] * 1049726001UL) + ((((uint)op[1] * 281380021UL) + ((uint)op[2] * 3298942106UL) + ((uint)op[3] * 1138226165UL) + ((uint)op[4] * 1318620905UL) + ((uint)op[5] * 130940849UL) + ((uint)op[6] * 2777209316UL) + ((uint)op[7] * 3461746829UL) + ((uint)op[8] * 1825549316UL) + ((uint)op[9] * 32802490UL) + ((uint)op[10] * 2090964479UL) + ((uint)op[11] * 789633811UL) + ((uint)op[12] * 3481490052UL)) * 2);
	tmp_q[1] = ((uint)op[0] * 3481490052UL) + ((uint)op[1] * 1049726001UL) + ((((uint)op[2] * 281380021UL) + ((uint)op[3] * 3298942106UL) + ((uint)op[4] * 1138226165UL) + ((uint)op[5] * 1318620905UL) + ((uint)op[6] * 130940849UL) + ((uint)op[7] * 2777209316UL) + ((uint)op[8] * 3461746829UL) + ((uint)op[9] * 1825549316UL) + ((uint)op[10] * 32802490UL) + ((uint)op[11] * 2090964479UL) + ((uint)op[12] * 789633811UL)) * 2);
	tmp_q[2] = ((uint)op[0] * 789633811UL) + ((uint)op[1] * 3481490052UL) + ((uint)op[2] * 1049726001UL) + ((((uint)op[3] * 281380021UL) + ((uint)op[4] * 3298942106UL) + ((uint)op[5] * 1138226165UL) + ((uint)op[6] * 1318620905UL) + ((uint)op[7] * 130940849UL) + ((uint)op[8] * 2777209316UL) + ((uint)op[9] * 3461746829UL) + ((uint)op[10] * 1825549316UL) + ((uint)op[11] * 32802490UL) + ((uint)op[12] * 2090964479UL)) * 2);
	tmp_q[3] = ((uint)op[0] * 2090964479UL) + ((uint)op[1] * 789633811UL) + ((uint)op[2] * 3481490052UL) + ((uint)op[3] * 1049726001UL) + ((((uint)op[4] * 281380021UL) + ((uint)op[5] * 3298942106UL) + ((uint)op[6] * 1138226165UL) + ((uint)op[7] * 1318620905UL) + ((uint)op[8] * 130940849UL) + ((uint)op[9] * 2777209316UL) + ((uint)op[10] * 3461746829UL) + ((uint)op[11] * 1825549316UL) + ((uint)op[12] * 32802490UL)) * 2);
	tmp_q[4] = ((uint)op[0] * 32802490UL) + ((uint)op[1] * 2090964479UL) + ((uint)op[2] * 789633811UL) + ((uint)op[3] * 3481490052UL) + ((uint)op[4] * 1049726001UL) + ((((uint)op[5] * 281380021UL) + ((uint)op[6] * 3298942106UL) + ((uint)op[7] * 1138226165UL) + ((uint)op[8] * 1318620905UL) + ((uint)op[9] * 130940849UL) + ((uint)op[10] * 2777209316UL) + ((uint)op[11] * 3461746829UL) + ((uint)op[12] * 1825549316UL)) * 2);
	tmp_q[5] = ((uint)op[0] * 1825549316UL) + ((uint)op[1] * 32802490UL) + ((uint)op[2] * 2090964479UL) + ((uint)op[3] * 789633811UL) + ((uint)op[4] * 3481490052UL) + ((uint)op[5] * 1049726001UL) + ((((uint)op[6] * 281380021UL) + ((uint)op[7] * 3298942106UL) + ((uint)op[8] * 1138226165UL) + ((uint)op[9] * 1318620905UL) + ((uint)op[10] * 130940849UL) + ((uint)op[11] * 2777209316UL) + ((uint)op[12] * 3461746829UL)) * 2);
	tmp_q[6] = ((uint)op[0] * 3461746829UL) + ((uint)op[1] * 1825549316UL) + ((uint)op[2] * 32802490UL) + ((uint)op[3] * 2090964479UL) + ((uint)op[4] * 789633811UL) + ((uint)op[5] * 3481490052UL) + ((uint)op[6] * 1049726001UL) + ((((uint)op[7] * 281380021UL) + ((uint)op[8] * 3298942106UL) + ((uint)op[9] * 1138226165UL) + ((uint)op[10] * 1318620905UL) + ((uint)op[11] * 130940849UL) + ((uint)op[12] * 2777209316UL)) * 2);
	tmp_q[7] = ((uint)op[0] * 2777209316UL) + ((uint)op[1] * 3461746829UL) + ((uint)op[2] * 1825549316UL) + ((uint)op[3] * 32802490UL) + ((uint)op[4] * 2090964479UL) + ((uint)op[5] * 789633811UL) + ((uint)op[6] * 3481490052UL) + ((uint)op[7] * 1049726001UL) + ((((uint)op[8] * 281380021UL) + ((uint)op[9] * 3298942106UL) + ((uint)op[10] * 1138226165UL) + ((uint)op[11] * 1318620905UL) + ((uint)op[12] * 130940849UL)) * 2);
	tmp_q[8] = ((uint)op[0] * 130940849UL) + ((uint)op[1] * 2777209316UL) + ((uint)op[2] * 3461746829UL) + ((uint)op[3] * 1825549316UL) + ((uint)op[4] * 32802490UL) + ((uint)op[5] * 2090964479UL) + ((uint)op[6] * 789633811UL) + ((uint)op[7] * 3481490052UL) + ((uint)op[8] * 1049726001UL) + ((((uint)op[9] * 281380021UL) + ((uint)op[10] * 3298942106UL) + ((uint)op[11] * 1138226165UL) + ((uint)op[12] * 1318620905UL)) * 2);
	tmp_q[9] = ((uint)op[0] * 1318620905UL) + ((uint)op[1] * 130940849UL) + ((uint)op[2] * 2777209316UL) + ((uint)op[3] * 3461746829UL) + ((uint)op[4] * 1825549316UL) + ((uint)op[5] * 32802490UL) + ((uint)op[6] * 2090964479UL) + ((uint)op[7] * 789633811UL) + ((uint)op[8] * 3481490052UL) + ((uint)op[9] * 1049726001UL) + ((((uint)op[10] * 281380021UL) + ((uint)op[11] * 3298942106UL) + ((uint)op[12] * 1138226165UL)) * 2);
	tmp_q[10] = ((uint)op[0] * 1138226165UL) + ((uint)op[1] * 1318620905UL) + ((uint)op[2] * 130940849UL) + ((uint)op[3] * 2777209316UL) + ((uint)op[4] * 3461746829UL) + ((uint)op[5] * 1825549316UL) + ((uint)op[6] * 32802490UL) + ((uint)op[7] * 2090964479UL) + ((uint)op[8] * 789633811UL) + ((uint)op[9] * 3481490052UL) + ((uint)op[10] * 1049726001UL) + ((((uint)op[11] * 281380021UL) + ((uint)op[12] * 3298942106UL)) * 2);
	tmp_q[11] = ((uint)op[0] * 3298942106UL) + ((uint)op[1] * 1138226165UL) + ((uint)op[2] * 1318620905UL) + ((uint)op[3] * 130940849UL) + ((uint)op[4] * 2777209316UL) + ((uint)op[5] * 3461746829UL) + ((uint)op[6] * 1825549316UL) + ((uint)op[7] * 32802490UL) + ((uint)op[8] * 2090964479UL) + ((uint)op[9] * 789633811UL) + ((uint)op[10] * 3481490052UL) + ((uint)op[11] * 1049726001UL) + ((uint)op[12] * 562760042UL);
	tmp_q[12] = ((uint)op[0] * 281380021UL) + ((uint)op[1] * 3298942106UL) + ((uint)op[2] * 1138226165UL) + ((uint)op[3] * 1318620905UL) + ((uint)op[4] * 130940849UL) + ((uint)op[5] * 2777209316UL) + ((uint)op[6] * 3461746829UL) + ((uint)op[7] * 1825549316UL) + ((uint)op[8] * 32802490UL) + ((uint)op[9] * 2090964479UL) + ((uint)op[10] * 789633811UL) + ((uint)op[11] * 3481490052UL) + ((uint)op[12] * 1049726001UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((llong)tmp_q[0] * 49379L) + ((((llong)tmp_q[1] * 235204L) - ((llong)tmp_q[2] * 340281L) - ((llong)tmp_q[3] * 231035L) - ((llong)tmp_q[4] * 559880L) + ((llong)tmp_q[5] * 76663L) + ((llong)tmp_q[6] * 456205L) + ((llong)tmp_q[7] * 225053L) + ((llong)tmp_q[8] * 100938L) - ((llong)tmp_q[9] * 328187L) + ((llong)tmp_q[10] * 82209L) - ((llong)tmp_q[11] * 484377L) + ((llong)tmp_q[12] * 336486L)) * 2);
	tmp_zero[1] = ((llong)tmp_q[0] * 336486L) - ((llong)tmp_q[1] * 49379L) + ((((llong)tmp_q[2] * 235204L) - ((llong)tmp_q[3] * 340281L) - ((llong)tmp_q[4] * 231035L) - ((llong)tmp_q[5] * 559880L) + ((llong)tmp_q[6] * 76663L) + ((llong)tmp_q[7] * 456205L) + ((llong)tmp_q[8] * 225053L) + ((llong)tmp_q[9] * 100938L) - ((llong)tmp_q[10] * 328187L) + ((llong)tmp_q[11] * 82209L) - ((llong)tmp_q[12] * 484377L)) * 2);
	tmp_zero[2] = -((llong)tmp_q[0] * 484377L) + ((llong)tmp_q[1] * 336486L) - ((llong)tmp_q[2] * 49379L) + ((((llong)tmp_q[3] * 235204L) - ((llong)tmp_q[4] * 340281L) - ((llong)tmp_q[5] * 231035L) - ((llong)tmp_q[6] * 559880L) + ((llong)tmp_q[7] * 76663L) + ((llong)tmp_q[8] * 456205L) + ((llong)tmp_q[9] * 225053L) + ((llong)tmp_q[10] * 100938L) - ((llong)tmp_q[11] * 328187L) + ((llong)tmp_q[12] * 82209L)) * 2);
	tmp_zero[3] = ((llong)tmp_q[0] * 82209L) - ((llong)tmp_q[1] * 484377L) + ((llong)tmp_q[2] * 336486L) - ((llong)tmp_q[3] * 49379L) + ((((llong)tmp_q[4] * 235204L) - ((llong)tmp_q[5] * 340281L) - ((llong)tmp_q[6] * 231035L) - ((llong)tmp_q[7] * 559880L) + ((llong)tmp_q[8] * 76663L) + ((llong)tmp_q[9] * 456205L) + ((llong)tmp_q[10] * 225053L) + ((llong)tmp_q[11] * 100938L) - ((llong)tmp_q[12] * 328187L)) * 2);
	tmp_zero[4] = -((llong)tmp_q[0] * 328187L) + ((llong)tmp_q[1] * 82209L) - ((llong)tmp_q[2] * 484377L) + ((llong)tmp_q[3] * 336486L) - ((llong)tmp_q[4] * 49379L) + ((((llong)tmp_q[5] * 235204L) - ((llong)tmp_q[6] * 340281L) - ((llong)tmp_q[7] * 231035L) - ((llong)tmp_q[8] * 559880L) + ((llong)tmp_q[9] * 76663L) + ((llong)tmp_q[10] * 456205L) + ((llong)tmp_q[11] * 225053L) + ((llong)tmp_q[12] * 100938L)) * 2);
	tmp_zero[5] = ((llong)tmp_q[0] * 100938L) - ((llong)tmp_q[1] * 328187L) + ((llong)tmp_q[2] * 82209L) - ((llong)tmp_q[3] * 484377L) + ((llong)tmp_q[4] * 336486L) - ((llong)tmp_q[5] * 49379L) + ((((llong)tmp_q[6] * 235204L) - ((llong)tmp_q[7] * 340281L) - ((llong)tmp_q[8] * 231035L) - ((llong)tmp_q[9] * 559880L) + ((llong)tmp_q[10] * 76663L) + ((llong)tmp_q[11] * 456205L) + ((llong)tmp_q[12] * 225053L)) * 2);
	tmp_zero[6] = ((llong)tmp_q[0] * 225053L) + ((llong)tmp_q[1] * 100938L) - ((llong)tmp_q[2] * 328187L) + ((llong)tmp_q[3] * 82209L) - ((llong)tmp_q[4] * 484377L) + ((llong)tmp_q[5] * 336486L) - ((llong)tmp_q[6] * 49379L) + ((((llong)tmp_q[7] * 235204L) - ((llong)tmp_q[8] * 340281L) - ((llong)tmp_q[9] * 231035L) - ((llong)tmp_q[10] * 559880L) + ((llong)tmp_q[11] * 76663L) + ((llong)tmp_q[12] * 456205L)) * 2);
	tmp_zero[7] = ((llong)tmp_q[0] * 456205L) + ((llong)tmp_q[1] * 225053L) + ((llong)tmp_q[2] * 100938L) - ((llong)tmp_q[3] * 328187L) + ((llong)tmp_q[4] * 82209L) - ((llong)tmp_q[5] * 484377L) + ((llong)tmp_q[6] * 336486L) - ((llong)tmp_q[7] * 49379L) + ((((llong)tmp_q[8] * 235204L) - ((llong)tmp_q[9] * 340281L) - ((llong)tmp_q[10] * 231035L) - ((llong)tmp_q[11] * 559880L) + ((llong)tmp_q[12] * 76663L)) * 2);
	tmp_zero[8] = ((llong)tmp_q[0] * 76663L) + ((llong)tmp_q[1] * 456205L) + ((llong)tmp_q[2] * 225053L) + ((llong)tmp_q[3] * 100938L) - ((llong)tmp_q[4] * 328187L) + ((llong)tmp_q[5] * 82209L) - ((llong)tmp_q[6] * 484377L) + ((llong)tmp_q[7] * 336486L) - ((llong)tmp_q[8] * 49379L) + ((((llong)tmp_q[9] * 235204L) - ((llong)tmp_q[10] * 340281L) - ((llong)tmp_q[11] * 231035L) - ((llong)tmp_q[12] * 559880L)) * 2);
	tmp_zero[9] = -((llong)tmp_q[0] * 559880L) + ((llong)tmp_q[1] * 76663L) + ((llong)tmp_q[2] * 456205L) + ((llong)tmp_q[3] * 225053L) + ((llong)tmp_q[4] * 100938L) - ((llong)tmp_q[5] * 328187L) + ((llong)tmp_q[6] * 82209L) - ((llong)tmp_q[7] * 484377L) + ((llong)tmp_q[8] * 336486L) - ((llong)tmp_q[9] * 49379L) + ((((llong)tmp_q[10] * 235204L) - ((llong)tmp_q[11] * 340281L) - ((llong)tmp_q[12] * 231035L)) * 2);
	tmp_zero[10] = -((llong)tmp_q[0] * 231035L) - ((llong)tmp_q[1] * 559880L) + ((llong)tmp_q[2] * 76663L) + ((llong)tmp_q[3] * 456205L) + ((llong)tmp_q[4] * 225053L) + ((llong)tmp_q[5] * 100938L) - ((llong)tmp_q[6] * 328187L) + ((llong)tmp_q[7] * 82209L) - ((llong)tmp_q[8] * 484377L) + ((llong)tmp_q[9] * 336486L) - ((llong)tmp_q[10] * 49379L) + ((((llong)tmp_q[11] * 235204L) - ((llong)tmp_q[12] * 340281L)) * 2);
	tmp_zero[11] = -((llong)tmp_q[0] * 340281L) - ((llong)tmp_q[1] * 231035L) - ((llong)tmp_q[2] * 559880L) + ((llong)tmp_q[3] * 76663L) + ((llong)tmp_q[4] * 456205L) + ((llong)tmp_q[5] * 225053L) + ((llong)tmp_q[6] * 100938L) - ((llong)tmp_q[7] * 328187L) + ((llong)tmp_q[8] * 82209L) - ((llong)tmp_q[9] * 484377L) + ((llong)tmp_q[10] * 336486L) - ((llong)tmp_q[11] * 49379L) + ((llong)tmp_q[12] * 470408L);
	tmp_zero[12] = ((llong)tmp_q[0] * 235204L) - ((llong)tmp_q[1] * 340281L) - ((llong)tmp_q[2] * 231035L) - ((llong)tmp_q[3] * 559880L) + ((llong)tmp_q[4] * 76663L) + ((llong)tmp_q[5] * 456205L) + ((llong)tmp_q[6] * 225053L) + ((llong)tmp_q[7] * 100938L) - ((llong)tmp_q[8] * 328187L) + ((llong)tmp_q[9] * 82209L) - ((llong)tmp_q[10] * 484377L) + ((llong)tmp_q[11] * 336486L) - ((llong)tmp_q[12] * 49379L);

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
	rop[12] = (op[12] + tmp_zero[12]) >> WORD_SIZE;
}

