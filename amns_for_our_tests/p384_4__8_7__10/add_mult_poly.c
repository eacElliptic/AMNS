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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16293154995155716429UL) + ((((uint64_t)op[1] * 13533995923480152711UL) + ((uint64_t)op[2] * 12381127729791310027UL) + ((uint64_t)op[3] * 4880860500674439569UL) + ((uint64_t)op[4] * 17317563474435168755UL) + ((uint64_t)op[5] * 9452832313073898122UL) + ((uint64_t)op[6] * 8395845340034823981UL) + ((uint64_t)op[7] * 7020775093386178821UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 7020775093386178821UL) + ((uint64_t)op[1] * 16293154995155716429UL) + ((((uint64_t)op[2] * 13533995923480152711UL) + ((uint64_t)op[3] * 12381127729791310027UL) + ((uint64_t)op[4] * 4880860500674439569UL) + ((uint64_t)op[5] * 17317563474435168755UL) + ((uint64_t)op[6] * 9452832313073898122UL) + ((uint64_t)op[7] * 8395845340034823981UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 8395845340034823981UL) + ((uint64_t)op[1] * 7020775093386178821UL) + ((uint64_t)op[2] * 16293154995155716429UL) + ((((uint64_t)op[3] * 13533995923480152711UL) + ((uint64_t)op[4] * 12381127729791310027UL) + ((uint64_t)op[5] * 4880860500674439569UL) + ((uint64_t)op[6] * 17317563474435168755UL) + ((uint64_t)op[7] * 9452832313073898122UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 9452832313073898122UL) + ((uint64_t)op[1] * 8395845340034823981UL) + ((uint64_t)op[2] * 7020775093386178821UL) + ((uint64_t)op[3] * 16293154995155716429UL) + ((((uint64_t)op[4] * 13533995923480152711UL) + ((uint64_t)op[5] * 12381127729791310027UL) + ((uint64_t)op[6] * 4880860500674439569UL) + ((uint64_t)op[7] * 17317563474435168755UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 17317563474435168755UL) + ((uint64_t)op[1] * 9452832313073898122UL) + ((uint64_t)op[2] * 8395845340034823981UL) + ((uint64_t)op[3] * 7020775093386178821UL) + ((uint64_t)op[4] * 16293154995155716429UL) + ((((uint64_t)op[5] * 13533995923480152711UL) + ((uint64_t)op[6] * 12381127729791310027UL) + ((uint64_t)op[7] * 4880860500674439569UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 4880860500674439569UL) + ((uint64_t)op[1] * 17317563474435168755UL) + ((uint64_t)op[2] * 9452832313073898122UL) + ((uint64_t)op[3] * 8395845340034823981UL) + ((uint64_t)op[4] * 7020775093386178821UL) + ((uint64_t)op[5] * 16293154995155716429UL) + ((((uint64_t)op[6] * 13533995923480152711UL) + ((uint64_t)op[7] * 12381127729791310027UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 12381127729791310027UL) + ((uint64_t)op[1] * 4880860500674439569UL) + ((uint64_t)op[2] * 17317563474435168755UL) + ((uint64_t)op[3] * 9452832313073898122UL) + ((uint64_t)op[4] * 8395845340034823981UL) + ((uint64_t)op[5] * 7020775093386178821UL) + ((uint64_t)op[6] * 16293154995155716429UL) + ((uint64_t)op[7] * 2504251095813310897UL);
	tmp_q[7] = ((uint64_t)op[0] * 13533995923480152711UL) + ((uint64_t)op[1] * 12381127729791310027UL) + ((uint64_t)op[2] * 4880860500674439569UL) + ((uint64_t)op[3] * 17317563474435168755UL) + ((uint64_t)op[4] * 9452832313073898122UL) + ((uint64_t)op[5] * 8395845340034823981UL) + ((uint64_t)op[6] * 7020775093386178821UL) + ((uint64_t)op[7] * 16293154995155716429UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 88795344353371L) + ((((int128)tmp_q[1] * 55318029920221L) + ((int128)tmp_q[2] * 113690752243981L) - ((int128)tmp_q[3] * 53661519766182L) + ((int128)tmp_q[4] * 146656173500179L) - ((int128)tmp_q[5] * 22969078315715L) - ((int128)tmp_q[6] * 93964571119661L) + ((int128)tmp_q[7] * 112014621705101L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 112014621705101L) + ((int128)tmp_q[1] * 88795344353371L) + ((((int128)tmp_q[2] * 55318029920221L) + ((int128)tmp_q[3] * 113690752243981L) - ((int128)tmp_q[4] * 53661519766182L) + ((int128)tmp_q[5] * 146656173500179L) - ((int128)tmp_q[6] * 22969078315715L) - ((int128)tmp_q[7] * 93964571119661L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 93964571119661L) + ((int128)tmp_q[1] * 112014621705101L) + ((int128)tmp_q[2] * 88795344353371L) + ((((int128)tmp_q[3] * 55318029920221L) + ((int128)tmp_q[4] * 113690752243981L) - ((int128)tmp_q[5] * 53661519766182L) + ((int128)tmp_q[6] * 146656173500179L) - ((int128)tmp_q[7] * 22969078315715L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 22969078315715L) - ((int128)tmp_q[1] * 93964571119661L) + ((int128)tmp_q[2] * 112014621705101L) + ((int128)tmp_q[3] * 88795344353371L) + ((((int128)tmp_q[4] * 55318029920221L) + ((int128)tmp_q[5] * 113690752243981L) - ((int128)tmp_q[6] * 53661519766182L) + ((int128)tmp_q[7] * 146656173500179L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 146656173500179L) - ((int128)tmp_q[1] * 22969078315715L) - ((int128)tmp_q[2] * 93964571119661L) + ((int128)tmp_q[3] * 112014621705101L) + ((int128)tmp_q[4] * 88795344353371L) + ((((int128)tmp_q[5] * 55318029920221L) + ((int128)tmp_q[6] * 113690752243981L) - ((int128)tmp_q[7] * 53661519766182L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 53661519766182L) + ((int128)tmp_q[1] * 146656173500179L) - ((int128)tmp_q[2] * 22969078315715L) - ((int128)tmp_q[3] * 93964571119661L) + ((int128)tmp_q[4] * 112014621705101L) + ((int128)tmp_q[5] * 88795344353371L) + ((((int128)tmp_q[6] * 55318029920221L) + ((int128)tmp_q[7] * 113690752243981L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 113690752243981L) - ((int128)tmp_q[1] * 53661519766182L) + ((int128)tmp_q[2] * 146656173500179L) - ((int128)tmp_q[3] * 22969078315715L) - ((int128)tmp_q[4] * 93964571119661L) + ((int128)tmp_q[5] * 112014621705101L) + ((int128)tmp_q[6] * 88795344353371L) + ((int128)tmp_q[7] * 387226209441547L);
	tmp_zero[7] = ((int128)tmp_q[0] * 55318029920221L) + ((int128)tmp_q[1] * 113690752243981L) - ((int128)tmp_q[2] * 53661519766182L) + ((int128)tmp_q[3] * 146656173500179L) - ((int128)tmp_q[4] * 22969078315715L) - ((int128)tmp_q[5] * 93964571119661L) + ((int128)tmp_q[6] * 112014621705101L) + ((int128)tmp_q[7] * 88795344353371L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

