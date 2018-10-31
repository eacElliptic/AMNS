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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4100389129211510281UL) + ((((uint64_t)op[1] * 7352338010449149434UL) + ((uint64_t)op[2] * 2248293533078716305UL) + ((uint64_t)op[3] * 15061395322929128710UL) + ((uint64_t)op[4] * 17031578839484112515UL) + ((uint64_t)op[5] * 13734516297553755886UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 13734516297553755886UL) + ((uint64_t)op[1] * 4100389129211510281UL) + ((((uint64_t)op[2] * 7352338010449149434UL) + ((uint64_t)op[3] * 2248293533078716305UL) + ((uint64_t)op[4] * 15061395322929128710UL) + ((uint64_t)op[5] * 17031578839484112515UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 17031578839484112515UL) + ((uint64_t)op[1] * 13734516297553755886UL) + ((uint64_t)op[2] * 4100389129211510281UL) + ((((uint64_t)op[3] * 7352338010449149434UL) + ((uint64_t)op[4] * 2248293533078716305UL) + ((uint64_t)op[5] * 15061395322929128710UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 15061395322929128710UL) + ((uint64_t)op[1] * 17031578839484112515UL) + ((uint64_t)op[2] * 13734516297553755886UL) + ((uint64_t)op[3] * 4100389129211510281UL) + ((((uint64_t)op[4] * 7352338010449149434UL) + ((uint64_t)op[5] * 2248293533078716305UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 2248293533078716305UL) + ((uint64_t)op[1] * 15061395322929128710UL) + ((uint64_t)op[2] * 17031578839484112515UL) + ((uint64_t)op[3] * 13734516297553755886UL) + ((uint64_t)op[4] * 4100389129211510281UL) + ((uint64_t)op[5] * 3742068052811252748UL);
	tmp_q[5] = ((uint64_t)op[0] * 7352338010449149434UL) + ((uint64_t)op[1] * 2248293533078716305UL) + ((uint64_t)op[2] * 15061395322929128710UL) + ((uint64_t)op[3] * 17031578839484112515UL) + ((uint64_t)op[4] * 13734516297553755886UL) + ((uint64_t)op[5] * 4100389129211510281UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1197068609L) - ((((int128)tmp_q[1] * 677482988L) + ((int128)tmp_q[2] * 46850884L) - ((int128)tmp_q[3] * 2269731546L) - ((int128)tmp_q[4] * 103851623L) - ((int128)tmp_q[5] * 2446964390L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 2446964390L) + ((int128)tmp_q[1] * 1197068609L) - ((((int128)tmp_q[2] * 677482988L) + ((int128)tmp_q[3] * 46850884L) - ((int128)tmp_q[4] * 2269731546L) - ((int128)tmp_q[5] * 103851623L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 103851623L) - ((int128)tmp_q[1] * 2446964390L) + ((int128)tmp_q[2] * 1197068609L) - ((((int128)tmp_q[3] * 677482988L) + ((int128)tmp_q[4] * 46850884L) - ((int128)tmp_q[5] * 2269731546L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2269731546L) - ((int128)tmp_q[1] * 103851623L) - ((int128)tmp_q[2] * 2446964390L) + ((int128)tmp_q[3] * 1197068609L) - ((((int128)tmp_q[4] * 677482988L) + ((int128)tmp_q[5] * 46850884L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 46850884L) - ((int128)tmp_q[1] * 2269731546L) - ((int128)tmp_q[2] * 103851623L) - ((int128)tmp_q[3] * 2446964390L) + ((int128)tmp_q[4] * 1197068609L) - ((int128)tmp_q[5] * 1354965976L);
	tmp_zero[5] = ((int128)tmp_q[0] * 677482988L) + ((int128)tmp_q[1] * 46850884L) - ((int128)tmp_q[2] * 2269731546L) - ((int128)tmp_q[3] * 103851623L) - ((int128)tmp_q[4] * 2446964390L) + ((int128)tmp_q[5] * 1197068609L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

