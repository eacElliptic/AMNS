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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + ((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + ((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + ((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + ((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + ((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + ((int128)pa[6] * pb[6]);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + ((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + ((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + ((int128)pa[6] * pa[6]);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2968790149170036589UL) + ((uint64_t)op[1] * 4095746028037246517UL) + ((uint64_t)op[2] * 12204949823985932922UL) + ((uint64_t)op[3] * 4138964285347794940UL) + ((uint64_t)op[4] * 12396601816369395832UL) + ((uint64_t)op[5] * 433829919101686902UL) + ((uint64_t)op[6] * 654606125407009529UL);
	tmp_q[1] = ((uint64_t)op[0] * 654606125407009529UL) + ((uint64_t)op[1] * 2968790149170036589UL) + ((uint64_t)op[2] * 4095746028037246517UL) + ((uint64_t)op[3] * 12204949823985932922UL) + ((uint64_t)op[4] * 4138964285347794940UL) + ((uint64_t)op[5] * 12396601816369395832UL) + ((uint64_t)op[6] * 433829919101686902UL);
	tmp_q[2] = ((uint64_t)op[0] * 433829919101686902UL) + ((uint64_t)op[1] * 654606125407009529UL) + ((uint64_t)op[2] * 2968790149170036589UL) + ((uint64_t)op[3] * 4095746028037246517UL) + ((uint64_t)op[4] * 12204949823985932922UL) + ((uint64_t)op[5] * 4138964285347794940UL) + ((uint64_t)op[6] * 12396601816369395832UL);
	tmp_q[3] = ((uint64_t)op[0] * 12396601816369395832UL) + ((uint64_t)op[1] * 433829919101686902UL) + ((uint64_t)op[2] * 654606125407009529UL) + ((uint64_t)op[3] * 2968790149170036589UL) + ((uint64_t)op[4] * 4095746028037246517UL) + ((uint64_t)op[5] * 12204949823985932922UL) + ((uint64_t)op[6] * 4138964285347794940UL);
	tmp_q[4] = ((uint64_t)op[0] * 4138964285347794940UL) + ((uint64_t)op[1] * 12396601816369395832UL) + ((uint64_t)op[2] * 433829919101686902UL) + ((uint64_t)op[3] * 654606125407009529UL) + ((uint64_t)op[4] * 2968790149170036589UL) + ((uint64_t)op[5] * 4095746028037246517UL) + ((uint64_t)op[6] * 12204949823985932922UL);
	tmp_q[5] = ((uint64_t)op[0] * 12204949823985932922UL) + ((uint64_t)op[1] * 4138964285347794940UL) + ((uint64_t)op[2] * 12396601816369395832UL) + ((uint64_t)op[3] * 433829919101686902UL) + ((uint64_t)op[4] * 654606125407009529UL) + ((uint64_t)op[5] * 2968790149170036589UL) + ((uint64_t)op[6] * 4095746028037246517UL);
	tmp_q[6] = ((uint64_t)op[0] * 4095746028037246517UL) + ((uint64_t)op[1] * 12204949823985932922UL) + ((uint64_t)op[2] * 4138964285347794940UL) + ((uint64_t)op[3] * 12396601816369395832UL) + ((uint64_t)op[4] * 433829919101686902UL) + ((uint64_t)op[5] * 654606125407009529UL) + ((uint64_t)op[6] * 2968790149170036589UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2820509949243L) + (-((int128)tmp_q[1] * 678162471525L) + ((int128)tmp_q[2] * 1882458016924L) - ((int128)tmp_q[3] * 965966272159L) + ((int128)tmp_q[4] * 177460504869L) + ((int128)tmp_q[5] * 4782079300980L) - ((int128)tmp_q[6] * 2377359129845L));
	tmp_zero[1] = -((int128)tmp_q[0] * 2377359129845L) - ((int128)tmp_q[1] * 2820509949243L) + (-((int128)tmp_q[2] * 678162471525L) + ((int128)tmp_q[3] * 1882458016924L) - ((int128)tmp_q[4] * 965966272159L) + ((int128)tmp_q[5] * 177460504869L) + ((int128)tmp_q[6] * 4782079300980L));
	tmp_zero[2] = ((int128)tmp_q[0] * 4782079300980L) - ((int128)tmp_q[1] * 2377359129845L) - ((int128)tmp_q[2] * 2820509949243L) + (-((int128)tmp_q[3] * 678162471525L) + ((int128)tmp_q[4] * 1882458016924L) - ((int128)tmp_q[5] * 965966272159L) + ((int128)tmp_q[6] * 177460504869L));
	tmp_zero[3] = ((int128)tmp_q[0] * 177460504869L) + ((int128)tmp_q[1] * 4782079300980L) - ((int128)tmp_q[2] * 2377359129845L) - ((int128)tmp_q[3] * 2820509949243L) + (-((int128)tmp_q[4] * 678162471525L) + ((int128)tmp_q[5] * 1882458016924L) - ((int128)tmp_q[6] * 965966272159L));
	tmp_zero[4] = -((int128)tmp_q[0] * 965966272159L) + ((int128)tmp_q[1] * 177460504869L) + ((int128)tmp_q[2] * 4782079300980L) - ((int128)tmp_q[3] * 2377359129845L) - ((int128)tmp_q[4] * 2820509949243L) + (-((int128)tmp_q[5] * 678162471525L) + ((int128)tmp_q[6] * 1882458016924L));
	tmp_zero[5] = ((int128)tmp_q[0] * 1882458016924L) - ((int128)tmp_q[1] * 965966272159L) + ((int128)tmp_q[2] * 177460504869L) + ((int128)tmp_q[3] * 4782079300980L) - ((int128)tmp_q[4] * 2377359129845L) - ((int128)tmp_q[5] * 2820509949243L) - ((int128)tmp_q[6] * 678162471525L);
	tmp_zero[6] = -((int128)tmp_q[0] * 678162471525L) + ((int128)tmp_q[1] * 1882458016924L) - ((int128)tmp_q[2] * 965966272159L) + ((int128)tmp_q[3] * 177460504869L) + ((int128)tmp_q[4] * 4782079300980L) - ((int128)tmp_q[5] * 2377359129845L) - ((int128)tmp_q[6] * 2820509949243L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

