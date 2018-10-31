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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17985314699620810593UL) + ((((uint64_t)op[1] * 17641371963274349591UL) + ((uint64_t)op[2] * 4009527718580235143UL) + ((uint64_t)op[3] * 2824393253248838470UL) + ((uint64_t)op[4] * 14482696666755457510UL) + ((uint64_t)op[5] * 3494674485647566529UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 3494674485647566529UL) + ((uint64_t)op[1] * 17985314699620810593UL) + ((((uint64_t)op[2] * 17641371963274349591UL) + ((uint64_t)op[3] * 4009527718580235143UL) + ((uint64_t)op[4] * 2824393253248838470UL) + ((uint64_t)op[5] * 14482696666755457510UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 14482696666755457510UL) + ((uint64_t)op[1] * 3494674485647566529UL) + ((uint64_t)op[2] * 17985314699620810593UL) + ((((uint64_t)op[3] * 17641371963274349591UL) + ((uint64_t)op[4] * 4009527718580235143UL) + ((uint64_t)op[5] * 2824393253248838470UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 2824393253248838470UL) + ((uint64_t)op[1] * 14482696666755457510UL) + ((uint64_t)op[2] * 3494674485647566529UL) + ((uint64_t)op[3] * 17985314699620810593UL) + ((((uint64_t)op[4] * 17641371963274349591UL) + ((uint64_t)op[5] * 4009527718580235143UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 4009527718580235143UL) + ((uint64_t)op[1] * 2824393253248838470UL) + ((uint64_t)op[2] * 14482696666755457510UL) + ((uint64_t)op[3] * 3494674485647566529UL) + ((uint64_t)op[4] * 17985314699620810593UL) + ((uint64_t)op[5] * 16835999852839147566UL);
	tmp_q[5] = ((uint64_t)op[0] * 17641371963274349591UL) + ((uint64_t)op[1] * 4009527718580235143UL) + ((uint64_t)op[2] * 2824393253248838470UL) + ((uint64_t)op[3] * 14482696666755457510UL) + ((uint64_t)op[4] * 3494674485647566529UL) + ((uint64_t)op[5] * 17985314699620810593UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2363881423L) + ((((int128)tmp_q[1] * 224517930L) - ((int128)tmp_q[2] * 324584250L) + ((int128)tmp_q[3] * 2656179597L) + ((int128)tmp_q[4] * 1019240915L) + ((int128)tmp_q[5] * 1871046961L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1871046961L) + ((int128)tmp_q[1] * 2363881423L) + ((((int128)tmp_q[2] * 224517930L) - ((int128)tmp_q[3] * 324584250L) + ((int128)tmp_q[4] * 2656179597L) + ((int128)tmp_q[5] * 1019240915L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1019240915L) + ((int128)tmp_q[1] * 1871046961L) + ((int128)tmp_q[2] * 2363881423L) + ((((int128)tmp_q[3] * 224517930L) - ((int128)tmp_q[4] * 324584250L) + ((int128)tmp_q[5] * 2656179597L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 2656179597L) + ((int128)tmp_q[1] * 1019240915L) + ((int128)tmp_q[2] * 1871046961L) + ((int128)tmp_q[3] * 2363881423L) + ((((int128)tmp_q[4] * 224517930L) - ((int128)tmp_q[5] * 324584250L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 324584250L) + ((int128)tmp_q[1] * 2656179597L) + ((int128)tmp_q[2] * 1019240915L) + ((int128)tmp_q[3] * 1871046961L) + ((int128)tmp_q[4] * 2363881423L) + ((int128)tmp_q[5] * 449035860L);
	tmp_zero[5] = ((int128)tmp_q[0] * 224517930L) - ((int128)tmp_q[1] * 324584250L) + ((int128)tmp_q[2] * 2656179597L) + ((int128)tmp_q[3] * 1019240915L) + ((int128)tmp_q[4] * 1871046961L) + ((int128)tmp_q[5] * 2363881423L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

