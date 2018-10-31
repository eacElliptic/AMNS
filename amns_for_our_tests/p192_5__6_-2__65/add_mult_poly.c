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
	tmp_q[0] = ((uint64_t)op[0] * 4187796153135610239UL) + ((((uint64_t)op[1] * 4921339207909015723UL) + ((uint64_t)op[2] * 8159638907829245933UL) + ((uint64_t)op[3] * 7028692788400141103UL) + ((uint64_t)op[4] * 10878155971202815464UL) + ((uint64_t)op[5] * 14434631647967697219UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 14434631647967697219UL) + ((uint64_t)op[1] * 4187796153135610239UL) + ((((uint64_t)op[2] * 4921339207909015723UL) + ((uint64_t)op[3] * 8159638907829245933UL) + ((uint64_t)op[4] * 7028692788400141103UL) + ((uint64_t)op[5] * 10878155971202815464UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 10878155971202815464UL) + ((uint64_t)op[1] * 14434631647967697219UL) + ((uint64_t)op[2] * 4187796153135610239UL) + ((((uint64_t)op[3] * 4921339207909015723UL) + ((uint64_t)op[4] * 8159638907829245933UL) + ((uint64_t)op[5] * 7028692788400141103UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 7028692788400141103UL) + ((uint64_t)op[1] * 10878155971202815464UL) + ((uint64_t)op[2] * 14434631647967697219UL) + ((uint64_t)op[3] * 4187796153135610239UL) + ((((uint64_t)op[4] * 4921339207909015723UL) + ((uint64_t)op[5] * 8159638907829245933UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 8159638907829245933UL) + ((uint64_t)op[1] * 7028692788400141103UL) + ((uint64_t)op[2] * 10878155971202815464UL) + ((uint64_t)op[3] * 14434631647967697219UL) + ((uint64_t)op[4] * 4187796153135610239UL) + ((uint64_t)op[5] * 8604065657891520170UL);
	tmp_q[5] = ((uint64_t)op[0] * 4921339207909015723UL) + ((uint64_t)op[1] * 8159638907829245933UL) + ((uint64_t)op[2] * 7028692788400141103UL) + ((uint64_t)op[3] * 10878155971202815464UL) + ((uint64_t)op[4] * 14434631647967697219UL) + ((uint64_t)op[5] * 4187796153135610239UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1060127287L) - ((((int128)tmp_q[1] * 497110195L) + ((int128)tmp_q[2] * 734375438L) + ((int128)tmp_q[3] * 425115556L) - ((int128)tmp_q[4] * 223777337L) - ((int128)tmp_q[5] * 3092538161L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3092538161L) + ((int128)tmp_q[1] * 1060127287L) - ((((int128)tmp_q[2] * 497110195L) + ((int128)tmp_q[3] * 734375438L) + ((int128)tmp_q[4] * 425115556L) - ((int128)tmp_q[5] * 223777337L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 223777337L) - ((int128)tmp_q[1] * 3092538161L) + ((int128)tmp_q[2] * 1060127287L) - ((((int128)tmp_q[3] * 497110195L) + ((int128)tmp_q[4] * 734375438L) + ((int128)tmp_q[5] * 425115556L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 425115556L) - ((int128)tmp_q[1] * 223777337L) - ((int128)tmp_q[2] * 3092538161L) + ((int128)tmp_q[3] * 1060127287L) - ((((int128)tmp_q[4] * 497110195L) + ((int128)tmp_q[5] * 734375438L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 734375438L) + ((int128)tmp_q[1] * 425115556L) - ((int128)tmp_q[2] * 223777337L) - ((int128)tmp_q[3] * 3092538161L) + ((int128)tmp_q[4] * 1060127287L) - ((int128)tmp_q[5] * 994220390L);
	tmp_zero[5] = ((int128)tmp_q[0] * 497110195L) + ((int128)tmp_q[1] * 734375438L) + ((int128)tmp_q[2] * 425115556L) - ((int128)tmp_q[3] * 223777337L) - ((int128)tmp_q[4] * 3092538161L) + ((int128)tmp_q[5] * 1060127287L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

