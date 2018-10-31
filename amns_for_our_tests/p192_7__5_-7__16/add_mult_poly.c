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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 541852517589249870UL) + ((((uint64_t)op[1] * 15567722612426336052UL) + ((uint64_t)op[2] * 291969749387011576UL) + ((uint64_t)op[3] * 17374336460169383582UL) + ((uint64_t)op[4] * 17558063862361670059UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 17558063862361670059UL) + ((uint64_t)op[1] * 541852517589249870UL) + ((((uint64_t)op[2] * 15567722612426336052UL) + ((uint64_t)op[3] * 291969749387011576UL) + ((uint64_t)op[4] * 17374336460169383582UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 17374336460169383582UL) + ((uint64_t)op[1] * 17558063862361670059UL) + ((uint64_t)op[2] * 541852517589249870UL) + ((((uint64_t)op[3] * 15567722612426336052UL) + ((uint64_t)op[4] * 291969749387011576UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 291969749387011576UL) + ((uint64_t)op[1] * 17374336460169383582UL) + ((uint64_t)op[2] * 17558063862361670059UL) + ((uint64_t)op[3] * 541852517589249870UL) + ((uint64_t)op[4] * 1706406155272957332UL);
	tmp_q[4] = ((uint64_t)op[0] * 15567722612426336052UL) + ((uint64_t)op[1] * 291969749387011576UL) + ((uint64_t)op[2] * 17374336460169383582UL) + ((uint64_t)op[3] * 17558063862361670059UL) + ((uint64_t)op[4] * 541852517589249870UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 102555668506L) - ((((int128)tmp_q[1] * 103340319197L) + ((int128)tmp_q[2] * 135991102070L) + ((int128)tmp_q[3] * 178398770976L) + ((int128)tmp_q[4] * 199705369396L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 199705369396L) - ((int128)tmp_q[1] * 102555668506L) - ((((int128)tmp_q[2] * 103340319197L) + ((int128)tmp_q[3] * 135991102070L) + ((int128)tmp_q[4] * 178398770976L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 178398770976L) + ((int128)tmp_q[1] * 199705369396L) - ((int128)tmp_q[2] * 102555668506L) - ((((int128)tmp_q[3] * 103340319197L) + ((int128)tmp_q[4] * 135991102070L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 135991102070L) + ((int128)tmp_q[1] * 178398770976L) + ((int128)tmp_q[2] * 199705369396L) - ((int128)tmp_q[3] * 102555668506L) - ((int128)tmp_q[4] * 723382234379L);
	tmp_zero[4] = ((int128)tmp_q[0] * 103340319197L) + ((int128)tmp_q[1] * 135991102070L) + ((int128)tmp_q[2] * 178398770976L) + ((int128)tmp_q[3] * 199705369396L) - ((int128)tmp_q[4] * 102555668506L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

