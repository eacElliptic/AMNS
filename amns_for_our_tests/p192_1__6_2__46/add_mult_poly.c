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
	tmp_q[0] = ((uint64_t)op[0] * 12245048511711231757UL) + ((((uint64_t)op[1] * 5928236917725341222UL) + ((uint64_t)op[2] * 2120867970304955432UL) + ((uint64_t)op[3] * 9840501184838633912UL) + ((uint64_t)op[4] * 3386043308364440882UL) + ((uint64_t)op[5] * 4236402436937045727UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 4236402436937045727UL) + ((uint64_t)op[1] * 12245048511711231757UL) + ((((uint64_t)op[2] * 5928236917725341222UL) + ((uint64_t)op[3] * 2120867970304955432UL) + ((uint64_t)op[4] * 9840501184838633912UL) + ((uint64_t)op[5] * 3386043308364440882UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 3386043308364440882UL) + ((uint64_t)op[1] * 4236402436937045727UL) + ((uint64_t)op[2] * 12245048511711231757UL) + ((((uint64_t)op[3] * 5928236917725341222UL) + ((uint64_t)op[4] * 2120867970304955432UL) + ((uint64_t)op[5] * 9840501184838633912UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 9840501184838633912UL) + ((uint64_t)op[1] * 3386043308364440882UL) + ((uint64_t)op[2] * 4236402436937045727UL) + ((uint64_t)op[3] * 12245048511711231757UL) + ((((uint64_t)op[4] * 5928236917725341222UL) + ((uint64_t)op[5] * 2120867970304955432UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 2120867970304955432UL) + ((uint64_t)op[1] * 9840501184838633912UL) + ((uint64_t)op[2] * 3386043308364440882UL) + ((uint64_t)op[3] * 4236402436937045727UL) + ((uint64_t)op[4] * 12245048511711231757UL) + ((uint64_t)op[5] * 11856473835450682444UL);
	tmp_q[5] = ((uint64_t)op[0] * 5928236917725341222UL) + ((uint64_t)op[1] * 2120867970304955432UL) + ((uint64_t)op[2] * 9840501184838633912UL) + ((uint64_t)op[3] * 3386043308364440882UL) + ((uint64_t)op[4] * 4236402436937045727UL) + ((uint64_t)op[5] * 12245048511711231757UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1043713591L) + ((((int128)tmp_q[1] * 1784341407L) + ((int128)tmp_q[2] * 230404123L) + ((int128)tmp_q[3] * 1449602745L) + ((int128)tmp_q[4] * 1937439419L) + ((int128)tmp_q[5] * 6069773L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 6069773L) - ((int128)tmp_q[1] * 1043713591L) + ((((int128)tmp_q[2] * 1784341407L) + ((int128)tmp_q[3] * 230404123L) + ((int128)tmp_q[4] * 1449602745L) + ((int128)tmp_q[5] * 1937439419L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1937439419L) + ((int128)tmp_q[1] * 6069773L) - ((int128)tmp_q[2] * 1043713591L) + ((((int128)tmp_q[3] * 1784341407L) + ((int128)tmp_q[4] * 230404123L) + ((int128)tmp_q[5] * 1449602745L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 1449602745L) + ((int128)tmp_q[1] * 1937439419L) + ((int128)tmp_q[2] * 6069773L) - ((int128)tmp_q[3] * 1043713591L) + ((((int128)tmp_q[4] * 1784341407L) + ((int128)tmp_q[5] * 230404123L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 230404123L) + ((int128)tmp_q[1] * 1449602745L) + ((int128)tmp_q[2] * 1937439419L) + ((int128)tmp_q[3] * 6069773L) - ((int128)tmp_q[4] * 1043713591L) + ((int128)tmp_q[5] * 3568682814L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1784341407L) + ((int128)tmp_q[1] * 230404123L) + ((int128)tmp_q[2] * 1449602745L) + ((int128)tmp_q[3] * 1937439419L) + ((int128)tmp_q[4] * 6069773L) - ((int128)tmp_q[5] * 1043713591L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

