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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14617475749330751917UL) + ((((uint64_t)op[1] * 11208600313900092258UL) + ((uint64_t)op[2] * 1791177471413012274UL) + ((uint64_t)op[3] * 14227300269943098882UL) + ((uint64_t)op[4] * 9062487476388170769UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 9062487476388170769UL) + ((uint64_t)op[1] * 14617475749330751917UL) + ((((uint64_t)op[2] * 11208600313900092258UL) + ((uint64_t)op[3] * 1791177471413012274UL) + ((uint64_t)op[4] * 14227300269943098882UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 14227300269943098882UL) + ((uint64_t)op[1] * 9062487476388170769UL) + ((uint64_t)op[2] * 14617475749330751917UL) + ((((uint64_t)op[3] * 11208600313900092258UL) + ((uint64_t)op[4] * 1791177471413012274UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 1791177471413012274UL) + ((uint64_t)op[1] * 14227300269943098882UL) + ((uint64_t)op[2] * 9062487476388170769UL) + ((uint64_t)op[3] * 14617475749330751917UL) + ((uint64_t)op[4] * 2564917857347020016UL);
	tmp_q[4] = ((uint64_t)op[0] * 11208600313900092258UL) + ((uint64_t)op[1] * 1791177471413012274UL) + ((uint64_t)op[2] * 14227300269943098882UL) + ((uint64_t)op[3] * 9062487476388170769UL) + ((uint64_t)op[4] * 14617475749330751917UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 80245928221L) - ((((int128)tmp_q[1] * 96509310211L) - ((int128)tmp_q[2] * 131937542953L) - ((int128)tmp_q[3] * 167042132515L) + ((int128)tmp_q[4] * 21015798001L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 21015798001L) - ((int128)tmp_q[1] * 80245928221L) - ((((int128)tmp_q[2] * 96509310211L) - ((int128)tmp_q[3] * 131937542953L) - ((int128)tmp_q[4] * 167042132515L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 167042132515L) + ((int128)tmp_q[1] * 21015798001L) - ((int128)tmp_q[2] * 80245928221L) - ((((int128)tmp_q[3] * 96509310211L) - ((int128)tmp_q[4] * 131937542953L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 131937542953L) - ((int128)tmp_q[1] * 167042132515L) + ((int128)tmp_q[2] * 21015798001L) - ((int128)tmp_q[3] * 80245928221L) - ((int128)tmp_q[4] * 772074481688L);
	tmp_zero[4] = ((int128)tmp_q[0] * 96509310211L) - ((int128)tmp_q[1] * 131937542953L) - ((int128)tmp_q[2] * 167042132515L) + ((int128)tmp_q[3] * 21015798001L) - ((int128)tmp_q[4] * 80245928221L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

