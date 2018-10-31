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
	tmp_q[0] = ((uint64_t)op[0] * 8475221856126237799UL) + ((((uint64_t)op[1] * 13174419169507106173UL) + ((uint64_t)op[2] * 16480908747901257306UL) + ((uint64_t)op[3] * 2347648242064620053UL) + ((uint64_t)op[4] * 3308954525266821762UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 3308954525266821762UL) + ((uint64_t)op[1] * 8475221856126237799UL) + ((((uint64_t)op[2] * 13174419169507106173UL) + ((uint64_t)op[3] * 16480908747901257306UL) + ((uint64_t)op[4] * 2347648242064620053UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 2347648242064620053UL) + ((uint64_t)op[1] * 3308954525266821762UL) + ((uint64_t)op[2] * 8475221856126237799UL) + ((((uint64_t)op[3] * 13174419169507106173UL) + ((uint64_t)op[4] * 16480908747901257306UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 16480908747901257306UL) + ((uint64_t)op[1] * 2347648242064620053UL) + ((uint64_t)op[2] * 3308954525266821762UL) + ((uint64_t)op[3] * 8475221856126237799UL) + ((uint64_t)op[4] * 12786181998014869UL);
	tmp_q[4] = ((uint64_t)op[0] * 13174419169507106173UL) + ((uint64_t)op[1] * 16480908747901257306UL) + ((uint64_t)op[2] * 2347648242064620053UL) + ((uint64_t)op[3] * 3308954525266821762UL) + ((uint64_t)op[4] * 8475221856126237799UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2858425367224L) - ((((int128)tmp_q[1] * 10725271849901L) + ((int128)tmp_q[2] * 18278133149617L) - ((int128)tmp_q[3] * 12682860799173L) + ((int128)tmp_q[4] * 4510642162954L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 4510642162954L) - ((int128)tmp_q[1] * 2858425367224L) - ((((int128)tmp_q[2] * 10725271849901L) + ((int128)tmp_q[3] * 18278133149617L) - ((int128)tmp_q[4] * 12682860799173L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 12682860799173L) + ((int128)tmp_q[1] * 4510642162954L) - ((int128)tmp_q[2] * 2858425367224L) - ((((int128)tmp_q[3] * 10725271849901L) + ((int128)tmp_q[4] * 18278133149617L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 18278133149617L) - ((int128)tmp_q[1] * 12682860799173L) + ((int128)tmp_q[2] * 4510642162954L) - ((int128)tmp_q[3] * 2858425367224L) - ((int128)tmp_q[4] * 75076902949307L);
	tmp_zero[4] = ((int128)tmp_q[0] * 10725271849901L) + ((int128)tmp_q[1] * 18278133149617L) - ((int128)tmp_q[2] * 12682860799173L) + ((int128)tmp_q[3] * 4510642162954L) - ((int128)tmp_q[4] * 2858425367224L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

