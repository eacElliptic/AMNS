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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3715695496837438259UL) + ((((uint64_t)op[1] * 12102729584031403704UL) + ((uint64_t)op[2] * 15395671585039837904UL) + ((uint64_t)op[3] * 14094686627234054257UL) + ((uint64_t)op[4] * 3071200044047868307UL) + ((uint64_t)op[5] * 6777007530064958100UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 6777007530064958100UL) + ((uint64_t)op[1] * 3715695496837438259UL) + ((((uint64_t)op[2] * 12102729584031403704UL) + ((uint64_t)op[3] * 15395671585039837904UL) + ((uint64_t)op[4] * 14094686627234054257UL) + ((uint64_t)op[5] * 3071200044047868307UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 3071200044047868307UL) + ((uint64_t)op[1] * 6777007530064958100UL) + ((uint64_t)op[2] * 3715695496837438259UL) + ((((uint64_t)op[3] * 12102729584031403704UL) + ((uint64_t)op[4] * 15395671585039837904UL) + ((uint64_t)op[5] * 14094686627234054257UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 14094686627234054257UL) + ((uint64_t)op[1] * 3071200044047868307UL) + ((uint64_t)op[2] * 6777007530064958100UL) + ((uint64_t)op[3] * 3715695496837438259UL) + ((((uint64_t)op[4] * 12102729584031403704UL) + ((uint64_t)op[5] * 15395671585039837904UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 15395671585039837904UL) + ((uint64_t)op[1] * 14094686627234054257UL) + ((uint64_t)op[2] * 3071200044047868307UL) + ((uint64_t)op[3] * 6777007530064958100UL) + ((uint64_t)op[4] * 3715695496837438259UL) + ((uint64_t)op[5] * 4588116303703471552UL);
	tmp_q[5] = ((uint64_t)op[0] * 12102729584031403704UL) + ((uint64_t)op[1] * 15395671585039837904UL) + ((uint64_t)op[2] * 14094686627234054257UL) + ((uint64_t)op[3] * 3071200044047868307UL) + ((uint64_t)op[4] * 6777007530064958100UL) + ((uint64_t)op[5] * 3715695496837438259UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 45061646708475L) + ((((int128)tmp_q[1] * 37360042679546L) + ((int128)tmp_q[2] * 71252819546677L) - ((int128)tmp_q[3] * 66785754315991L) + ((int128)tmp_q[4] * 54496883643611L) + ((int128)tmp_q[5] * 24052866488108L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 24052866488108L) - ((int128)tmp_q[1] * 45061646708475L) + ((((int128)tmp_q[2] * 37360042679546L) + ((int128)tmp_q[3] * 71252819546677L) - ((int128)tmp_q[4] * 66785754315991L) + ((int128)tmp_q[5] * 54496883643611L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 54496883643611L) + ((int128)tmp_q[1] * 24052866488108L) - ((int128)tmp_q[2] * 45061646708475L) + ((((int128)tmp_q[3] * 37360042679546L) + ((int128)tmp_q[4] * 71252819546677L) - ((int128)tmp_q[5] * 66785754315991L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 66785754315991L) + ((int128)tmp_q[1] * 54496883643611L) + ((int128)tmp_q[2] * 24052866488108L) - ((int128)tmp_q[3] * 45061646708475L) + ((((int128)tmp_q[4] * 37360042679546L) + ((int128)tmp_q[5] * 71252819546677L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 71252819546677L) - ((int128)tmp_q[1] * 66785754315991L) + ((int128)tmp_q[2] * 54496883643611L) + ((int128)tmp_q[3] * 24052866488108L) - ((int128)tmp_q[4] * 45061646708475L) + ((int128)tmp_q[5] * 298880341436368L);
	tmp_zero[5] = ((int128)tmp_q[0] * 37360042679546L) + ((int128)tmp_q[1] * 71252819546677L) - ((int128)tmp_q[2] * 66785754315991L) + ((int128)tmp_q[3] * 54496883643611L) + ((int128)tmp_q[4] * 24052866488108L) - ((int128)tmp_q[5] * 45061646708475L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

