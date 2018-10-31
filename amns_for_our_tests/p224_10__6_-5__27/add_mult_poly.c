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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8424133393911328010UL) + ((((uint64_t)op[1] * 14031588572060339720UL) + ((uint64_t)op[2] * 5140516678921490265UL) + ((uint64_t)op[3] * 2443914767955133598UL) + ((uint64_t)op[4] * 1578648027859839330UL) + ((uint64_t)op[5] * 17500498346195573178UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 17500498346195573178UL) + ((uint64_t)op[1] * 8424133393911328010UL) + ((((uint64_t)op[2] * 14031588572060339720UL) + ((uint64_t)op[3] * 5140516678921490265UL) + ((uint64_t)op[4] * 2443914767955133598UL) + ((uint64_t)op[5] * 1578648027859839330UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1578648027859839330UL) + ((uint64_t)op[1] * 17500498346195573178UL) + ((uint64_t)op[2] * 8424133393911328010UL) + ((((uint64_t)op[3] * 14031588572060339720UL) + ((uint64_t)op[4] * 5140516678921490265UL) + ((uint64_t)op[5] * 2443914767955133598UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 2443914767955133598UL) + ((uint64_t)op[1] * 1578648027859839330UL) + ((uint64_t)op[2] * 17500498346195573178UL) + ((uint64_t)op[3] * 8424133393911328010UL) + ((((uint64_t)op[4] * 14031588572060339720UL) + ((uint64_t)op[5] * 5140516678921490265UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 5140516678921490265UL) + ((uint64_t)op[1] * 2443914767955133598UL) + ((uint64_t)op[2] * 1578648027859839330UL) + ((uint64_t)op[3] * 17500498346195573178UL) + ((uint64_t)op[4] * 8424133393911328010UL) + ((uint64_t)op[5] * 3629033434536507864UL);
	tmp_q[5] = ((uint64_t)op[0] * 14031588572060339720UL) + ((uint64_t)op[1] * 5140516678921490265UL) + ((uint64_t)op[2] * 2443914767955133598UL) + ((uint64_t)op[3] * 1578648027859839330UL) + ((uint64_t)op[4] * 17500498346195573178UL) + ((uint64_t)op[5] * 8424133393911328010UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 68994029734L) - ((-((int128)tmp_q[1] * 69718552150L) + ((int128)tmp_q[2] * 12239787334L) - ((int128)tmp_q[3] * 112466272096L) - ((int128)tmp_q[4] * 154542648159L) - ((int128)tmp_q[5] * 11949998590L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 11949998590L) + ((int128)tmp_q[1] * 68994029734L) - ((-((int128)tmp_q[2] * 69718552150L) + ((int128)tmp_q[3] * 12239787334L) - ((int128)tmp_q[4] * 112466272096L) - ((int128)tmp_q[5] * 154542648159L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 154542648159L) - ((int128)tmp_q[1] * 11949998590L) + ((int128)tmp_q[2] * 68994029734L) - ((-((int128)tmp_q[3] * 69718552150L) + ((int128)tmp_q[4] * 12239787334L) - ((int128)tmp_q[5] * 112466272096L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 112466272096L) - ((int128)tmp_q[1] * 154542648159L) - ((int128)tmp_q[2] * 11949998590L) + ((int128)tmp_q[3] * 68994029734L) - ((-((int128)tmp_q[4] * 69718552150L) + ((int128)tmp_q[5] * 12239787334L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 12239787334L) - ((int128)tmp_q[1] * 112466272096L) - ((int128)tmp_q[2] * 154542648159L) - ((int128)tmp_q[3] * 11949998590L) + ((int128)tmp_q[4] * 68994029734L) + ((int128)tmp_q[5] * 348592760750L);
	tmp_zero[5] = -((int128)tmp_q[0] * 69718552150L) + ((int128)tmp_q[1] * 12239787334L) - ((int128)tmp_q[2] * 112466272096L) - ((int128)tmp_q[3] * 154542648159L) - ((int128)tmp_q[4] * 11949998590L) + ((int128)tmp_q[5] * 68994029734L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

