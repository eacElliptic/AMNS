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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6013453867113144171UL) + ((((uint64_t)op[1] * 28237443943436938UL) + ((uint64_t)op[2] * 11440255356102298554UL) + ((uint64_t)op[3] * 737028494404468880UL) + ((uint64_t)op[4] * 15422406159952311192UL) + ((uint64_t)op[5] * 11595665265325304226UL) + ((uint64_t)op[6] * 15461418571372666421UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 15461418571372666421UL) + ((uint64_t)op[1] * 6013453867113144171UL) + ((((uint64_t)op[2] * 28237443943436938UL) + ((uint64_t)op[3] * 11440255356102298554UL) + ((uint64_t)op[4] * 737028494404468880UL) + ((uint64_t)op[5] * 15422406159952311192UL) + ((uint64_t)op[6] * 11595665265325304226UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 11595665265325304226UL) + ((uint64_t)op[1] * 15461418571372666421UL) + ((uint64_t)op[2] * 6013453867113144171UL) + ((((uint64_t)op[3] * 28237443943436938UL) + ((uint64_t)op[4] * 11440255356102298554UL) + ((uint64_t)op[5] * 737028494404468880UL) + ((uint64_t)op[6] * 15422406159952311192UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 15422406159952311192UL) + ((uint64_t)op[1] * 11595665265325304226UL) + ((uint64_t)op[2] * 15461418571372666421UL) + ((uint64_t)op[3] * 6013453867113144171UL) + ((((uint64_t)op[4] * 28237443943436938UL) + ((uint64_t)op[5] * 11440255356102298554UL) + ((uint64_t)op[6] * 737028494404468880UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 737028494404468880UL) + ((uint64_t)op[1] * 15422406159952311192UL) + ((uint64_t)op[2] * 11595665265325304226UL) + ((uint64_t)op[3] * 15461418571372666421UL) + ((uint64_t)op[4] * 6013453867113144171UL) + ((((uint64_t)op[5] * 28237443943436938UL) + ((uint64_t)op[6] * 11440255356102298554UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 11440255356102298554UL) + ((uint64_t)op[1] * 737028494404468880UL) + ((uint64_t)op[2] * 15422406159952311192UL) + ((uint64_t)op[3] * 11595665265325304226UL) + ((uint64_t)op[4] * 15461418571372666421UL) + ((uint64_t)op[5] * 6013453867113144171UL) + ((uint64_t)op[6] * 112949775773747752UL);
	tmp_q[6] = ((uint64_t)op[0] * 28237443943436938UL) + ((uint64_t)op[1] * 11440255356102298554UL) + ((uint64_t)op[2] * 737028494404468880UL) + ((uint64_t)op[3] * 15422406159952311192UL) + ((uint64_t)op[4] * 11595665265325304226UL) + ((uint64_t)op[5] * 15461418571372666421UL) + ((uint64_t)op[6] * 6013453867113144171UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 56994042727L) + ((((int128)tmp_q[1] * 57090534785L) - ((int128)tmp_q[2] * 43754604921L) + ((int128)tmp_q[3] * 10562419859L) + ((int128)tmp_q[4] * 10957444941L) + ((int128)tmp_q[5] * 54500245899L) + ((int128)tmp_q[6] * 12554064001L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 12554064001L) - ((int128)tmp_q[1] * 56994042727L) + ((((int128)tmp_q[2] * 57090534785L) - ((int128)tmp_q[3] * 43754604921L) + ((int128)tmp_q[4] * 10562419859L) + ((int128)tmp_q[5] * 10957444941L) + ((int128)tmp_q[6] * 54500245899L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 54500245899L) + ((int128)tmp_q[1] * 12554064001L) - ((int128)tmp_q[2] * 56994042727L) + ((((int128)tmp_q[3] * 57090534785L) - ((int128)tmp_q[4] * 43754604921L) + ((int128)tmp_q[5] * 10562419859L) + ((int128)tmp_q[6] * 10957444941L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 10957444941L) + ((int128)tmp_q[1] * 54500245899L) + ((int128)tmp_q[2] * 12554064001L) - ((int128)tmp_q[3] * 56994042727L) + ((((int128)tmp_q[4] * 57090534785L) - ((int128)tmp_q[5] * 43754604921L) + ((int128)tmp_q[6] * 10562419859L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 10562419859L) + ((int128)tmp_q[1] * 10957444941L) + ((int128)tmp_q[2] * 54500245899L) + ((int128)tmp_q[3] * 12554064001L) - ((int128)tmp_q[4] * 56994042727L) + ((((int128)tmp_q[5] * 57090534785L) - ((int128)tmp_q[6] * 43754604921L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 43754604921L) + ((int128)tmp_q[1] * 10562419859L) + ((int128)tmp_q[2] * 10957444941L) + ((int128)tmp_q[3] * 54500245899L) + ((int128)tmp_q[4] * 12554064001L) - ((int128)tmp_q[5] * 56994042727L) + ((int128)tmp_q[6] * 228362139140L);
	tmp_zero[6] = ((int128)tmp_q[0] * 57090534785L) - ((int128)tmp_q[1] * 43754604921L) + ((int128)tmp_q[2] * 10562419859L) + ((int128)tmp_q[3] * 10957444941L) + ((int128)tmp_q[4] * 54500245899L) + ((int128)tmp_q[5] * 12554064001L) - ((int128)tmp_q[6] * 56994042727L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

