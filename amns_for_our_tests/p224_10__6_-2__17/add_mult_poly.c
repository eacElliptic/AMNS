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
	tmp_q[0] = ((uint64_t)op[0] * 17740737518988789911UL) + ((((uint64_t)op[1] * 6786653434676270088UL) + ((uint64_t)op[2] * 15025569982877118022UL) + ((uint64_t)op[3] * 16135997692946310157UL) + ((uint64_t)op[4] * 2412456884085593928UL) + ((uint64_t)op[5] * 10102378073999687720UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 10102378073999687720UL) + ((uint64_t)op[1] * 17740737518988789911UL) + ((((uint64_t)op[2] * 6786653434676270088UL) + ((uint64_t)op[3] * 15025569982877118022UL) + ((uint64_t)op[4] * 16135997692946310157UL) + ((uint64_t)op[5] * 2412456884085593928UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2412456884085593928UL) + ((uint64_t)op[1] * 10102378073999687720UL) + ((uint64_t)op[2] * 17740737518988789911UL) + ((((uint64_t)op[3] * 6786653434676270088UL) + ((uint64_t)op[4] * 15025569982877118022UL) + ((uint64_t)op[5] * 16135997692946310157UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 16135997692946310157UL) + ((uint64_t)op[1] * 2412456884085593928UL) + ((uint64_t)op[2] * 10102378073999687720UL) + ((uint64_t)op[3] * 17740737518988789911UL) + ((((uint64_t)op[4] * 6786653434676270088UL) + ((uint64_t)op[5] * 15025569982877118022UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 15025569982877118022UL) + ((uint64_t)op[1] * 16135997692946310157UL) + ((uint64_t)op[2] * 2412456884085593928UL) + ((uint64_t)op[3] * 10102378073999687720UL) + ((uint64_t)op[4] * 17740737518988789911UL) + ((uint64_t)op[5] * 4873437204357011440UL);
	tmp_q[5] = ((uint64_t)op[0] * 6786653434676270088UL) + ((uint64_t)op[1] * 15025569982877118022UL) + ((uint64_t)op[2] * 16135997692946310157UL) + ((uint64_t)op[3] * 2412456884085593928UL) + ((uint64_t)op[4] * 10102378073999687720UL) + ((uint64_t)op[5] * 17740737518988789911UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 75406472221L) - ((((int128)tmp_q[1] * 42934987440L) - ((int128)tmp_q[2] * 12413541574L) + ((int128)tmp_q[3] * 14919175951L) + ((int128)tmp_q[4] * 12143799504L) + ((int128)tmp_q[5] * 113515028576L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 113515028576L) - ((int128)tmp_q[1] * 75406472221L) - ((((int128)tmp_q[2] * 42934987440L) - ((int128)tmp_q[3] * 12413541574L) + ((int128)tmp_q[4] * 14919175951L) + ((int128)tmp_q[5] * 12143799504L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 12143799504L) + ((int128)tmp_q[1] * 113515028576L) - ((int128)tmp_q[2] * 75406472221L) - ((((int128)tmp_q[3] * 42934987440L) - ((int128)tmp_q[4] * 12413541574L) + ((int128)tmp_q[5] * 14919175951L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 14919175951L) + ((int128)tmp_q[1] * 12143799504L) + ((int128)tmp_q[2] * 113515028576L) - ((int128)tmp_q[3] * 75406472221L) - ((((int128)tmp_q[4] * 42934987440L) - ((int128)tmp_q[5] * 12413541574L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 12413541574L) + ((int128)tmp_q[1] * 14919175951L) + ((int128)tmp_q[2] * 12143799504L) + ((int128)tmp_q[3] * 113515028576L) - ((int128)tmp_q[4] * 75406472221L) - ((int128)tmp_q[5] * 85869974880L);
	tmp_zero[5] = ((int128)tmp_q[0] * 42934987440L) - ((int128)tmp_q[1] * 12413541574L) + ((int128)tmp_q[2] * 14919175951L) + ((int128)tmp_q[3] * 12143799504L) + ((int128)tmp_q[4] * 113515028576L) - ((int128)tmp_q[5] * 75406472221L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

