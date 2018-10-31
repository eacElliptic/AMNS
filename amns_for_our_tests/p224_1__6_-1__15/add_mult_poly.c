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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16501328160651383777UL) + ((((uint64_t)op[1] * 12827555012385284180UL) + ((uint64_t)op[2] * 1985945785608092665UL) + ((uint64_t)op[3] * 14702767068665144442UL) + ((uint64_t)op[4] * 14146863575974876062UL) + ((uint64_t)op[5] * 2960314648851010357UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 2960314648851010357UL) + ((uint64_t)op[1] * 16501328160651383777UL) + ((((uint64_t)op[2] * 12827555012385284180UL) + ((uint64_t)op[3] * 1985945785608092665UL) + ((uint64_t)op[4] * 14702767068665144442UL) + ((uint64_t)op[5] * 14146863575974876062UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 14146863575974876062UL) + ((uint64_t)op[1] * 2960314648851010357UL) + ((uint64_t)op[2] * 16501328160651383777UL) + ((((uint64_t)op[3] * 12827555012385284180UL) + ((uint64_t)op[4] * 1985945785608092665UL) + ((uint64_t)op[5] * 14702767068665144442UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 14702767068665144442UL) + ((uint64_t)op[1] * 14146863575974876062UL) + ((uint64_t)op[2] * 2960314648851010357UL) + ((uint64_t)op[3] * 16501328160651383777UL) + ((((uint64_t)op[4] * 12827555012385284180UL) + ((uint64_t)op[5] * 1985945785608092665UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 1985945785608092665UL) + ((uint64_t)op[1] * 14702767068665144442UL) + ((uint64_t)op[2] * 14146863575974876062UL) + ((uint64_t)op[3] * 2960314648851010357UL) + ((uint64_t)op[4] * 16501328160651383777UL) + ((uint64_t)op[5] * 5619189061324267436UL);
	tmp_q[5] = ((uint64_t)op[0] * 12827555012385284180UL) + ((uint64_t)op[1] * 1985945785608092665UL) + ((uint64_t)op[2] * 14702767068665144442UL) + ((uint64_t)op[3] * 14146863575974876062UL) + ((uint64_t)op[4] * 2960314648851010357UL) + ((uint64_t)op[5] * 16501328160651383777UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1201296798897013L) - (-((int128)tmp_q[1] * 22775465338234822L) + ((int128)tmp_q[2] * 33322420335829497L) - ((int128)tmp_q[3] * 9911721338919124L) + ((int128)tmp_q[4] * 34523717134726506L) + ((int128)tmp_q[5] * 12863743999315697L));
	tmp_zero[1] = ((int128)tmp_q[0] * 12863743999315697L) + ((int128)tmp_q[1] * 1201296798897013L) - (-((int128)tmp_q[2] * 22775465338234822L) + ((int128)tmp_q[3] * 33322420335829497L) - ((int128)tmp_q[4] * 9911721338919124L) + ((int128)tmp_q[5] * 34523717134726506L));
	tmp_zero[2] = ((int128)tmp_q[0] * 34523717134726506L) + ((int128)tmp_q[1] * 12863743999315697L) + ((int128)tmp_q[2] * 1201296798897013L) - (-((int128)tmp_q[3] * 22775465338234822L) + ((int128)tmp_q[4] * 33322420335829497L) - ((int128)tmp_q[5] * 9911721338919124L));
	tmp_zero[3] = -((int128)tmp_q[0] * 9911721338919124L) + ((int128)tmp_q[1] * 34523717134726506L) + ((int128)tmp_q[2] * 12863743999315697L) + ((int128)tmp_q[3] * 1201296798897013L) - (-((int128)tmp_q[4] * 22775465338234822L) + ((int128)tmp_q[5] * 33322420335829497L));
	tmp_zero[4] = ((int128)tmp_q[0] * 33322420335829497L) - ((int128)tmp_q[1] * 9911721338919124L) + ((int128)tmp_q[2] * 34523717134726506L) + ((int128)tmp_q[3] * 12863743999315697L) + ((int128)tmp_q[4] * 1201296798897013L) + ((int128)tmp_q[5] * 22775465338234822L);
	tmp_zero[5] = -((int128)tmp_q[0] * 22775465338234822L) + ((int128)tmp_q[1] * 33322420335829497L) - ((int128)tmp_q[2] * 9911721338919124L) + ((int128)tmp_q[3] * 34523717134726506L) + ((int128)tmp_q[4] * 12863743999315697L) + ((int128)tmp_q[5] * 1201296798897013L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

