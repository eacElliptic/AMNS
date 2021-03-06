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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10404632948360274676UL) + ((((uint64_t)op[1] * 9295201629628824626UL) + ((uint64_t)op[2] * 216028180068067394UL) + ((uint64_t)op[3] * 4806064102166829882UL) + ((uint64_t)op[4] * 11078033927222974377UL) + ((uint64_t)op[5] * 759983527995208239UL) + ((uint64_t)op[6] * 16058229768117021811UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16058229768117021811UL) + ((uint64_t)op[1] * 10404632948360274676UL) + ((((uint64_t)op[2] * 9295201629628824626UL) + ((uint64_t)op[3] * 216028180068067394UL) + ((uint64_t)op[4] * 4806064102166829882UL) + ((uint64_t)op[5] * 11078033927222974377UL) + ((uint64_t)op[6] * 759983527995208239UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 759983527995208239UL) + ((uint64_t)op[1] * 16058229768117021811UL) + ((uint64_t)op[2] * 10404632948360274676UL) + ((((uint64_t)op[3] * 9295201629628824626UL) + ((uint64_t)op[4] * 216028180068067394UL) + ((uint64_t)op[5] * 4806064102166829882UL) + ((uint64_t)op[6] * 11078033927222974377UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 11078033927222974377UL) + ((uint64_t)op[1] * 759983527995208239UL) + ((uint64_t)op[2] * 16058229768117021811UL) + ((uint64_t)op[3] * 10404632948360274676UL) + ((((uint64_t)op[4] * 9295201629628824626UL) + ((uint64_t)op[5] * 216028180068067394UL) + ((uint64_t)op[6] * 4806064102166829882UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 4806064102166829882UL) + ((uint64_t)op[1] * 11078033927222974377UL) + ((uint64_t)op[2] * 759983527995208239UL) + ((uint64_t)op[3] * 16058229768117021811UL) + ((uint64_t)op[4] * 10404632948360274676UL) + ((((uint64_t)op[5] * 9295201629628824626UL) + ((uint64_t)op[6] * 216028180068067394UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 216028180068067394UL) + ((uint64_t)op[1] * 4806064102166829882UL) + ((uint64_t)op[2] * 11078033927222974377UL) + ((uint64_t)op[3] * 759983527995208239UL) + ((uint64_t)op[4] * 16058229768117021811UL) + ((uint64_t)op[5] * 10404632948360274676UL) + ((uint64_t)op[6] * 9007883258532629354UL);
	tmp_q[6] = ((uint64_t)op[0] * 9295201629628824626UL) + ((uint64_t)op[1] * 216028180068067394UL) + ((uint64_t)op[2] * 4806064102166829882UL) + ((uint64_t)op[3] * 11078033927222974377UL) + ((uint64_t)op[4] * 759983527995208239UL) + ((uint64_t)op[5] * 16058229768117021811UL) + ((uint64_t)op[6] * 10404632948360274676UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 52669427398L) - ((((int128)tmp_q[1] * 10197253589L) - ((int128)tmp_q[2] * 20271055481L) + ((int128)tmp_q[3] * 58331311239L) + ((int128)tmp_q[4] * 9832118170L) + ((int128)tmp_q[5] * 2476005829L) + ((int128)tmp_q[6] * 47308276887L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 47308276887L) - ((int128)tmp_q[1] * 52669427398L) - ((((int128)tmp_q[2] * 10197253589L) - ((int128)tmp_q[3] * 20271055481L) + ((int128)tmp_q[4] * 58331311239L) + ((int128)tmp_q[5] * 9832118170L) + ((int128)tmp_q[6] * 2476005829L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2476005829L) + ((int128)tmp_q[1] * 47308276887L) - ((int128)tmp_q[2] * 52669427398L) - ((((int128)tmp_q[3] * 10197253589L) - ((int128)tmp_q[4] * 20271055481L) + ((int128)tmp_q[5] * 58331311239L) + ((int128)tmp_q[6] * 9832118170L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 9832118170L) + ((int128)tmp_q[1] * 2476005829L) + ((int128)tmp_q[2] * 47308276887L) - ((int128)tmp_q[3] * 52669427398L) - ((((int128)tmp_q[4] * 10197253589L) - ((int128)tmp_q[5] * 20271055481L) + ((int128)tmp_q[6] * 58331311239L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 58331311239L) + ((int128)tmp_q[1] * 9832118170L) + ((int128)tmp_q[2] * 2476005829L) + ((int128)tmp_q[3] * 47308276887L) - ((int128)tmp_q[4] * 52669427398L) - ((((int128)tmp_q[5] * 10197253589L) - ((int128)tmp_q[6] * 20271055481L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 20271055481L) + ((int128)tmp_q[1] * 58331311239L) + ((int128)tmp_q[2] * 9832118170L) + ((int128)tmp_q[3] * 2476005829L) + ((int128)tmp_q[4] * 47308276887L) - ((int128)tmp_q[5] * 52669427398L) - ((int128)tmp_q[6] * 30591760767L);
	tmp_zero[6] = ((int128)tmp_q[0] * 10197253589L) - ((int128)tmp_q[1] * 20271055481L) + ((int128)tmp_q[2] * 58331311239L) + ((int128)tmp_q[3] * 9832118170L) + ((int128)tmp_q[4] * 2476005829L) + ((int128)tmp_q[5] * 47308276887L) - ((int128)tmp_q[6] * 52669427398L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

