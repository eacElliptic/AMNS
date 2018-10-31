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
	tmp_q[0] = ((uint64_t)op[0] * 3533462319222442014UL) + ((((uint64_t)op[1] * 15665710081600427048UL) + ((uint64_t)op[2] * 14398242937726658089UL) + ((uint64_t)op[3] * 2652916518172900234UL) + ((uint64_t)op[4] * 4493482622192519192UL) + ((uint64_t)op[5] * 263968011730393540UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 263968011730393540UL) + ((uint64_t)op[1] * 3533462319222442014UL) + ((((uint64_t)op[2] * 15665710081600427048UL) + ((uint64_t)op[3] * 14398242937726658089UL) + ((uint64_t)op[4] * 2652916518172900234UL) + ((uint64_t)op[5] * 4493482622192519192UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 4493482622192519192UL) + ((uint64_t)op[1] * 263968011730393540UL) + ((uint64_t)op[2] * 3533462319222442014UL) + ((((uint64_t)op[3] * 15665710081600427048UL) + ((uint64_t)op[4] * 14398242937726658089UL) + ((uint64_t)op[5] * 2652916518172900234UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 2652916518172900234UL) + ((uint64_t)op[1] * 4493482622192519192UL) + ((uint64_t)op[2] * 263968011730393540UL) + ((uint64_t)op[3] * 3533462319222442014UL) + ((((uint64_t)op[4] * 15665710081600427048UL) + ((uint64_t)op[5] * 14398242937726658089UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 14398242937726658089UL) + ((uint64_t)op[1] * 2652916518172900234UL) + ((uint64_t)op[2] * 4493482622192519192UL) + ((uint64_t)op[3] * 263968011730393540UL) + ((uint64_t)op[4] * 3533462319222442014UL) + ((uint64_t)op[5] * 13905169960545622840UL);
	tmp_q[5] = ((uint64_t)op[0] * 15665710081600427048UL) + ((uint64_t)op[1] * 14398242937726658089UL) + ((uint64_t)op[2] * 2652916518172900234UL) + ((uint64_t)op[3] * 4493482622192519192UL) + ((uint64_t)op[4] * 263968011730393540UL) + ((uint64_t)op[5] * 3533462319222442014UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4720536856L) - ((((int128)tmp_q[1] * 19033664500L) + ((int128)tmp_q[2] * 85047552638L) + ((int128)tmp_q[3] * 99547132960L) - ((int128)tmp_q[4] * 35012715387L) + ((int128)tmp_q[5] * 113244579150L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 113244579150L) + ((int128)tmp_q[1] * 4720536856L) - ((((int128)tmp_q[2] * 19033664500L) + ((int128)tmp_q[3] * 85047552638L) + ((int128)tmp_q[4] * 99547132960L) - ((int128)tmp_q[5] * 35012715387L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 35012715387L) + ((int128)tmp_q[1] * 113244579150L) + ((int128)tmp_q[2] * 4720536856L) - ((((int128)tmp_q[3] * 19033664500L) + ((int128)tmp_q[4] * 85047552638L) + ((int128)tmp_q[5] * 99547132960L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 99547132960L) - ((int128)tmp_q[1] * 35012715387L) + ((int128)tmp_q[2] * 113244579150L) + ((int128)tmp_q[3] * 4720536856L) - ((((int128)tmp_q[4] * 19033664500L) + ((int128)tmp_q[5] * 85047552638L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 85047552638L) + ((int128)tmp_q[1] * 99547132960L) - ((int128)tmp_q[2] * 35012715387L) + ((int128)tmp_q[3] * 113244579150L) + ((int128)tmp_q[4] * 4720536856L) - ((int128)tmp_q[5] * 95168322500L);
	tmp_zero[5] = ((int128)tmp_q[0] * 19033664500L) + ((int128)tmp_q[1] * 85047552638L) + ((int128)tmp_q[2] * 99547132960L) - ((int128)tmp_q[3] * 35012715387L) + ((int128)tmp_q[4] * 113244579150L) + ((int128)tmp_q[5] * 4720536856L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

