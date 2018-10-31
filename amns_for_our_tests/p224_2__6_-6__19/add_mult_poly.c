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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4817354972920368763UL) + ((((uint64_t)op[1] * 1259982272013502782UL) + ((uint64_t)op[2] * 2645640152882856394UL) + ((uint64_t)op[3] * 10803778160660387465UL) + ((uint64_t)op[4] * 15108230011247219638UL) + ((uint64_t)op[5] * 14642156543533888575UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 14642156543533888575UL) + ((uint64_t)op[1] * 4817354972920368763UL) + ((((uint64_t)op[2] * 1259982272013502782UL) + ((uint64_t)op[3] * 2645640152882856394UL) + ((uint64_t)op[4] * 10803778160660387465UL) + ((uint64_t)op[5] * 15108230011247219638UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 15108230011247219638UL) + ((uint64_t)op[1] * 14642156543533888575UL) + ((uint64_t)op[2] * 4817354972920368763UL) + ((((uint64_t)op[3] * 1259982272013502782UL) + ((uint64_t)op[4] * 2645640152882856394UL) + ((uint64_t)op[5] * 10803778160660387465UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 10803778160660387465UL) + ((uint64_t)op[1] * 15108230011247219638UL) + ((uint64_t)op[2] * 14642156543533888575UL) + ((uint64_t)op[3] * 4817354972920368763UL) + ((((uint64_t)op[4] * 1259982272013502782UL) + ((uint64_t)op[5] * 2645640152882856394UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 2645640152882856394UL) + ((uint64_t)op[1] * 10803778160660387465UL) + ((uint64_t)op[2] * 15108230011247219638UL) + ((uint64_t)op[3] * 14642156543533888575UL) + ((uint64_t)op[4] * 4817354972920368763UL) + ((uint64_t)op[5] * 10886850441628534924UL);
	tmp_q[5] = ((uint64_t)op[0] * 1259982272013502782UL) + ((uint64_t)op[1] * 2645640152882856394UL) + ((uint64_t)op[2] * 10803778160660387465UL) + ((uint64_t)op[3] * 15108230011247219638UL) + ((uint64_t)op[4] * 14642156543533888575UL) + ((uint64_t)op[5] * 4817354972920368763UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 58888950655L) - ((((int128)tmp_q[1] * 29428073534L) - ((int128)tmp_q[2] * 29378474413L) + ((int128)tmp_q[3] * 80899944670L) - ((int128)tmp_q[4] * 37452933795L) + ((int128)tmp_q[5] * 82365359365L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 82365359365L) - ((int128)tmp_q[1] * 58888950655L) - ((((int128)tmp_q[2] * 29428073534L) - ((int128)tmp_q[3] * 29378474413L) + ((int128)tmp_q[4] * 80899944670L) - ((int128)tmp_q[5] * 37452933795L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 37452933795L) + ((int128)tmp_q[1] * 82365359365L) - ((int128)tmp_q[2] * 58888950655L) - ((((int128)tmp_q[3] * 29428073534L) - ((int128)tmp_q[4] * 29378474413L) + ((int128)tmp_q[5] * 80899944670L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 80899944670L) - ((int128)tmp_q[1] * 37452933795L) + ((int128)tmp_q[2] * 82365359365L) - ((int128)tmp_q[3] * 58888950655L) - ((((int128)tmp_q[4] * 29428073534L) - ((int128)tmp_q[5] * 29378474413L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 29378474413L) + ((int128)tmp_q[1] * 80899944670L) - ((int128)tmp_q[2] * 37452933795L) + ((int128)tmp_q[3] * 82365359365L) - ((int128)tmp_q[4] * 58888950655L) - ((int128)tmp_q[5] * 176568441204L);
	tmp_zero[5] = ((int128)tmp_q[0] * 29428073534L) - ((int128)tmp_q[1] * 29378474413L) + ((int128)tmp_q[2] * 80899944670L) - ((int128)tmp_q[3] * 37452933795L) + ((int128)tmp_q[4] * 82365359365L) - ((int128)tmp_q[5] * 58888950655L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

