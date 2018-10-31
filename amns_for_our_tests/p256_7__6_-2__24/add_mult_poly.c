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
	tmp_q[0] = ((uint64_t)op[0] * 3910625841485880599UL) + ((((uint64_t)op[1] * 8866670322834024305UL) + ((uint64_t)op[2] * 3803786514392401876UL) + ((uint64_t)op[3] * 4610794903108971736UL) + ((uint64_t)op[4] * 2842287504724998510UL) + ((uint64_t)op[5] * 11869179699193161128UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 11869179699193161128UL) + ((uint64_t)op[1] * 3910625841485880599UL) + ((((uint64_t)op[2] * 8866670322834024305UL) + ((uint64_t)op[3] * 3803786514392401876UL) + ((uint64_t)op[4] * 4610794903108971736UL) + ((uint64_t)op[5] * 2842287504724998510UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2842287504724998510UL) + ((uint64_t)op[1] * 11869179699193161128UL) + ((uint64_t)op[2] * 3910625841485880599UL) + ((((uint64_t)op[3] * 8866670322834024305UL) + ((uint64_t)op[4] * 3803786514392401876UL) + ((uint64_t)op[5] * 4610794903108971736UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 4610794903108971736UL) + ((uint64_t)op[1] * 2842287504724998510UL) + ((uint64_t)op[2] * 11869179699193161128UL) + ((uint64_t)op[3] * 3910625841485880599UL) + ((((uint64_t)op[4] * 8866670322834024305UL) + ((uint64_t)op[5] * 3803786514392401876UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 3803786514392401876UL) + ((uint64_t)op[1] * 4610794903108971736UL) + ((uint64_t)op[2] * 2842287504724998510UL) + ((uint64_t)op[3] * 11869179699193161128UL) + ((uint64_t)op[4] * 3910625841485880599UL) + ((uint64_t)op[5] * 713403428041503006UL);
	tmp_q[5] = ((uint64_t)op[0] * 8866670322834024305UL) + ((uint64_t)op[1] * 3803786514392401876UL) + ((uint64_t)op[2] * 4610794903108971736UL) + ((uint64_t)op[3] * 2842287504724998510UL) + ((uint64_t)op[4] * 11869179699193161128UL) + ((uint64_t)op[5] * 3910625841485880599UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3953161647601L) - ((((int128)tmp_q[1] * 2129607044881L) - ((int128)tmp_q[2] * 1835111005002L) - ((int128)tmp_q[3] * 2822460356076L) + ((int128)tmp_q[4] * 2583755172822L) - ((int128)tmp_q[5] * 3609229947808L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3609229947808L) + ((int128)tmp_q[1] * 3953161647601L) - ((((int128)tmp_q[2] * 2129607044881L) - ((int128)tmp_q[3] * 1835111005002L) - ((int128)tmp_q[4] * 2822460356076L) + ((int128)tmp_q[5] * 2583755172822L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 2583755172822L) - ((int128)tmp_q[1] * 3609229947808L) + ((int128)tmp_q[2] * 3953161647601L) - ((((int128)tmp_q[3] * 2129607044881L) - ((int128)tmp_q[4] * 1835111005002L) - ((int128)tmp_q[5] * 2822460356076L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2822460356076L) + ((int128)tmp_q[1] * 2583755172822L) - ((int128)tmp_q[2] * 3609229947808L) + ((int128)tmp_q[3] * 3953161647601L) - ((((int128)tmp_q[4] * 2129607044881L) - ((int128)tmp_q[5] * 1835111005002L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1835111005002L) - ((int128)tmp_q[1] * 2822460356076L) + ((int128)tmp_q[2] * 2583755172822L) - ((int128)tmp_q[3] * 3609229947808L) + ((int128)tmp_q[4] * 3953161647601L) - ((int128)tmp_q[5] * 4259214089762L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2129607044881L) - ((int128)tmp_q[1] * 1835111005002L) - ((int128)tmp_q[2] * 2822460356076L) + ((int128)tmp_q[3] * 2583755172822L) - ((int128)tmp_q[4] * 3609229947808L) + ((int128)tmp_q[5] * 3953161647601L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

