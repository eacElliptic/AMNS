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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10784050099739641021UL) + ((((uint64_t)op[1] * 2255633894185717448UL) + ((uint64_t)op[2] * 8735007387696146306UL) + ((uint64_t)op[3] * 1146659682499932707UL) + ((uint64_t)op[4] * 9212879237597363052UL) + ((uint64_t)op[5] * 225542354071717181UL) + ((uint64_t)op[6] * 13199206665590782225UL) + ((uint64_t)op[7] * 16786944494442192043UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 16786944494442192043UL) + ((uint64_t)op[1] * 10784050099739641021UL) + ((((uint64_t)op[2] * 2255633894185717448UL) + ((uint64_t)op[3] * 8735007387696146306UL) + ((uint64_t)op[4] * 1146659682499932707UL) + ((uint64_t)op[5] * 9212879237597363052UL) + ((uint64_t)op[6] * 225542354071717181UL) + ((uint64_t)op[7] * 13199206665590782225UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 13199206665590782225UL) + ((uint64_t)op[1] * 16786944494442192043UL) + ((uint64_t)op[2] * 10784050099739641021UL) + ((((uint64_t)op[3] * 2255633894185717448UL) + ((uint64_t)op[4] * 8735007387696146306UL) + ((uint64_t)op[5] * 1146659682499932707UL) + ((uint64_t)op[6] * 9212879237597363052UL) + ((uint64_t)op[7] * 225542354071717181UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 225542354071717181UL) + ((uint64_t)op[1] * 13199206665590782225UL) + ((uint64_t)op[2] * 16786944494442192043UL) + ((uint64_t)op[3] * 10784050099739641021UL) + ((((uint64_t)op[4] * 2255633894185717448UL) + ((uint64_t)op[5] * 8735007387696146306UL) + ((uint64_t)op[6] * 1146659682499932707UL) + ((uint64_t)op[7] * 9212879237597363052UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 9212879237597363052UL) + ((uint64_t)op[1] * 225542354071717181UL) + ((uint64_t)op[2] * 13199206665590782225UL) + ((uint64_t)op[3] * 16786944494442192043UL) + ((uint64_t)op[4] * 10784050099739641021UL) + ((((uint64_t)op[5] * 2255633894185717448UL) + ((uint64_t)op[6] * 8735007387696146306UL) + ((uint64_t)op[7] * 1146659682499932707UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 1146659682499932707UL) + ((uint64_t)op[1] * 9212879237597363052UL) + ((uint64_t)op[2] * 225542354071717181UL) + ((uint64_t)op[3] * 13199206665590782225UL) + ((uint64_t)op[4] * 16786944494442192043UL) + ((uint64_t)op[5] * 10784050099739641021UL) + ((((uint64_t)op[6] * 2255633894185717448UL) + ((uint64_t)op[7] * 8735007387696146306UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 8735007387696146306UL) + ((uint64_t)op[1] * 1146659682499932707UL) + ((uint64_t)op[2] * 9212879237597363052UL) + ((uint64_t)op[3] * 225542354071717181UL) + ((uint64_t)op[4] * 13199206665590782225UL) + ((uint64_t)op[5] * 16786944494442192043UL) + ((uint64_t)op[6] * 10784050099739641021UL) + ((uint64_t)op[7] * 7168574602780964376UL);
	tmp_q[7] = ((uint64_t)op[0] * 2255633894185717448UL) + ((uint64_t)op[1] * 8735007387696146306UL) + ((uint64_t)op[2] * 1146659682499932707UL) + ((uint64_t)op[3] * 9212879237597363052UL) + ((uint64_t)op[4] * 225542354071717181UL) + ((uint64_t)op[5] * 13199206665590782225UL) + ((uint64_t)op[6] * 16786944494442192043UL) + ((uint64_t)op[7] * 10784050099739641021UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 165904740762209L) - ((((int128)tmp_q[1] * 62111579610476L) - ((int128)tmp_q[2] * 151622209691075L) - ((int128)tmp_q[3] * 6493661565959L) - ((int128)tmp_q[4] * 25652912948860L) - ((int128)tmp_q[5] * 61639842363622L) + ((int128)tmp_q[6] * 25469012134120L) + ((int128)tmp_q[7] * 137154458786384L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 137154458786384L) - ((int128)tmp_q[1] * 165904740762209L) - ((((int128)tmp_q[2] * 62111579610476L) - ((int128)tmp_q[3] * 151622209691075L) - ((int128)tmp_q[4] * 6493661565959L) - ((int128)tmp_q[5] * 25652912948860L) - ((int128)tmp_q[6] * 61639842363622L) + ((int128)tmp_q[7] * 25469012134120L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 25469012134120L) + ((int128)tmp_q[1] * 137154458786384L) - ((int128)tmp_q[2] * 165904740762209L) - ((((int128)tmp_q[3] * 62111579610476L) - ((int128)tmp_q[4] * 151622209691075L) - ((int128)tmp_q[5] * 6493661565959L) - ((int128)tmp_q[6] * 25652912948860L) - ((int128)tmp_q[7] * 61639842363622L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 61639842363622L) + ((int128)tmp_q[1] * 25469012134120L) + ((int128)tmp_q[2] * 137154458786384L) - ((int128)tmp_q[3] * 165904740762209L) - ((((int128)tmp_q[4] * 62111579610476L) - ((int128)tmp_q[5] * 151622209691075L) - ((int128)tmp_q[6] * 6493661565959L) - ((int128)tmp_q[7] * 25652912948860L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 25652912948860L) - ((int128)tmp_q[1] * 61639842363622L) + ((int128)tmp_q[2] * 25469012134120L) + ((int128)tmp_q[3] * 137154458786384L) - ((int128)tmp_q[4] * 165904740762209L) - ((((int128)tmp_q[5] * 62111579610476L) - ((int128)tmp_q[6] * 151622209691075L) - ((int128)tmp_q[7] * 6493661565959L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 6493661565959L) - ((int128)tmp_q[1] * 25652912948860L) - ((int128)tmp_q[2] * 61639842363622L) + ((int128)tmp_q[3] * 25469012134120L) + ((int128)tmp_q[4] * 137154458786384L) - ((int128)tmp_q[5] * 165904740762209L) - ((((int128)tmp_q[6] * 62111579610476L) - ((int128)tmp_q[7] * 151622209691075L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 151622209691075L) - ((int128)tmp_q[1] * 6493661565959L) - ((int128)tmp_q[2] * 25652912948860L) - ((int128)tmp_q[3] * 61639842363622L) + ((int128)tmp_q[4] * 25469012134120L) + ((int128)tmp_q[5] * 137154458786384L) - ((int128)tmp_q[6] * 165904740762209L) - ((int128)tmp_q[7] * 310557898052380L);
	tmp_zero[7] = ((int128)tmp_q[0] * 62111579610476L) - ((int128)tmp_q[1] * 151622209691075L) - ((int128)tmp_q[2] * 6493661565959L) - ((int128)tmp_q[3] * 25652912948860L) - ((int128)tmp_q[4] * 61639842363622L) + ((int128)tmp_q[5] * 25469012134120L) + ((int128)tmp_q[6] * 137154458786384L) - ((int128)tmp_q[7] * 165904740762209L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

