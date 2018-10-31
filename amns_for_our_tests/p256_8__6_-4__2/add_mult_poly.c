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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17063656555731682565UL) + ((((uint64_t)op[1] * 9246393160275671650UL) + ((uint64_t)op[2] * 6104390915319896689UL) + ((uint64_t)op[3] * 13122629497368839886UL) + ((uint64_t)op[4] * 7718877475540000672UL) + ((uint64_t)op[5] * 5688508637058477063UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 5688508637058477063UL) + ((uint64_t)op[1] * 17063656555731682565UL) + ((((uint64_t)op[2] * 9246393160275671650UL) + ((uint64_t)op[3] * 6104390915319896689UL) + ((uint64_t)op[4] * 13122629497368839886UL) + ((uint64_t)op[5] * 7718877475540000672UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 7718877475540000672UL) + ((uint64_t)op[1] * 5688508637058477063UL) + ((uint64_t)op[2] * 17063656555731682565UL) + ((((uint64_t)op[3] * 9246393160275671650UL) + ((uint64_t)op[4] * 6104390915319896689UL) + ((uint64_t)op[5] * 13122629497368839886UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 13122629497368839886UL) + ((uint64_t)op[1] * 7718877475540000672UL) + ((uint64_t)op[2] * 5688508637058477063UL) + ((uint64_t)op[3] * 17063656555731682565UL) + ((((uint64_t)op[4] * 9246393160275671650UL) + ((uint64_t)op[5] * 6104390915319896689UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 6104390915319896689UL) + ((uint64_t)op[1] * 13122629497368839886UL) + ((uint64_t)op[2] * 7718877475540000672UL) + ((uint64_t)op[3] * 5688508637058477063UL) + ((uint64_t)op[4] * 17063656555731682565UL) + ((uint64_t)op[5] * 18354659580025968248UL);
	tmp_q[5] = ((uint64_t)op[0] * 9246393160275671650UL) + ((uint64_t)op[1] * 6104390915319896689UL) + ((uint64_t)op[2] * 13122629497368839886UL) + ((uint64_t)op[3] * 7718877475540000672UL) + ((uint64_t)op[4] * 5688508637058477063UL) + ((uint64_t)op[5] * 17063656555731682565UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 300087228997L) - ((-((int128)tmp_q[1] * 2190727212583L) - ((int128)tmp_q[2] * 3299473112040L) + ((int128)tmp_q[3] * 4197192561645L) - ((int128)tmp_q[4] * 1402364855665L) - ((int128)tmp_q[5] * 2813242274301L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 2813242274301L) - ((int128)tmp_q[1] * 300087228997L) - ((-((int128)tmp_q[2] * 2190727212583L) - ((int128)tmp_q[3] * 3299473112040L) + ((int128)tmp_q[4] * 4197192561645L) - ((int128)tmp_q[5] * 1402364855665L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 1402364855665L) - ((int128)tmp_q[1] * 2813242274301L) - ((int128)tmp_q[2] * 300087228997L) - ((-((int128)tmp_q[3] * 2190727212583L) - ((int128)tmp_q[4] * 3299473112040L) + ((int128)tmp_q[5] * 4197192561645L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 4197192561645L) - ((int128)tmp_q[1] * 1402364855665L) - ((int128)tmp_q[2] * 2813242274301L) - ((int128)tmp_q[3] * 300087228997L) - ((-((int128)tmp_q[4] * 2190727212583L) - ((int128)tmp_q[5] * 3299473112040L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 3299473112040L) + ((int128)tmp_q[1] * 4197192561645L) - ((int128)tmp_q[2] * 1402364855665L) - ((int128)tmp_q[3] * 2813242274301L) - ((int128)tmp_q[4] * 300087228997L) + ((int128)tmp_q[5] * 8762908850332L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2190727212583L) - ((int128)tmp_q[1] * 3299473112040L) + ((int128)tmp_q[2] * 4197192561645L) - ((int128)tmp_q[3] * 1402364855665L) - ((int128)tmp_q[4] * 2813242274301L) - ((int128)tmp_q[5] * 300087228997L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

