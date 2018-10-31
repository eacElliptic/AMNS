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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6543897660486689569UL) + ((((uint64_t)op[1] * 11354102283748031673UL) + ((uint64_t)op[2] * 11605163401546896190UL) + ((uint64_t)op[3] * 15885973952608620966UL) + ((uint64_t)op[4] * 814064101681862423UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 814064101681862423UL) + ((uint64_t)op[1] * 6543897660486689569UL) + ((((uint64_t)op[2] * 11354102283748031673UL) + ((uint64_t)op[3] * 11605163401546896190UL) + ((uint64_t)op[4] * 15885973952608620966UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 15885973952608620966UL) + ((uint64_t)op[1] * 814064101681862423UL) + ((uint64_t)op[2] * 6543897660486689569UL) + ((((uint64_t)op[3] * 11354102283748031673UL) + ((uint64_t)op[4] * 11605163401546896190UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 11605163401546896190UL) + ((uint64_t)op[1] * 15885973952608620966UL) + ((uint64_t)op[2] * 814064101681862423UL) + ((uint64_t)op[3] * 6543897660486689569UL) + ((uint64_t)op[4] * 14185283579923039886UL);
	tmp_q[4] = ((uint64_t)op[0] * 11354102283748031673UL) + ((uint64_t)op[1] * 11605163401546896190UL) + ((uint64_t)op[2] * 15885973952608620966UL) + ((uint64_t)op[3] * 814064101681862423UL) + ((uint64_t)op[4] * 6543897660486689569UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 54200865529L) - ((-((int128)tmp_q[1] * 69471286590L) - ((int128)tmp_q[2] * 57958752425L) + ((int128)tmp_q[3] * 252612188315L) - ((int128)tmp_q[4] * 113919251709L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 113919251709L) + ((int128)tmp_q[1] * 54200865529L) - ((-((int128)tmp_q[2] * 69471286590L) - ((int128)tmp_q[3] * 57958752425L) + ((int128)tmp_q[4] * 252612188315L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 252612188315L) - ((int128)tmp_q[1] * 113919251709L) + ((int128)tmp_q[2] * 54200865529L) - ((-((int128)tmp_q[3] * 69471286590L) - ((int128)tmp_q[4] * 57958752425L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 57958752425L) + ((int128)tmp_q[1] * 252612188315L) - ((int128)tmp_q[2] * 113919251709L) + ((int128)tmp_q[3] * 54200865529L) + ((int128)tmp_q[4] * 138942573180L);
	tmp_zero[4] = -((int128)tmp_q[0] * 69471286590L) - ((int128)tmp_q[1] * 57958752425L) + ((int128)tmp_q[2] * 252612188315L) - ((int128)tmp_q[3] * 113919251709L) + ((int128)tmp_q[4] * 54200865529L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

