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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1928531154395413065UL) + ((((uint64_t)op[1] * 11452728137014091956UL) + ((uint64_t)op[2] * 8127050212441235016UL) + ((uint64_t)op[3] * 13405698020349902087UL) + ((uint64_t)op[4] * 17424890848275528133UL) + ((uint64_t)op[5] * 6449371278222096354UL) + ((uint64_t)op[6] * 1557958617802259182UL) + ((uint64_t)op[7] * 4894510289828809588UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 4894510289828809588UL) + ((uint64_t)op[1] * 1928531154395413065UL) + ((((uint64_t)op[2] * 11452728137014091956UL) + ((uint64_t)op[3] * 8127050212441235016UL) + ((uint64_t)op[4] * 13405698020349902087UL) + ((uint64_t)op[5] * 17424890848275528133UL) + ((uint64_t)op[6] * 6449371278222096354UL) + ((uint64_t)op[7] * 1557958617802259182UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 1557958617802259182UL) + ((uint64_t)op[1] * 4894510289828809588UL) + ((uint64_t)op[2] * 1928531154395413065UL) + ((((uint64_t)op[3] * 11452728137014091956UL) + ((uint64_t)op[4] * 8127050212441235016UL) + ((uint64_t)op[5] * 13405698020349902087UL) + ((uint64_t)op[6] * 17424890848275528133UL) + ((uint64_t)op[7] * 6449371278222096354UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 6449371278222096354UL) + ((uint64_t)op[1] * 1557958617802259182UL) + ((uint64_t)op[2] * 4894510289828809588UL) + ((uint64_t)op[3] * 1928531154395413065UL) + ((((uint64_t)op[4] * 11452728137014091956UL) + ((uint64_t)op[5] * 8127050212441235016UL) + ((uint64_t)op[6] * 13405698020349902087UL) + ((uint64_t)op[7] * 17424890848275528133UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 17424890848275528133UL) + ((uint64_t)op[1] * 6449371278222096354UL) + ((uint64_t)op[2] * 1557958617802259182UL) + ((uint64_t)op[3] * 4894510289828809588UL) + ((uint64_t)op[4] * 1928531154395413065UL) + ((((uint64_t)op[5] * 11452728137014091956UL) + ((uint64_t)op[6] * 8127050212441235016UL) + ((uint64_t)op[7] * 13405698020349902087UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 13405698020349902087UL) + ((uint64_t)op[1] * 17424890848275528133UL) + ((uint64_t)op[2] * 6449371278222096354UL) + ((uint64_t)op[3] * 1557958617802259182UL) + ((uint64_t)op[4] * 4894510289828809588UL) + ((uint64_t)op[5] * 1928531154395413065UL) + ((((uint64_t)op[6] * 11452728137014091956UL) + ((uint64_t)op[7] * 8127050212441235016UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 8127050212441235016UL) + ((uint64_t)op[1] * 13405698020349902087UL) + ((uint64_t)op[2] * 17424890848275528133UL) + ((uint64_t)op[3] * 6449371278222096354UL) + ((uint64_t)op[4] * 1557958617802259182UL) + ((uint64_t)op[5] * 4894510289828809588UL) + ((uint64_t)op[6] * 1928531154395413065UL) + ((uint64_t)op[7] * 6382120664260437228UL);
	tmp_q[7] = ((uint64_t)op[0] * 11452728137014091956UL) + ((uint64_t)op[1] * 8127050212441235016UL) + ((uint64_t)op[2] * 13405698020349902087UL) + ((uint64_t)op[3] * 17424890848275528133UL) + ((uint64_t)op[4] * 6449371278222096354UL) + ((uint64_t)op[5] * 1557958617802259182UL) + ((uint64_t)op[6] * 4894510289828809588UL) + ((uint64_t)op[7] * 1928531154395413065UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 886231200172L) + ((-((int128)tmp_q[1] * 83858902226536L) - ((int128)tmp_q[2] * 77280455277279L) - ((int128)tmp_q[3] * 9581712622158L) + ((int128)tmp_q[4] * 12180954715130L) - ((int128)tmp_q[5] * 124300498625825L) - ((int128)tmp_q[6] * 13810132935217L) + ((int128)tmp_q[7] * 158598614436118L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 158598614436118L) - ((int128)tmp_q[1] * 886231200172L) + ((-((int128)tmp_q[2] * 83858902226536L) - ((int128)tmp_q[3] * 77280455277279L) - ((int128)tmp_q[4] * 9581712622158L) + ((int128)tmp_q[5] * 12180954715130L) - ((int128)tmp_q[6] * 124300498625825L) - ((int128)tmp_q[7] * 13810132935217L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 13810132935217L) + ((int128)tmp_q[1] * 158598614436118L) - ((int128)tmp_q[2] * 886231200172L) + ((-((int128)tmp_q[3] * 83858902226536L) - ((int128)tmp_q[4] * 77280455277279L) - ((int128)tmp_q[5] * 9581712622158L) + ((int128)tmp_q[6] * 12180954715130L) - ((int128)tmp_q[7] * 124300498625825L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 124300498625825L) - ((int128)tmp_q[1] * 13810132935217L) + ((int128)tmp_q[2] * 158598614436118L) - ((int128)tmp_q[3] * 886231200172L) + ((-((int128)tmp_q[4] * 83858902226536L) - ((int128)tmp_q[5] * 77280455277279L) - ((int128)tmp_q[6] * 9581712622158L) + ((int128)tmp_q[7] * 12180954715130L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 12180954715130L) - ((int128)tmp_q[1] * 124300498625825L) - ((int128)tmp_q[2] * 13810132935217L) + ((int128)tmp_q[3] * 158598614436118L) - ((int128)tmp_q[4] * 886231200172L) + ((-((int128)tmp_q[5] * 83858902226536L) - ((int128)tmp_q[6] * 77280455277279L) - ((int128)tmp_q[7] * 9581712622158L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 9581712622158L) + ((int128)tmp_q[1] * 12180954715130L) - ((int128)tmp_q[2] * 124300498625825L) - ((int128)tmp_q[3] * 13810132935217L) + ((int128)tmp_q[4] * 158598614436118L) - ((int128)tmp_q[5] * 886231200172L) + ((-((int128)tmp_q[6] * 83858902226536L) - ((int128)tmp_q[7] * 77280455277279L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 77280455277279L) - ((int128)tmp_q[1] * 9581712622158L) + ((int128)tmp_q[2] * 12180954715130L) - ((int128)tmp_q[3] * 124300498625825L) - ((int128)tmp_q[4] * 13810132935217L) + ((int128)tmp_q[5] * 158598614436118L) - ((int128)tmp_q[6] * 886231200172L) - ((int128)tmp_q[7] * 587012315585752L);
	tmp_zero[7] = -((int128)tmp_q[0] * 83858902226536L) - ((int128)tmp_q[1] * 77280455277279L) - ((int128)tmp_q[2] * 9581712622158L) + ((int128)tmp_q[3] * 12180954715130L) - ((int128)tmp_q[4] * 124300498625825L) - ((int128)tmp_q[5] * 13810132935217L) + ((int128)tmp_q[6] * 158598614436118L) - ((int128)tmp_q[7] * 886231200172L);

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

