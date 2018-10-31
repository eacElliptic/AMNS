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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16699456423906297011UL) + ((((uint64_t)op[1] * 16882037927047396494UL) + ((uint64_t)op[2] * 5064942502019342525UL) + ((uint64_t)op[3] * 5448271459658808049UL) + ((uint64_t)op[4] * 17167885498952964414UL) + ((uint64_t)op[5] * 4990998245263221594UL) + ((uint64_t)op[6] * 7890244193801106025UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 7890244193801106025UL) + ((uint64_t)op[1] * 16699456423906297011UL) + ((((uint64_t)op[2] * 16882037927047396494UL) + ((uint64_t)op[3] * 5064942502019342525UL) + ((uint64_t)op[4] * 5448271459658808049UL) + ((uint64_t)op[5] * 17167885498952964414UL) + ((uint64_t)op[6] * 4990998245263221594UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 4990998245263221594UL) + ((uint64_t)op[1] * 7890244193801106025UL) + ((uint64_t)op[2] * 16699456423906297011UL) + ((((uint64_t)op[3] * 16882037927047396494UL) + ((uint64_t)op[4] * 5064942502019342525UL) + ((uint64_t)op[5] * 5448271459658808049UL) + ((uint64_t)op[6] * 17167885498952964414UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 17167885498952964414UL) + ((uint64_t)op[1] * 4990998245263221594UL) + ((uint64_t)op[2] * 7890244193801106025UL) + ((uint64_t)op[3] * 16699456423906297011UL) + ((((uint64_t)op[4] * 16882037927047396494UL) + ((uint64_t)op[5] * 5064942502019342525UL) + ((uint64_t)op[6] * 5448271459658808049UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 5448271459658808049UL) + ((uint64_t)op[1] * 17167885498952964414UL) + ((uint64_t)op[2] * 4990998245263221594UL) + ((uint64_t)op[3] * 7890244193801106025UL) + ((uint64_t)op[4] * 16699456423906297011UL) + ((((uint64_t)op[5] * 16882037927047396494UL) + ((uint64_t)op[6] * 5064942502019342525UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 5064942502019342525UL) + ((uint64_t)op[1] * 5448271459658808049UL) + ((uint64_t)op[2] * 17167885498952964414UL) + ((uint64_t)op[3] * 4990998245263221594UL) + ((uint64_t)op[4] * 7890244193801106025UL) + ((uint64_t)op[5] * 16699456423906297011UL) + ((uint64_t)op[6] * 12517649173297240976UL);
	tmp_q[6] = ((uint64_t)op[0] * 16882037927047396494UL) + ((uint64_t)op[1] * 5064942502019342525UL) + ((uint64_t)op[2] * 5448271459658808049UL) + ((uint64_t)op[3] * 17167885498952964414UL) + ((uint64_t)op[4] * 4990998245263221594UL) + ((uint64_t)op[5] * 7890244193801106025UL) + ((uint64_t)op[6] * 16699456423906297011UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 25956347189L) - ((((int128)tmp_q[1] * 52693138586L) + ((int128)tmp_q[2] * 24422003878L) + ((int128)tmp_q[3] * 4696639620L) - ((int128)tmp_q[4] * 38208918477L) + ((int128)tmp_q[5] * 59640429567L) - ((int128)tmp_q[6] * 7595735319L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 7595735319L) + ((int128)tmp_q[1] * 25956347189L) - ((((int128)tmp_q[2] * 52693138586L) + ((int128)tmp_q[3] * 24422003878L) + ((int128)tmp_q[4] * 4696639620L) - ((int128)tmp_q[5] * 38208918477L) + ((int128)tmp_q[6] * 59640429567L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 59640429567L) - ((int128)tmp_q[1] * 7595735319L) + ((int128)tmp_q[2] * 25956347189L) - ((((int128)tmp_q[3] * 52693138586L) + ((int128)tmp_q[4] * 24422003878L) + ((int128)tmp_q[5] * 4696639620L) - ((int128)tmp_q[6] * 38208918477L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 38208918477L) + ((int128)tmp_q[1] * 59640429567L) - ((int128)tmp_q[2] * 7595735319L) + ((int128)tmp_q[3] * 25956347189L) - ((((int128)tmp_q[4] * 52693138586L) + ((int128)tmp_q[5] * 24422003878L) + ((int128)tmp_q[6] * 4696639620L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 4696639620L) - ((int128)tmp_q[1] * 38208918477L) + ((int128)tmp_q[2] * 59640429567L) - ((int128)tmp_q[3] * 7595735319L) + ((int128)tmp_q[4] * 25956347189L) - ((((int128)tmp_q[5] * 52693138586L) + ((int128)tmp_q[6] * 24422003878L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 24422003878L) + ((int128)tmp_q[1] * 4696639620L) - ((int128)tmp_q[2] * 38208918477L) + ((int128)tmp_q[3] * 59640429567L) - ((int128)tmp_q[4] * 7595735319L) + ((int128)tmp_q[5] * 25956347189L) - ((int128)tmp_q[6] * 421545108688L);
	tmp_zero[6] = ((int128)tmp_q[0] * 52693138586L) + ((int128)tmp_q[1] * 24422003878L) + ((int128)tmp_q[2] * 4696639620L) - ((int128)tmp_q[3] * 38208918477L) + ((int128)tmp_q[4] * 59640429567L) - ((int128)tmp_q[5] * 7595735319L) + ((int128)tmp_q[6] * 25956347189L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

