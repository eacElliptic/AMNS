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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13485411751727628859UL) + ((((uint64_t)op[1] * 12327683149522424072UL) + ((uint64_t)op[2] * 5266965531185242337UL) + ((uint64_t)op[3] * 18266668688124882787UL) + ((uint64_t)op[4] * 1768565902401614667UL) + ((uint64_t)op[5] * 2152917719180212607UL) + ((uint64_t)op[6] * 1944629471025519892UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 1944629471025519892UL) + ((uint64_t)op[1] * 13485411751727628859UL) + ((((uint64_t)op[2] * 12327683149522424072UL) + ((uint64_t)op[3] * 5266965531185242337UL) + ((uint64_t)op[4] * 18266668688124882787UL) + ((uint64_t)op[5] * 1768565902401614667UL) + ((uint64_t)op[6] * 2152917719180212607UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 2152917719180212607UL) + ((uint64_t)op[1] * 1944629471025519892UL) + ((uint64_t)op[2] * 13485411751727628859UL) + ((((uint64_t)op[3] * 12327683149522424072UL) + ((uint64_t)op[4] * 5266965531185242337UL) + ((uint64_t)op[5] * 18266668688124882787UL) + ((uint64_t)op[6] * 1768565902401614667UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 1768565902401614667UL) + ((uint64_t)op[1] * 2152917719180212607UL) + ((uint64_t)op[2] * 1944629471025519892UL) + ((uint64_t)op[3] * 13485411751727628859UL) + ((((uint64_t)op[4] * 12327683149522424072UL) + ((uint64_t)op[5] * 5266965531185242337UL) + ((uint64_t)op[6] * 18266668688124882787UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 18266668688124882787UL) + ((uint64_t)op[1] * 1768565902401614667UL) + ((uint64_t)op[2] * 2152917719180212607UL) + ((uint64_t)op[3] * 1944629471025519892UL) + ((uint64_t)op[4] * 13485411751727628859UL) + ((((uint64_t)op[5] * 12327683149522424072UL) + ((uint64_t)op[6] * 5266965531185242337UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 5266965531185242337UL) + ((uint64_t)op[1] * 18266668688124882787UL) + ((uint64_t)op[2] * 1768565902401614667UL) + ((uint64_t)op[3] * 2152917719180212607UL) + ((uint64_t)op[4] * 1944629471025519892UL) + ((uint64_t)op[5] * 13485411751727628859UL) + ((uint64_t)op[6] * 89561301148168984UL);
	tmp_q[6] = ((uint64_t)op[0] * 12327683149522424072UL) + ((uint64_t)op[1] * 5266965531185242337UL) + ((uint64_t)op[2] * 18266668688124882787UL) + ((uint64_t)op[3] * 1768565902401614667UL) + ((uint64_t)op[4] * 2152917719180212607UL) + ((uint64_t)op[5] * 1944629471025519892UL) + ((uint64_t)op[6] * 13485411751727628859UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 22974718903L) + ((((int128)tmp_q[1] * 35728307406L) + ((int128)tmp_q[2] * 6139862054L) - ((int128)tmp_q[3] * 5414528123L) + ((int128)tmp_q[4] * 55544032171L) - ((int128)tmp_q[5] * 34789274932L) - ((int128)tmp_q[6] * 52250425442L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 52250425442L) + ((int128)tmp_q[1] * 22974718903L) + ((((int128)tmp_q[2] * 35728307406L) + ((int128)tmp_q[3] * 6139862054L) - ((int128)tmp_q[4] * 5414528123L) + ((int128)tmp_q[5] * 55544032171L) - ((int128)tmp_q[6] * 34789274932L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 34789274932L) - ((int128)tmp_q[1] * 52250425442L) + ((int128)tmp_q[2] * 22974718903L) + ((((int128)tmp_q[3] * 35728307406L) + ((int128)tmp_q[4] * 6139862054L) - ((int128)tmp_q[5] * 5414528123L) + ((int128)tmp_q[6] * 55544032171L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 55544032171L) - ((int128)tmp_q[1] * 34789274932L) - ((int128)tmp_q[2] * 52250425442L) + ((int128)tmp_q[3] * 22974718903L) + ((((int128)tmp_q[4] * 35728307406L) + ((int128)tmp_q[5] * 6139862054L) - ((int128)tmp_q[6] * 5414528123L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 5414528123L) + ((int128)tmp_q[1] * 55544032171L) - ((int128)tmp_q[2] * 34789274932L) - ((int128)tmp_q[3] * 52250425442L) + ((int128)tmp_q[4] * 22974718903L) + ((((int128)tmp_q[5] * 35728307406L) + ((int128)tmp_q[6] * 6139862054L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 6139862054L) - ((int128)tmp_q[1] * 5414528123L) + ((int128)tmp_q[2] * 55544032171L) - ((int128)tmp_q[3] * 34789274932L) - ((int128)tmp_q[4] * 52250425442L) + ((int128)tmp_q[5] * 22974718903L) + ((int128)tmp_q[6] * 107184922218L);
	tmp_zero[6] = ((int128)tmp_q[0] * 35728307406L) + ((int128)tmp_q[1] * 6139862054L) - ((int128)tmp_q[2] * 5414528123L) + ((int128)tmp_q[3] * 55544032171L) - ((int128)tmp_q[4] * 34789274932L) - ((int128)tmp_q[5] * 52250425442L) + ((int128)tmp_q[6] * 22974718903L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

