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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10267894324072613377UL) + ((((uint64_t)op[1] * 4147259510870461806UL) + ((uint64_t)op[2] * 12177068667683431813UL) + ((uint64_t)op[3] * 2398662181628872264UL) + ((uint64_t)op[4] * 13174311866102213863UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 13174311866102213863UL) + ((uint64_t)op[1] * 10267894324072613377UL) + ((((uint64_t)op[2] * 4147259510870461806UL) + ((uint64_t)op[3] * 12177068667683431813UL) + ((uint64_t)op[4] * 2398662181628872264UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 2398662181628872264UL) + ((uint64_t)op[1] * 13174311866102213863UL) + ((uint64_t)op[2] * 10267894324072613377UL) + ((((uint64_t)op[3] * 4147259510870461806UL) + ((uint64_t)op[4] * 12177068667683431813UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 12177068667683431813UL) + ((uint64_t)op[1] * 2398662181628872264UL) + ((uint64_t)op[2] * 13174311866102213863UL) + ((uint64_t)op[3] * 10267894324072613377UL) + ((uint64_t)op[4] * 16157190593066794202UL);
	tmp_q[4] = ((uint64_t)op[0] * 4147259510870461806UL) + ((uint64_t)op[1] * 12177068667683431813UL) + ((uint64_t)op[2] * 2398662181628872264UL) + ((uint64_t)op[3] * 13174311866102213863UL) + ((uint64_t)op[4] * 10267894324072613377UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 177268662894L) - ((((int128)tmp_q[1] * 68431118552L) - ((int128)tmp_q[2] * 123704812217L) - ((int128)tmp_q[3] * 28381495823L) + ((int128)tmp_q[4] * 92310048823L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 92310048823L) - ((int128)tmp_q[1] * 177268662894L) - ((((int128)tmp_q[2] * 68431118552L) - ((int128)tmp_q[3] * 123704812217L) - ((int128)tmp_q[4] * 28381495823L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 28381495823L) + ((int128)tmp_q[1] * 92310048823L) - ((int128)tmp_q[2] * 177268662894L) - ((((int128)tmp_q[3] * 68431118552L) - ((int128)tmp_q[4] * 123704812217L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 123704812217L) - ((int128)tmp_q[1] * 28381495823L) + ((int128)tmp_q[2] * 92310048823L) - ((int128)tmp_q[3] * 177268662894L) - ((int128)tmp_q[4] * 342155592760L);
	tmp_zero[4] = ((int128)tmp_q[0] * 68431118552L) - ((int128)tmp_q[1] * 123704812217L) - ((int128)tmp_q[2] * 28381495823L) + ((int128)tmp_q[3] * 92310048823L) - ((int128)tmp_q[4] * 177268662894L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

