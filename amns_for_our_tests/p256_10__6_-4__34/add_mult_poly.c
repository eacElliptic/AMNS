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
	tmp_q[0] = ((uint64_t)op[0] * 14632954802491846113UL) + ((((uint64_t)op[1] * 15578335808725696243UL) + ((uint64_t)op[2] * 6065648274945352460UL) + ((uint64_t)op[3] * 16095490962372666641UL) + ((uint64_t)op[4] * 6121351022196209543UL) + ((uint64_t)op[5] * 13967241133428022967UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 13967241133428022967UL) + ((uint64_t)op[1] * 14632954802491846113UL) + ((((uint64_t)op[2] * 15578335808725696243UL) + ((uint64_t)op[3] * 6065648274945352460UL) + ((uint64_t)op[4] * 16095490962372666641UL) + ((uint64_t)op[5] * 6121351022196209543UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 6121351022196209543UL) + ((uint64_t)op[1] * 13967241133428022967UL) + ((uint64_t)op[2] * 14632954802491846113UL) + ((((uint64_t)op[3] * 15578335808725696243UL) + ((uint64_t)op[4] * 6065648274945352460UL) + ((uint64_t)op[5] * 16095490962372666641UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 16095490962372666641UL) + ((uint64_t)op[1] * 6121351022196209543UL) + ((uint64_t)op[2] * 13967241133428022967UL) + ((uint64_t)op[3] * 14632954802491846113UL) + ((((uint64_t)op[4] * 15578335808725696243UL) + ((uint64_t)op[5] * 6065648274945352460UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 6065648274945352460UL) + ((uint64_t)op[1] * 16095490962372666641UL) + ((uint64_t)op[2] * 6121351022196209543UL) + ((uint64_t)op[3] * 13967241133428022967UL) + ((uint64_t)op[4] * 14632954802491846113UL) + ((uint64_t)op[5] * 11473633059935421492UL);
	tmp_q[5] = ((uint64_t)op[0] * 15578335808725696243UL) + ((uint64_t)op[1] * 6065648274945352460UL) + ((uint64_t)op[2] * 16095490962372666641UL) + ((uint64_t)op[3] * 6121351022196209543UL) + ((uint64_t)op[4] * 13967241133428022967UL) + ((uint64_t)op[5] * 14632954802491846113UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1480730806233L) - ((-((int128)tmp_q[1] * 626956318172L) + ((int128)tmp_q[2] * 3562199235773L) + ((int128)tmp_q[3] * 4948332377106L) - ((int128)tmp_q[4] * 1174140353190L) + ((int128)tmp_q[5] * 4867805273787L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 4867805273787L) - ((int128)tmp_q[1] * 1480730806233L) - ((-((int128)tmp_q[2] * 626956318172L) + ((int128)tmp_q[3] * 3562199235773L) + ((int128)tmp_q[4] * 4948332377106L) - ((int128)tmp_q[5] * 1174140353190L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 1174140353190L) + ((int128)tmp_q[1] * 4867805273787L) - ((int128)tmp_q[2] * 1480730806233L) - ((-((int128)tmp_q[3] * 626956318172L) + ((int128)tmp_q[4] * 3562199235773L) + ((int128)tmp_q[5] * 4948332377106L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 4948332377106L) - ((int128)tmp_q[1] * 1174140353190L) + ((int128)tmp_q[2] * 4867805273787L) - ((int128)tmp_q[3] * 1480730806233L) - ((-((int128)tmp_q[4] * 626956318172L) + ((int128)tmp_q[5] * 3562199235773L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 3562199235773L) + ((int128)tmp_q[1] * 4948332377106L) - ((int128)tmp_q[2] * 1174140353190L) + ((int128)tmp_q[3] * 4867805273787L) - ((int128)tmp_q[4] * 1480730806233L) + ((int128)tmp_q[5] * 2507825272688L);
	tmp_zero[5] = -((int128)tmp_q[0] * 626956318172L) + ((int128)tmp_q[1] * 3562199235773L) + ((int128)tmp_q[2] * 4948332377106L) - ((int128)tmp_q[3] * 1174140353190L) + ((int128)tmp_q[4] * 4867805273787L) - ((int128)tmp_q[5] * 1480730806233L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

