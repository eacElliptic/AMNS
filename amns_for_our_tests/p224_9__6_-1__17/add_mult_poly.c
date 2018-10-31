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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11328228932280001837UL) + ((((uint64_t)op[1] * 12769308476718936333UL) + ((uint64_t)op[2] * 11527209125458483320UL) + ((uint64_t)op[3] * 920124495266367041UL) + ((uint64_t)op[4] * 5493796576600083636UL) + ((uint64_t)op[5] * 2257149721972381944UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 2257149721972381944UL) + ((uint64_t)op[1] * 11328228932280001837UL) + ((((uint64_t)op[2] * 12769308476718936333UL) + ((uint64_t)op[3] * 11527209125458483320UL) + ((uint64_t)op[4] * 920124495266367041UL) + ((uint64_t)op[5] * 5493796576600083636UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 5493796576600083636UL) + ((uint64_t)op[1] * 2257149721972381944UL) + ((uint64_t)op[2] * 11328228932280001837UL) + ((((uint64_t)op[3] * 12769308476718936333UL) + ((uint64_t)op[4] * 11527209125458483320UL) + ((uint64_t)op[5] * 920124495266367041UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 920124495266367041UL) + ((uint64_t)op[1] * 5493796576600083636UL) + ((uint64_t)op[2] * 2257149721972381944UL) + ((uint64_t)op[3] * 11328228932280001837UL) + ((((uint64_t)op[4] * 12769308476718936333UL) + ((uint64_t)op[5] * 11527209125458483320UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 11527209125458483320UL) + ((uint64_t)op[1] * 920124495266367041UL) + ((uint64_t)op[2] * 5493796576600083636UL) + ((uint64_t)op[3] * 2257149721972381944UL) + ((uint64_t)op[4] * 11328228932280001837UL) + ((uint64_t)op[5] * 5677435596990615283UL);
	tmp_q[5] = ((uint64_t)op[0] * 12769308476718936333UL) + ((uint64_t)op[1] * 11527209125458483320UL) + ((uint64_t)op[2] * 920124495266367041UL) + ((uint64_t)op[3] * 5493796576600083636UL) + ((uint64_t)op[4] * 2257149721972381944UL) + ((uint64_t)op[5] * 11328228932280001837UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5894773318506378L) - (-((int128)tmp_q[1] * 41589546992545483L) + ((int128)tmp_q[2] * 14428876396087386L) - ((int128)tmp_q[3] * 18887375538047620L) + ((int128)tmp_q[4] * 20323649714593765L) + ((int128)tmp_q[5] * 22702171454497867L));
	tmp_zero[1] = ((int128)tmp_q[0] * 22702171454497867L) + ((int128)tmp_q[1] * 5894773318506378L) - (-((int128)tmp_q[2] * 41589546992545483L) + ((int128)tmp_q[3] * 14428876396087386L) - ((int128)tmp_q[4] * 18887375538047620L) + ((int128)tmp_q[5] * 20323649714593765L));
	tmp_zero[2] = ((int128)tmp_q[0] * 20323649714593765L) + ((int128)tmp_q[1] * 22702171454497867L) + ((int128)tmp_q[2] * 5894773318506378L) - (-((int128)tmp_q[3] * 41589546992545483L) + ((int128)tmp_q[4] * 14428876396087386L) - ((int128)tmp_q[5] * 18887375538047620L));
	tmp_zero[3] = -((int128)tmp_q[0] * 18887375538047620L) + ((int128)tmp_q[1] * 20323649714593765L) + ((int128)tmp_q[2] * 22702171454497867L) + ((int128)tmp_q[3] * 5894773318506378L) - (-((int128)tmp_q[4] * 41589546992545483L) + ((int128)tmp_q[5] * 14428876396087386L));
	tmp_zero[4] = ((int128)tmp_q[0] * 14428876396087386L) - ((int128)tmp_q[1] * 18887375538047620L) + ((int128)tmp_q[2] * 20323649714593765L) + ((int128)tmp_q[3] * 22702171454497867L) + ((int128)tmp_q[4] * 5894773318506378L) + ((int128)tmp_q[5] * 41589546992545483L);
	tmp_zero[5] = -((int128)tmp_q[0] * 41589546992545483L) + ((int128)tmp_q[1] * 14428876396087386L) - ((int128)tmp_q[2] * 18887375538047620L) + ((int128)tmp_q[3] * 20323649714593765L) + ((int128)tmp_q[4] * 22702171454497867L) + ((int128)tmp_q[5] * 5894773318506378L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

