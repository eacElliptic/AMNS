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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4664658345487559701UL) + ((((uint64_t)op[1] * 1871794861150078872UL) + ((uint64_t)op[2] * 5076152652981121950UL) + ((uint64_t)op[3] * 718875892311984093UL) + ((uint64_t)op[4] * 17858806931274470650UL) + ((uint64_t)op[5] * 6531428913276521393UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 6531428913276521393UL) + ((uint64_t)op[1] * 4664658345487559701UL) + ((((uint64_t)op[2] * 1871794861150078872UL) + ((uint64_t)op[3] * 5076152652981121950UL) + ((uint64_t)op[4] * 718875892311984093UL) + ((uint64_t)op[5] * 17858806931274470650UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 17858806931274470650UL) + ((uint64_t)op[1] * 6531428913276521393UL) + ((uint64_t)op[2] * 4664658345487559701UL) + ((((uint64_t)op[3] * 1871794861150078872UL) + ((uint64_t)op[4] * 5076152652981121950UL) + ((uint64_t)op[5] * 718875892311984093UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 718875892311984093UL) + ((uint64_t)op[1] * 17858806931274470650UL) + ((uint64_t)op[2] * 6531428913276521393UL) + ((uint64_t)op[3] * 4664658345487559701UL) + ((((uint64_t)op[4] * 1871794861150078872UL) + ((uint64_t)op[5] * 5076152652981121950UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 5076152652981121950UL) + ((uint64_t)op[1] * 718875892311984093UL) + ((uint64_t)op[2] * 17858806931274470650UL) + ((uint64_t)op[3] * 6531428913276521393UL) + ((uint64_t)op[4] * 4664658345487559701UL) + ((uint64_t)op[5] * 13102564028050552104UL);
	tmp_q[5] = ((uint64_t)op[0] * 1871794861150078872UL) + ((uint64_t)op[1] * 5076152652981121950UL) + ((uint64_t)op[2] * 718875892311984093UL) + ((uint64_t)op[3] * 17858806931274470650UL) + ((uint64_t)op[4] * 6531428913276521393UL) + ((uint64_t)op[5] * 4664658345487559701UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1223170263552L) + ((((int128)tmp_q[1] * 1605316050831L) + ((int128)tmp_q[2] * 446783680299L) - ((int128)tmp_q[3] * 1621972848820L) - ((int128)tmp_q[4] * 4258808734620L) - ((int128)tmp_q[5] * 2470277787055L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 2470277787055L) + ((int128)tmp_q[1] * 1223170263552L) + ((((int128)tmp_q[2] * 1605316050831L) + ((int128)tmp_q[3] * 446783680299L) - ((int128)tmp_q[4] * 1621972848820L) - ((int128)tmp_q[5] * 4258808734620L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 4258808734620L) - ((int128)tmp_q[1] * 2470277787055L) + ((int128)tmp_q[2] * 1223170263552L) + ((((int128)tmp_q[3] * 1605316050831L) + ((int128)tmp_q[4] * 446783680299L) - ((int128)tmp_q[5] * 1621972848820L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 1621972848820L) - ((int128)tmp_q[1] * 4258808734620L) - ((int128)tmp_q[2] * 2470277787055L) + ((int128)tmp_q[3] * 1223170263552L) + ((((int128)tmp_q[4] * 1605316050831L) + ((int128)tmp_q[5] * 446783680299L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 446783680299L) - ((int128)tmp_q[1] * 1621972848820L) - ((int128)tmp_q[2] * 4258808734620L) - ((int128)tmp_q[3] * 2470277787055L) + ((int128)tmp_q[4] * 1223170263552L) + ((int128)tmp_q[5] * 11237212355817L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1605316050831L) + ((int128)tmp_q[1] * 446783680299L) - ((int128)tmp_q[2] * 1621972848820L) - ((int128)tmp_q[3] * 4258808734620L) - ((int128)tmp_q[4] * 2470277787055L) + ((int128)tmp_q[5] * 1223170263552L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

