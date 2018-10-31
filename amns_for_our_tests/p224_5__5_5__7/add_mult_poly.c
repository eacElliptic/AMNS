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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18114174157671867153UL) + ((((uint64_t)op[1] * 9844831334997105866UL) + ((uint64_t)op[2] * 1665871752822393707UL) + ((uint64_t)op[3] * 1501150661744454193UL) + ((uint64_t)op[4] * 11694318740217316340UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 11694318740217316340UL) + ((uint64_t)op[1] * 18114174157671867153UL) + ((((uint64_t)op[2] * 9844831334997105866UL) + ((uint64_t)op[3] * 1665871752822393707UL) + ((uint64_t)op[4] * 1501150661744454193UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 1501150661744454193UL) + ((uint64_t)op[1] * 11694318740217316340UL) + ((uint64_t)op[2] * 18114174157671867153UL) + ((((uint64_t)op[3] * 9844831334997105866UL) + ((uint64_t)op[4] * 1665871752822393707UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 1665871752822393707UL) + ((uint64_t)op[1] * 1501150661744454193UL) + ((uint64_t)op[2] * 11694318740217316340UL) + ((uint64_t)op[3] * 18114174157671867153UL) + ((uint64_t)op[4] * 12330668527566426098UL);
	tmp_q[4] = ((uint64_t)op[0] * 9844831334997105866UL) + ((uint64_t)op[1] * 1665871752822393707UL) + ((uint64_t)op[2] * 1501150661744454193UL) + ((uint64_t)op[3] * 11694318740217316340UL) + ((uint64_t)op[4] * 18114174157671867153UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2255310445987L) + ((-((int128)tmp_q[1] * 16699786093943L) + ((int128)tmp_q[2] * 21278422845272L) - ((int128)tmp_q[3] * 8696241007800L) - ((int128)tmp_q[4] * 6029662073897L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 6029662073897L) - ((int128)tmp_q[1] * 2255310445987L) + ((-((int128)tmp_q[2] * 16699786093943L) + ((int128)tmp_q[3] * 21278422845272L) - ((int128)tmp_q[4] * 8696241007800L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 8696241007800L) - ((int128)tmp_q[1] * 6029662073897L) - ((int128)tmp_q[2] * 2255310445987L) + ((-((int128)tmp_q[3] * 16699786093943L) + ((int128)tmp_q[4] * 21278422845272L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 21278422845272L) - ((int128)tmp_q[1] * 8696241007800L) - ((int128)tmp_q[2] * 6029662073897L) - ((int128)tmp_q[3] * 2255310445987L) - ((int128)tmp_q[4] * 83498930469715L);
	tmp_zero[4] = -((int128)tmp_q[0] * 16699786093943L) + ((int128)tmp_q[1] * 21278422845272L) - ((int128)tmp_q[2] * 8696241007800L) - ((int128)tmp_q[3] * 6029662073897L) - ((int128)tmp_q[4] * 2255310445987L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

