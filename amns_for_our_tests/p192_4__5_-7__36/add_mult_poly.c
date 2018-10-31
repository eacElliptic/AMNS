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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13807273087295615040UL) + ((((uint64_t)op[1] * 5238950903379897927UL) + ((uint64_t)op[2] * 3696578203269834912UL) + ((uint64_t)op[3] * 10943876834641849847UL) + ((uint64_t)op[4] * 16110340301201372303UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 16110340301201372303UL) + ((uint64_t)op[1] * 13807273087295615040UL) + ((((uint64_t)op[2] * 5238950903379897927UL) + ((uint64_t)op[3] * 3696578203269834912UL) + ((uint64_t)op[4] * 10943876834641849847UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 10943876834641849847UL) + ((uint64_t)op[1] * 16110340301201372303UL) + ((uint64_t)op[2] * 13807273087295615040UL) + ((((uint64_t)op[3] * 5238950903379897927UL) + ((uint64_t)op[4] * 3696578203269834912UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 3696578203269834912UL) + ((uint64_t)op[1] * 10943876834641849847UL) + ((uint64_t)op[2] * 16110340301201372303UL) + ((uint64_t)op[3] * 13807273087295615040UL) + ((uint64_t)op[4] * 220831823759817743UL);
	tmp_q[4] = ((uint64_t)op[0] * 5238950903379897927UL) + ((uint64_t)op[1] * 3696578203269834912UL) + ((uint64_t)op[2] * 10943876834641849847UL) + ((uint64_t)op[3] * 16110340301201372303UL) + ((uint64_t)op[4] * 13807273087295615040UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 119355461901L) - ((((int128)tmp_q[1] * 217041929734L) + ((int128)tmp_q[2] * 90186881190L) - ((int128)tmp_q[3] * 47341138547L) - ((int128)tmp_q[4] * 79809183979L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 79809183979L) + ((int128)tmp_q[1] * 119355461901L) - ((((int128)tmp_q[2] * 217041929734L) + ((int128)tmp_q[3] * 90186881190L) - ((int128)tmp_q[4] * 47341138547L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 47341138547L) - ((int128)tmp_q[1] * 79809183979L) + ((int128)tmp_q[2] * 119355461901L) - ((((int128)tmp_q[3] * 217041929734L) + ((int128)tmp_q[4] * 90186881190L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 90186881190L) - ((int128)tmp_q[1] * 47341138547L) - ((int128)tmp_q[2] * 79809183979L) + ((int128)tmp_q[3] * 119355461901L) - ((int128)tmp_q[4] * 1519293508138L);
	tmp_zero[4] = ((int128)tmp_q[0] * 217041929734L) + ((int128)tmp_q[1] * 90186881190L) - ((int128)tmp_q[2] * 47341138547L) - ((int128)tmp_q[3] * 79809183979L) + ((int128)tmp_q[4] * 119355461901L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

