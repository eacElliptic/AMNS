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
	tmp_q[0] = ((uint64_t)op[0] * 17770082284859449245UL) + ((((uint64_t)op[1] * 17997397546819325845UL) + ((uint64_t)op[2] * 15270539488762462126UL) + ((uint64_t)op[3] * 8974994113041499571UL) + ((uint64_t)op[4] * 17480507358758609282UL) + ((uint64_t)op[5] * 11688878168603761036UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 11688878168603761036UL) + ((uint64_t)op[1] * 17770082284859449245UL) + ((((uint64_t)op[2] * 17997397546819325845UL) + ((uint64_t)op[3] * 15270539488762462126UL) + ((uint64_t)op[4] * 8974994113041499571UL) + ((uint64_t)op[5] * 17480507358758609282UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 17480507358758609282UL) + ((uint64_t)op[1] * 11688878168603761036UL) + ((uint64_t)op[2] * 17770082284859449245UL) + ((((uint64_t)op[3] * 17997397546819325845UL) + ((uint64_t)op[4] * 15270539488762462126UL) + ((uint64_t)op[5] * 8974994113041499571UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 8974994113041499571UL) + ((uint64_t)op[1] * 17480507358758609282UL) + ((uint64_t)op[2] * 11688878168603761036UL) + ((uint64_t)op[3] * 17770082284859449245UL) + ((((uint64_t)op[4] * 17997397546819325845UL) + ((uint64_t)op[5] * 15270539488762462126UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 15270539488762462126UL) + ((uint64_t)op[1] * 8974994113041499571UL) + ((uint64_t)op[2] * 17480507358758609282UL) + ((uint64_t)op[3] * 11688878168603761036UL) + ((uint64_t)op[4] * 17770082284859449245UL) + ((uint64_t)op[5] * 1797386107560903084UL);
	tmp_q[5] = ((uint64_t)op[0] * 17997397546819325845UL) + ((uint64_t)op[1] * 15270539488762462126UL) + ((uint64_t)op[2] * 8974994113041499571UL) + ((uint64_t)op[3] * 17480507358758609282UL) + ((uint64_t)op[4] * 11688878168603761036UL) + ((uint64_t)op[5] * 17770082284859449245UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 103796917617L) - ((((int128)tmp_q[1] * 17144614453L) + ((int128)tmp_q[2] * 53853933342L) - ((int128)tmp_q[3] * 31868975681L) - ((int128)tmp_q[4] * 116191461358L) + ((int128)tmp_q[5] * 38996987388L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 38996987388L) - ((int128)tmp_q[1] * 103796917617L) - ((((int128)tmp_q[2] * 17144614453L) + ((int128)tmp_q[3] * 53853933342L) - ((int128)tmp_q[4] * 31868975681L) - ((int128)tmp_q[5] * 116191461358L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 116191461358L) + ((int128)tmp_q[1] * 38996987388L) - ((int128)tmp_q[2] * 103796917617L) - ((((int128)tmp_q[3] * 17144614453L) + ((int128)tmp_q[4] * 53853933342L) - ((int128)tmp_q[5] * 31868975681L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 31868975681L) - ((int128)tmp_q[1] * 116191461358L) + ((int128)tmp_q[2] * 38996987388L) - ((int128)tmp_q[3] * 103796917617L) - ((((int128)tmp_q[4] * 17144614453L) + ((int128)tmp_q[5] * 53853933342L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 53853933342L) - ((int128)tmp_q[1] * 31868975681L) - ((int128)tmp_q[2] * 116191461358L) + ((int128)tmp_q[3] * 38996987388L) - ((int128)tmp_q[4] * 103796917617L) - ((int128)tmp_q[5] * 68578457812L);
	tmp_zero[5] = ((int128)tmp_q[0] * 17144614453L) + ((int128)tmp_q[1] * 53853933342L) - ((int128)tmp_q[2] * 31868975681L) - ((int128)tmp_q[3] * 116191461358L) + ((int128)tmp_q[4] * 38996987388L) - ((int128)tmp_q[5] * 103796917617L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

