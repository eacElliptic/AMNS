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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7739017959692478316UL) + ((((uint64_t)op[1] * 17324124284988689557UL) + ((uint64_t)op[2] * 6784552652182239995UL) + ((uint64_t)op[3] * 375355679208724421UL) + ((uint64_t)op[4] * 3647049000892604615UL) + ((uint64_t)op[5] * 10432664229507985997UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 10432664229507985997UL) + ((uint64_t)op[1] * 7739017959692478316UL) + ((((uint64_t)op[2] * 17324124284988689557UL) + ((uint64_t)op[3] * 6784552652182239995UL) + ((uint64_t)op[4] * 375355679208724421UL) + ((uint64_t)op[5] * 3647049000892604615UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 3647049000892604615UL) + ((uint64_t)op[1] * 10432664229507985997UL) + ((uint64_t)op[2] * 7739017959692478316UL) + ((((uint64_t)op[3] * 17324124284988689557UL) + ((uint64_t)op[4] * 6784552652182239995UL) + ((uint64_t)op[5] * 375355679208724421UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 375355679208724421UL) + ((uint64_t)op[1] * 3647049000892604615UL) + ((uint64_t)op[2] * 10432664229507985997UL) + ((uint64_t)op[3] * 7739017959692478316UL) + ((((uint64_t)op[4] * 17324124284988689557UL) + ((uint64_t)op[5] * 6784552652182239995UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 6784552652182239995UL) + ((uint64_t)op[1] * 375355679208724421UL) + ((uint64_t)op[2] * 3647049000892604615UL) + ((uint64_t)op[3] * 10432664229507985997UL) + ((uint64_t)op[4] * 7739017959692478316UL) + ((uint64_t)op[5] * 5613098943604310295UL);
	tmp_q[5] = ((uint64_t)op[0] * 17324124284988689557UL) + ((uint64_t)op[1] * 6784552652182239995UL) + ((uint64_t)op[2] * 375355679208724421UL) + ((uint64_t)op[3] * 3647049000892604615UL) + ((uint64_t)op[4] * 10432664229507985997UL) + ((uint64_t)op[5] * 7739017959692478316UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1558661244L) - ((-((int128)tmp_q[1] * 2771191791L) + ((int128)tmp_q[2] * 929956053L) + ((int128)tmp_q[3] * 68219577L) - ((int128)tmp_q[4] * 2007232951L) - ((int128)tmp_q[5] * 185904023L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 185904023L) - ((int128)tmp_q[1] * 1558661244L) - ((-((int128)tmp_q[2] * 2771191791L) + ((int128)tmp_q[3] * 929956053L) + ((int128)tmp_q[4] * 68219577L) - ((int128)tmp_q[5] * 2007232951L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 2007232951L) - ((int128)tmp_q[1] * 185904023L) - ((int128)tmp_q[2] * 1558661244L) - ((-((int128)tmp_q[3] * 2771191791L) + ((int128)tmp_q[4] * 929956053L) + ((int128)tmp_q[5] * 68219577L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 68219577L) - ((int128)tmp_q[1] * 2007232951L) - ((int128)tmp_q[2] * 185904023L) - ((int128)tmp_q[3] * 1558661244L) - ((-((int128)tmp_q[4] * 2771191791L) + ((int128)tmp_q[5] * 929956053L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 929956053L) + ((int128)tmp_q[1] * 68219577L) - ((int128)tmp_q[2] * 2007232951L) - ((int128)tmp_q[3] * 185904023L) - ((int128)tmp_q[4] * 1558661244L) + ((int128)tmp_q[5] * 13855958955L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2771191791L) + ((int128)tmp_q[1] * 929956053L) + ((int128)tmp_q[2] * 68219577L) - ((int128)tmp_q[3] * 2007232951L) - ((int128)tmp_q[4] * 185904023L) - ((int128)tmp_q[5] * 1558661244L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

