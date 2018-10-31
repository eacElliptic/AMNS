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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11499951378022298619UL) + ((((uint64_t)op[1] * 12569413876983351304UL) + ((uint64_t)op[2] * 13197185022238243863UL) + ((uint64_t)op[3] * 6826532954814190826UL) + ((uint64_t)op[4] * 17318558190257132933UL) + ((uint64_t)op[5] * 8192503701286449806UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 8192503701286449806UL) + ((uint64_t)op[1] * 11499951378022298619UL) + ((((uint64_t)op[2] * 12569413876983351304UL) + ((uint64_t)op[3] * 13197185022238243863UL) + ((uint64_t)op[4] * 6826532954814190826UL) + ((uint64_t)op[5] * 17318558190257132933UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 17318558190257132933UL) + ((uint64_t)op[1] * 8192503701286449806UL) + ((uint64_t)op[2] * 11499951378022298619UL) + ((((uint64_t)op[3] * 12569413876983351304UL) + ((uint64_t)op[4] * 13197185022238243863UL) + ((uint64_t)op[5] * 6826532954814190826UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 6826532954814190826UL) + ((uint64_t)op[1] * 17318558190257132933UL) + ((uint64_t)op[2] * 8192503701286449806UL) + ((uint64_t)op[3] * 11499951378022298619UL) + ((((uint64_t)op[4] * 12569413876983351304UL) + ((uint64_t)op[5] * 13197185022238243863UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 13197185022238243863UL) + ((uint64_t)op[1] * 6826532954814190826UL) + ((uint64_t)op[2] * 17318558190257132933UL) + ((uint64_t)op[3] * 8192503701286449806UL) + ((uint64_t)op[4] * 11499951378022298619UL) + ((uint64_t)op[5] * 8321590647319052352UL);
	tmp_q[5] = ((uint64_t)op[0] * 12569413876983351304UL) + ((uint64_t)op[1] * 13197185022238243863UL) + ((uint64_t)op[2] * 6826532954814190826UL) + ((uint64_t)op[3] * 17318558190257132933UL) + ((uint64_t)op[4] * 8192503701286449806UL) + ((uint64_t)op[5] * 11499951378022298619UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 20792129022171L) + ((((int128)tmp_q[1] * 50926238520226L) - ((int128)tmp_q[2] * 5061850630528L) - ((int128)tmp_q[3] * 83794408375986L) + ((int128)tmp_q[4] * 44115183359609L) + ((int128)tmp_q[5] * 29165644557934L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 29165644557934L) - ((int128)tmp_q[1] * 20792129022171L) + ((((int128)tmp_q[2] * 50926238520226L) - ((int128)tmp_q[3] * 5061850630528L) - ((int128)tmp_q[4] * 83794408375986L) + ((int128)tmp_q[5] * 44115183359609L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 44115183359609L) + ((int128)tmp_q[1] * 29165644557934L) - ((int128)tmp_q[2] * 20792129022171L) + ((((int128)tmp_q[3] * 50926238520226L) - ((int128)tmp_q[4] * 5061850630528L) - ((int128)tmp_q[5] * 83794408375986L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 83794408375986L) + ((int128)tmp_q[1] * 44115183359609L) + ((int128)tmp_q[2] * 29165644557934L) - ((int128)tmp_q[3] * 20792129022171L) + ((((int128)tmp_q[4] * 50926238520226L) - ((int128)tmp_q[5] * 5061850630528L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 5061850630528L) - ((int128)tmp_q[1] * 83794408375986L) + ((int128)tmp_q[2] * 44115183359609L) + ((int128)tmp_q[3] * 29165644557934L) - ((int128)tmp_q[4] * 20792129022171L) + ((int128)tmp_q[5] * 407409908161808L);
	tmp_zero[5] = ((int128)tmp_q[0] * 50926238520226L) - ((int128)tmp_q[1] * 5061850630528L) - ((int128)tmp_q[2] * 83794408375986L) + ((int128)tmp_q[3] * 44115183359609L) + ((int128)tmp_q[4] * 29165644557934L) - ((int128)tmp_q[5] * 20792129022171L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

