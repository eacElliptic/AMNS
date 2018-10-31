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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 770394644257536378UL) + ((((uint64_t)op[1] * 9899834421608492143UL) + ((uint64_t)op[2] * 17725403233001841257UL) + ((uint64_t)op[3] * 17726083336730706744UL) + ((uint64_t)op[4] * 7505347849104820664UL) + ((uint64_t)op[5] * 13877710500065099983UL) + ((uint64_t)op[6] * 16293031365414220272UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 16293031365414220272UL) + ((uint64_t)op[1] * 770394644257536378UL) + ((((uint64_t)op[2] * 9899834421608492143UL) + ((uint64_t)op[3] * 17725403233001841257UL) + ((uint64_t)op[4] * 17726083336730706744UL) + ((uint64_t)op[5] * 7505347849104820664UL) + ((uint64_t)op[6] * 13877710500065099983UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13877710500065099983UL) + ((uint64_t)op[1] * 16293031365414220272UL) + ((uint64_t)op[2] * 770394644257536378UL) + ((((uint64_t)op[3] * 9899834421608492143UL) + ((uint64_t)op[4] * 17725403233001841257UL) + ((uint64_t)op[5] * 17726083336730706744UL) + ((uint64_t)op[6] * 7505347849104820664UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 7505347849104820664UL) + ((uint64_t)op[1] * 13877710500065099983UL) + ((uint64_t)op[2] * 16293031365414220272UL) + ((uint64_t)op[3] * 770394644257536378UL) + ((((uint64_t)op[4] * 9899834421608492143UL) + ((uint64_t)op[5] * 17725403233001841257UL) + ((uint64_t)op[6] * 17726083336730706744UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 17726083336730706744UL) + ((uint64_t)op[1] * 7505347849104820664UL) + ((uint64_t)op[2] * 13877710500065099983UL) + ((uint64_t)op[3] * 16293031365414220272UL) + ((uint64_t)op[4] * 770394644257536378UL) + ((((uint64_t)op[5] * 9899834421608492143UL) + ((uint64_t)op[6] * 17725403233001841257UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 17725403233001841257UL) + ((uint64_t)op[1] * 17726083336730706744UL) + ((uint64_t)op[2] * 7505347849104820664UL) + ((uint64_t)op[3] * 13877710500065099983UL) + ((uint64_t)op[4] * 16293031365414220272UL) + ((uint64_t)op[5] * 770394644257536378UL) + ((uint64_t)op[6] * 12605683960623357483UL);
	tmp_q[6] = ((uint64_t)op[0] * 9899834421608492143UL) + ((uint64_t)op[1] * 17725403233001841257UL) + ((uint64_t)op[2] * 17726083336730706744UL) + ((uint64_t)op[3] * 7505347849104820664UL) + ((uint64_t)op[4] * 13877710500065099983UL) + ((uint64_t)op[5] * 16293031365414220272UL) + ((uint64_t)op[6] * 770394644257536378UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3783640219L) + ((-((int128)tmp_q[1] * 44139697576L) - ((int128)tmp_q[2] * 23100733407L) - ((int128)tmp_q[3] * 27644193090L) - ((int128)tmp_q[4] * 606608315L) + ((int128)tmp_q[5] * 53943600399L) + ((int128)tmp_q[6] * 49280666113L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 49280666113L) + ((int128)tmp_q[1] * 3783640219L) + ((-((int128)tmp_q[2] * 44139697576L) - ((int128)tmp_q[3] * 23100733407L) - ((int128)tmp_q[4] * 27644193090L) - ((int128)tmp_q[5] * 606608315L) + ((int128)tmp_q[6] * 53943600399L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 53943600399L) + ((int128)tmp_q[1] * 49280666113L) + ((int128)tmp_q[2] * 3783640219L) + ((-((int128)tmp_q[3] * 44139697576L) - ((int128)tmp_q[4] * 23100733407L) - ((int128)tmp_q[5] * 27644193090L) - ((int128)tmp_q[6] * 606608315L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 606608315L) + ((int128)tmp_q[1] * 53943600399L) + ((int128)tmp_q[2] * 49280666113L) + ((int128)tmp_q[3] * 3783640219L) + ((-((int128)tmp_q[4] * 44139697576L) - ((int128)tmp_q[5] * 23100733407L) - ((int128)tmp_q[6] * 27644193090L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 27644193090L) - ((int128)tmp_q[1] * 606608315L) + ((int128)tmp_q[2] * 53943600399L) + ((int128)tmp_q[3] * 49280666113L) + ((int128)tmp_q[4] * 3783640219L) + ((-((int128)tmp_q[5] * 44139697576L) - ((int128)tmp_q[6] * 23100733407L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 23100733407L) - ((int128)tmp_q[1] * 27644193090L) - ((int128)tmp_q[2] * 606608315L) + ((int128)tmp_q[3] * 53943600399L) + ((int128)tmp_q[4] * 49280666113L) + ((int128)tmp_q[5] * 3783640219L) - ((int128)tmp_q[6] * 220698487880L);
	tmp_zero[6] = -((int128)tmp_q[0] * 44139697576L) - ((int128)tmp_q[1] * 23100733407L) - ((int128)tmp_q[2] * 27644193090L) - ((int128)tmp_q[3] * 606608315L) + ((int128)tmp_q[4] * 53943600399L) + ((int128)tmp_q[5] * 49280666113L) + ((int128)tmp_q[6] * 3783640219L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

