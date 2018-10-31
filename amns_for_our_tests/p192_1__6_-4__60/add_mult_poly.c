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
	tmp_q[0] = ((uint64_t)op[0] * 4431761109422774153UL) + ((((uint64_t)op[1] * 8115196317886924197UL) + ((uint64_t)op[2] * 16953861090145866032UL) + ((uint64_t)op[3] * 3791169927398002113UL) + ((uint64_t)op[4] * 656929817170253670UL) + ((uint64_t)op[5] * 9290904232011442815UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 9290904232011442815UL) + ((uint64_t)op[1] * 4431761109422774153UL) + ((((uint64_t)op[2] * 8115196317886924197UL) + ((uint64_t)op[3] * 16953861090145866032UL) + ((uint64_t)op[4] * 3791169927398002113UL) + ((uint64_t)op[5] * 656929817170253670UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 656929817170253670UL) + ((uint64_t)op[1] * 9290904232011442815UL) + ((uint64_t)op[2] * 4431761109422774153UL) + ((((uint64_t)op[3] * 8115196317886924197UL) + ((uint64_t)op[4] * 16953861090145866032UL) + ((uint64_t)op[5] * 3791169927398002113UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 3791169927398002113UL) + ((uint64_t)op[1] * 656929817170253670UL) + ((uint64_t)op[2] * 9290904232011442815UL) + ((uint64_t)op[3] * 4431761109422774153UL) + ((((uint64_t)op[4] * 8115196317886924197UL) + ((uint64_t)op[5] * 16953861090145866032UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 16953861090145866032UL) + ((uint64_t)op[1] * 3791169927398002113UL) + ((uint64_t)op[2] * 656929817170253670UL) + ((uint64_t)op[3] * 9290904232011442815UL) + ((uint64_t)op[4] * 4431761109422774153UL) + ((uint64_t)op[5] * 4432702875871406444UL);
	tmp_q[5] = ((uint64_t)op[0] * 8115196317886924197UL) + ((uint64_t)op[1] * 16953861090145866032UL) + ((uint64_t)op[2] * 3791169927398002113UL) + ((uint64_t)op[3] * 656929817170253670UL) + ((uint64_t)op[4] * 9290904232011442815UL) + ((uint64_t)op[5] * 4431761109422774153UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2142165135L) - ((-((int128)tmp_q[1] * 1488962073L) + ((int128)tmp_q[2] * 1052034731L) - ((int128)tmp_q[3] * 2010983172L) + ((int128)tmp_q[4] * 136558953L) - ((int128)tmp_q[5] * 1484971745L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 1484971745L) + ((int128)tmp_q[1] * 2142165135L) - ((-((int128)tmp_q[2] * 1488962073L) + ((int128)tmp_q[3] * 1052034731L) - ((int128)tmp_q[4] * 2010983172L) + ((int128)tmp_q[5] * 136558953L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 136558953L) - ((int128)tmp_q[1] * 1484971745L) + ((int128)tmp_q[2] * 2142165135L) - ((-((int128)tmp_q[3] * 1488962073L) + ((int128)tmp_q[4] * 1052034731L) - ((int128)tmp_q[5] * 2010983172L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 2010983172L) + ((int128)tmp_q[1] * 136558953L) - ((int128)tmp_q[2] * 1484971745L) + ((int128)tmp_q[3] * 2142165135L) - ((-((int128)tmp_q[4] * 1488962073L) + ((int128)tmp_q[5] * 1052034731L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 1052034731L) - ((int128)tmp_q[1] * 2010983172L) + ((int128)tmp_q[2] * 136558953L) - ((int128)tmp_q[3] * 1484971745L) + ((int128)tmp_q[4] * 2142165135L) + ((int128)tmp_q[5] * 5955848292L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1488962073L) + ((int128)tmp_q[1] * 1052034731L) - ((int128)tmp_q[2] * 2010983172L) + ((int128)tmp_q[3] * 136558953L) - ((int128)tmp_q[4] * 1484971745L) + ((int128)tmp_q[5] * 2142165135L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

