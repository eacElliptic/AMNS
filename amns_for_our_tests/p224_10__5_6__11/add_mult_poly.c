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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2245332225466513347UL) + ((((uint64_t)op[1] * 10616681384960371817UL) + ((uint64_t)op[2] * 16040829461017354903UL) + ((uint64_t)op[3] * 11247942867118418561UL) + ((uint64_t)op[4] * 10682835023227501825UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 10682835023227501825UL) + ((uint64_t)op[1] * 2245332225466513347UL) + ((((uint64_t)op[2] * 10616681384960371817UL) + ((uint64_t)op[3] * 16040829461017354903UL) + ((uint64_t)op[4] * 11247942867118418561UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 11247942867118418561UL) + ((uint64_t)op[1] * 10682835023227501825UL) + ((uint64_t)op[2] * 2245332225466513347UL) + ((((uint64_t)op[3] * 10616681384960371817UL) + ((uint64_t)op[4] * 16040829461017354903UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 16040829461017354903UL) + ((uint64_t)op[1] * 11247942867118418561UL) + ((uint64_t)op[2] * 10682835023227501825UL) + ((uint64_t)op[3] * 2245332225466513347UL) + ((uint64_t)op[4] * 8359856088633576054UL);
	tmp_q[4] = ((uint64_t)op[0] * 10616681384960371817UL) + ((uint64_t)op[1] * 16040829461017354903UL) + ((uint64_t)op[2] * 11247942867118418561UL) + ((uint64_t)op[3] * 10682835023227501825UL) + ((uint64_t)op[4] * 2245332225466513347UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 9745496897727L) + ((-((int128)tmp_q[1] * 222122512988L) + ((int128)tmp_q[2] * 17423387609402L) + ((int128)tmp_q[3] * 13821851695110L) + ((int128)tmp_q[4] * 450542665995L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 450542665995L) + ((int128)tmp_q[1] * 9745496897727L) + ((-((int128)tmp_q[2] * 222122512988L) + ((int128)tmp_q[3] * 17423387609402L) + ((int128)tmp_q[4] * 13821851695110L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 13821851695110L) + ((int128)tmp_q[1] * 450542665995L) + ((int128)tmp_q[2] * 9745496897727L) + ((-((int128)tmp_q[3] * 222122512988L) + ((int128)tmp_q[4] * 17423387609402L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 17423387609402L) + ((int128)tmp_q[1] * 13821851695110L) + ((int128)tmp_q[2] * 450542665995L) + ((int128)tmp_q[3] * 9745496897727L) - ((int128)tmp_q[4] * 1332735077928L);
	tmp_zero[4] = -((int128)tmp_q[0] * 222122512988L) + ((int128)tmp_q[1] * 17423387609402L) + ((int128)tmp_q[2] * 13821851695110L) + ((int128)tmp_q[3] * 450542665995L) + ((int128)tmp_q[4] * 9745496897727L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

