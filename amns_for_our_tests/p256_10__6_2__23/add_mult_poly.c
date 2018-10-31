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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7819930993785315841UL) + ((((uint64_t)op[1] * 10787106121651287333UL) + ((uint64_t)op[2] * 13060004903866626680UL) + ((uint64_t)op[3] * 8156200463887534810UL) + ((uint64_t)op[4] * 6787165478596582077UL) + ((uint64_t)op[5] * 2312172895617522163UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 2312172895617522163UL) + ((uint64_t)op[1] * 7819930993785315841UL) + ((((uint64_t)op[2] * 10787106121651287333UL) + ((uint64_t)op[3] * 13060004903866626680UL) + ((uint64_t)op[4] * 8156200463887534810UL) + ((uint64_t)op[5] * 6787165478596582077UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 6787165478596582077UL) + ((uint64_t)op[1] * 2312172895617522163UL) + ((uint64_t)op[2] * 7819930993785315841UL) + ((((uint64_t)op[3] * 10787106121651287333UL) + ((uint64_t)op[4] * 13060004903866626680UL) + ((uint64_t)op[5] * 8156200463887534810UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 8156200463887534810UL) + ((uint64_t)op[1] * 6787165478596582077UL) + ((uint64_t)op[2] * 2312172895617522163UL) + ((uint64_t)op[3] * 7819930993785315841UL) + ((((uint64_t)op[4] * 10787106121651287333UL) + ((uint64_t)op[5] * 13060004903866626680UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 13060004903866626680UL) + ((uint64_t)op[1] * 8156200463887534810UL) + ((uint64_t)op[2] * 6787165478596582077UL) + ((uint64_t)op[3] * 2312172895617522163UL) + ((uint64_t)op[4] * 7819930993785315841UL) + ((uint64_t)op[5] * 3127468169593023050UL);
	tmp_q[5] = ((uint64_t)op[0] * 10787106121651287333UL) + ((uint64_t)op[1] * 13060004903866626680UL) + ((uint64_t)op[2] * 8156200463887534810UL) + ((uint64_t)op[3] * 6787165478596582077UL) + ((uint64_t)op[4] * 2312172895617522163UL) + ((uint64_t)op[5] * 7819930993785315841UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 181731197569L) + ((((int128)tmp_q[1] * 656354451375L) - ((int128)tmp_q[2] * 3884699439211L) - ((int128)tmp_q[3] * 354273377663L) - ((int128)tmp_q[4] * 2346467812388L) - ((int128)tmp_q[5] * 708158485541L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 708158485541L) + ((int128)tmp_q[1] * 181731197569L) + ((((int128)tmp_q[2] * 656354451375L) - ((int128)tmp_q[3] * 3884699439211L) - ((int128)tmp_q[4] * 354273377663L) - ((int128)tmp_q[5] * 2346467812388L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 2346467812388L) - ((int128)tmp_q[1] * 708158485541L) + ((int128)tmp_q[2] * 181731197569L) + ((((int128)tmp_q[3] * 656354451375L) - ((int128)tmp_q[4] * 3884699439211L) - ((int128)tmp_q[5] * 354273377663L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 354273377663L) - ((int128)tmp_q[1] * 2346467812388L) - ((int128)tmp_q[2] * 708158485541L) + ((int128)tmp_q[3] * 181731197569L) + ((((int128)tmp_q[4] * 656354451375L) - ((int128)tmp_q[5] * 3884699439211L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 3884699439211L) - ((int128)tmp_q[1] * 354273377663L) - ((int128)tmp_q[2] * 2346467812388L) - ((int128)tmp_q[3] * 708158485541L) + ((int128)tmp_q[4] * 181731197569L) + ((int128)tmp_q[5] * 1312708902750L);
	tmp_zero[5] = ((int128)tmp_q[0] * 656354451375L) - ((int128)tmp_q[1] * 3884699439211L) - ((int128)tmp_q[2] * 354273377663L) - ((int128)tmp_q[3] * 2346467812388L) - ((int128)tmp_q[4] * 708158485541L) + ((int128)tmp_q[5] * 181731197569L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

