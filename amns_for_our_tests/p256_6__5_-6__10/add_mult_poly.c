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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5809160881786028891UL) + ((((uint64_t)op[1] * 13146885663826085357UL) + ((uint64_t)op[2] * 1722234353254358194UL) + ((uint64_t)op[3] * 14989815789129363828UL) + ((uint64_t)op[4] * 10580047948882452637UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 10580047948882452637UL) + ((uint64_t)op[1] * 5809160881786028891UL) + ((((uint64_t)op[2] * 13146885663826085357UL) + ((uint64_t)op[3] * 1722234353254358194UL) + ((uint64_t)op[4] * 14989815789129363828UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 14989815789129363828UL) + ((uint64_t)op[1] * 10580047948882452637UL) + ((uint64_t)op[2] * 5809160881786028891UL) + ((((uint64_t)op[3] * 13146885663826085357UL) + ((uint64_t)op[4] * 1722234353254358194UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 1722234353254358194UL) + ((uint64_t)op[1] * 14989815789129363828UL) + ((uint64_t)op[2] * 10580047948882452637UL) + ((uint64_t)op[3] * 5809160881786028891UL) + ((uint64_t)op[4] * 13352406385591245938UL);
	tmp_q[4] = ((uint64_t)op[0] * 13146885663826085357UL) + ((uint64_t)op[1] * 1722234353254358194UL) + ((uint64_t)op[2] * 14989815789129363828UL) + ((uint64_t)op[3] * 10580047948882452637UL) + ((uint64_t)op[4] * 5809160881786028891UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1461880424685217L) - ((((int128)tmp_q[1] * 1486579628267562L) - ((int128)tmp_q[2] * 458449968195883L) + ((int128)tmp_q[3] * 1119059807217263L) + ((int128)tmp_q[4] * 106375322864705L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 106375322864705L) - ((int128)tmp_q[1] * 1461880424685217L) - ((((int128)tmp_q[2] * 1486579628267562L) - ((int128)tmp_q[3] * 458449968195883L) + ((int128)tmp_q[4] * 1119059807217263L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 1119059807217263L) + ((int128)tmp_q[1] * 106375322864705L) - ((int128)tmp_q[2] * 1461880424685217L) - ((((int128)tmp_q[3] * 1486579628267562L) - ((int128)tmp_q[4] * 458449968195883L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 458449968195883L) + ((int128)tmp_q[1] * 1119059807217263L) + ((int128)tmp_q[2] * 106375322864705L) - ((int128)tmp_q[3] * 1461880424685217L) - ((int128)tmp_q[4] * 8919477769605372L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1486579628267562L) - ((int128)tmp_q[1] * 458449968195883L) + ((int128)tmp_q[2] * 1119059807217263L) + ((int128)tmp_q[3] * 106375322864705L) - ((int128)tmp_q[4] * 1461880424685217L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

