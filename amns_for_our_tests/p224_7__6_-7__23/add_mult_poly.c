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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10521125701040736494UL) + ((((uint64_t)op[1] * 11280804208994989261UL) + ((uint64_t)op[2] * 13963967372579492751UL) + ((uint64_t)op[3] * 9949377925830247772UL) + ((uint64_t)op[4] * 17061016913264105423UL) + ((uint64_t)op[5] * 16709127374315262144UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 16709127374315262144UL) + ((uint64_t)op[1] * 10521125701040736494UL) + ((((uint64_t)op[2] * 11280804208994989261UL) + ((uint64_t)op[3] * 13963967372579492751UL) + ((uint64_t)op[4] * 9949377925830247772UL) + ((uint64_t)op[5] * 17061016913264105423UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 17061016913264105423UL) + ((uint64_t)op[1] * 16709127374315262144UL) + ((uint64_t)op[2] * 10521125701040736494UL) + ((((uint64_t)op[3] * 11280804208994989261UL) + ((uint64_t)op[4] * 13963967372579492751UL) + ((uint64_t)op[5] * 9949377925830247772UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 9949377925830247772UL) + ((uint64_t)op[1] * 17061016913264105423UL) + ((uint64_t)op[2] * 16709127374315262144UL) + ((uint64_t)op[3] * 10521125701040736494UL) + ((((uint64_t)op[4] * 11280804208994989261UL) + ((uint64_t)op[5] * 13963967372579492751UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 13963967372579492751UL) + ((uint64_t)op[1] * 9949377925830247772UL) + ((uint64_t)op[2] * 17061016913264105423UL) + ((uint64_t)op[3] * 16709127374315262144UL) + ((uint64_t)op[4] * 10521125701040736494UL) + ((uint64_t)op[5] * 13268090905582833253UL);
	tmp_q[5] = ((uint64_t)op[0] * 11280804208994989261UL) + ((uint64_t)op[1] * 13963967372579492751UL) + ((uint64_t)op[2] * 9949377925830247772UL) + ((uint64_t)op[3] * 17061016913264105423UL) + ((uint64_t)op[4] * 16709127374315262144UL) + ((uint64_t)op[5] * 10521125701040736494UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 102436406321L) - ((((int128)tmp_q[1] * 66821299752L) - ((int128)tmp_q[2] * 99120499324L) + ((int128)tmp_q[3] * 49503898625L) - ((int128)tmp_q[4] * 19359518983L) - ((int128)tmp_q[5] * 39968927202L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 39968927202L) - ((int128)tmp_q[1] * 102436406321L) - ((((int128)tmp_q[2] * 66821299752L) - ((int128)tmp_q[3] * 99120499324L) + ((int128)tmp_q[4] * 49503898625L) - ((int128)tmp_q[5] * 19359518983L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 19359518983L) - ((int128)tmp_q[1] * 39968927202L) - ((int128)tmp_q[2] * 102436406321L) - ((((int128)tmp_q[3] * 66821299752L) - ((int128)tmp_q[4] * 99120499324L) + ((int128)tmp_q[5] * 49503898625L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 49503898625L) - ((int128)tmp_q[1] * 19359518983L) - ((int128)tmp_q[2] * 39968927202L) - ((int128)tmp_q[3] * 102436406321L) - ((((int128)tmp_q[4] * 66821299752L) - ((int128)tmp_q[5] * 99120499324L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 99120499324L) + ((int128)tmp_q[1] * 49503898625L) - ((int128)tmp_q[2] * 19359518983L) - ((int128)tmp_q[3] * 39968927202L) - ((int128)tmp_q[4] * 102436406321L) - ((int128)tmp_q[5] * 467749098264L);
	tmp_zero[5] = ((int128)tmp_q[0] * 66821299752L) - ((int128)tmp_q[1] * 99120499324L) + ((int128)tmp_q[2] * 49503898625L) - ((int128)tmp_q[3] * 19359518983L) - ((int128)tmp_q[4] * 39968927202L) - ((int128)tmp_q[5] * 102436406321L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

