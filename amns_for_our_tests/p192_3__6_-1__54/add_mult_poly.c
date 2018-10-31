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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12732819988388589337UL) + ((((uint64_t)op[1] * 2105514175937775181UL) + ((uint64_t)op[2] * 3560746091600662547UL) + ((uint64_t)op[3] * 13432326373923409027UL) + ((uint64_t)op[4] * 4941723573091066274UL) + ((uint64_t)op[5] * 9907831884623360645UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 9907831884623360645UL) + ((uint64_t)op[1] * 12732819988388589337UL) + ((((uint64_t)op[2] * 2105514175937775181UL) + ((uint64_t)op[3] * 3560746091600662547UL) + ((uint64_t)op[4] * 13432326373923409027UL) + ((uint64_t)op[5] * 4941723573091066274UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 4941723573091066274UL) + ((uint64_t)op[1] * 9907831884623360645UL) + ((uint64_t)op[2] * 12732819988388589337UL) + ((((uint64_t)op[3] * 2105514175937775181UL) + ((uint64_t)op[4] * 3560746091600662547UL) + ((uint64_t)op[5] * 13432326373923409027UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 13432326373923409027UL) + ((uint64_t)op[1] * 4941723573091066274UL) + ((uint64_t)op[2] * 9907831884623360645UL) + ((uint64_t)op[3] * 12732819988388589337UL) + ((((uint64_t)op[4] * 2105514175937775181UL) + ((uint64_t)op[5] * 3560746091600662547UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 3560746091600662547UL) + ((uint64_t)op[1] * 13432326373923409027UL) + ((uint64_t)op[2] * 4941723573091066274UL) + ((uint64_t)op[3] * 9907831884623360645UL) + ((uint64_t)op[4] * 12732819988388589337UL) + ((uint64_t)op[5] * 16341229897771776435UL);
	tmp_q[5] = ((uint64_t)op[0] * 2105514175937775181UL) + ((uint64_t)op[1] * 3560746091600662547UL) + ((uint64_t)op[2] * 13432326373923409027UL) + ((uint64_t)op[3] * 4941723573091066274UL) + ((uint64_t)op[4] * 9907831884623360645UL) + ((uint64_t)op[5] * 12732819988388589337UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 25198163024455L) - (-((int128)tmp_q[1] * 34407772882727L) - ((int128)tmp_q[2] * 137235661266854L) + ((int128)tmp_q[3] * 120501804324415L) - ((int128)tmp_q[4] * 162433824291307L) + ((int128)tmp_q[5] * 154909577207145L));
	tmp_zero[1] = ((int128)tmp_q[0] * 154909577207145L) - ((int128)tmp_q[1] * 25198163024455L) - (-((int128)tmp_q[2] * 34407772882727L) - ((int128)tmp_q[3] * 137235661266854L) + ((int128)tmp_q[4] * 120501804324415L) - ((int128)tmp_q[5] * 162433824291307L));
	tmp_zero[2] = -((int128)tmp_q[0] * 162433824291307L) + ((int128)tmp_q[1] * 154909577207145L) - ((int128)tmp_q[2] * 25198163024455L) - (-((int128)tmp_q[3] * 34407772882727L) - ((int128)tmp_q[4] * 137235661266854L) + ((int128)tmp_q[5] * 120501804324415L));
	tmp_zero[3] = ((int128)tmp_q[0] * 120501804324415L) - ((int128)tmp_q[1] * 162433824291307L) + ((int128)tmp_q[2] * 154909577207145L) - ((int128)tmp_q[3] * 25198163024455L) - (-((int128)tmp_q[4] * 34407772882727L) - ((int128)tmp_q[5] * 137235661266854L));
	tmp_zero[4] = -((int128)tmp_q[0] * 137235661266854L) + ((int128)tmp_q[1] * 120501804324415L) - ((int128)tmp_q[2] * 162433824291307L) + ((int128)tmp_q[3] * 154909577207145L) - ((int128)tmp_q[4] * 25198163024455L) + ((int128)tmp_q[5] * 34407772882727L);
	tmp_zero[5] = -((int128)tmp_q[0] * 34407772882727L) - ((int128)tmp_q[1] * 137235661266854L) + ((int128)tmp_q[2] * 120501804324415L) - ((int128)tmp_q[3] * 162433824291307L) + ((int128)tmp_q[4] * 154909577207145L) - ((int128)tmp_q[5] * 25198163024455L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

