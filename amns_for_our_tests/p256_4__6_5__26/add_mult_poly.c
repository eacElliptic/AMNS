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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 649004168493624797UL) + ((((uint64_t)op[1] * 10498333618907959607UL) + ((uint64_t)op[2] * 15033928591955064065UL) + ((uint64_t)op[3] * 11538416242778337869UL) + ((uint64_t)op[4] * 7353446789147583596UL) + ((uint64_t)op[5] * 2642269354282590125UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 2642269354282590125UL) + ((uint64_t)op[1] * 649004168493624797UL) + ((((uint64_t)op[2] * 10498333618907959607UL) + ((uint64_t)op[3] * 15033928591955064065UL) + ((uint64_t)op[4] * 11538416242778337869UL) + ((uint64_t)op[5] * 7353446789147583596UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 7353446789147583596UL) + ((uint64_t)op[1] * 2642269354282590125UL) + ((uint64_t)op[2] * 649004168493624797UL) + ((((uint64_t)op[3] * 10498333618907959607UL) + ((uint64_t)op[4] * 15033928591955064065UL) + ((uint64_t)op[5] * 11538416242778337869UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 11538416242778337869UL) + ((uint64_t)op[1] * 7353446789147583596UL) + ((uint64_t)op[2] * 2642269354282590125UL) + ((uint64_t)op[3] * 649004168493624797UL) + ((((uint64_t)op[4] * 10498333618907959607UL) + ((uint64_t)op[5] * 15033928591955064065UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 15033928591955064065UL) + ((uint64_t)op[1] * 11538416242778337869UL) + ((uint64_t)op[2] * 7353446789147583596UL) + ((uint64_t)op[3] * 2642269354282590125UL) + ((uint64_t)op[4] * 649004168493624797UL) + ((uint64_t)op[5] * 15598179947120694803UL);
	tmp_q[5] = ((uint64_t)op[0] * 10498333618907959607UL) + ((uint64_t)op[1] * 15033928591955064065UL) + ((uint64_t)op[2] * 11538416242778337869UL) + ((uint64_t)op[3] * 7353446789147583596UL) + ((uint64_t)op[4] * 2642269354282590125UL) + ((uint64_t)op[5] * 649004168493624797UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1949045595289L) + ((((int128)tmp_q[1] * 6161514114523L) + ((int128)tmp_q[2] * 4265703394496L) + ((int128)tmp_q[3] * 1890440451891L) + ((int128)tmp_q[4] * 4573441916695L) + ((int128)tmp_q[5] * 983687878817L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 983687878817L) - ((int128)tmp_q[1] * 1949045595289L) + ((((int128)tmp_q[2] * 6161514114523L) + ((int128)tmp_q[3] * 4265703394496L) + ((int128)tmp_q[4] * 1890440451891L) + ((int128)tmp_q[5] * 4573441916695L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 4573441916695L) + ((int128)tmp_q[1] * 983687878817L) - ((int128)tmp_q[2] * 1949045595289L) + ((((int128)tmp_q[3] * 6161514114523L) + ((int128)tmp_q[4] * 4265703394496L) + ((int128)tmp_q[5] * 1890440451891L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1890440451891L) + ((int128)tmp_q[1] * 4573441916695L) + ((int128)tmp_q[2] * 983687878817L) - ((int128)tmp_q[3] * 1949045595289L) + ((((int128)tmp_q[4] * 6161514114523L) + ((int128)tmp_q[5] * 4265703394496L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 4265703394496L) + ((int128)tmp_q[1] * 1890440451891L) + ((int128)tmp_q[2] * 4573441916695L) + ((int128)tmp_q[3] * 983687878817L) - ((int128)tmp_q[4] * 1949045595289L) + ((int128)tmp_q[5] * 30807570572615L);
	tmp_zero[5] = ((int128)tmp_q[0] * 6161514114523L) + ((int128)tmp_q[1] * 4265703394496L) + ((int128)tmp_q[2] * 1890440451891L) + ((int128)tmp_q[3] * 4573441916695L) + ((int128)tmp_q[4] * 983687878817L) - ((int128)tmp_q[5] * 1949045595289L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

