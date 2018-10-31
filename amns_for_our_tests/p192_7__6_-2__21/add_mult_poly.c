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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 950810701173480383UL) + ((((uint64_t)op[1] * 12026704580619975407UL) + ((uint64_t)op[2] * 10940822869836911902UL) + ((uint64_t)op[3] * 9175100114726909362UL) + ((uint64_t)op[4] * 5727366632989300505UL) + ((uint64_t)op[5] * 17336541244018628812UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 17336541244018628812UL) + ((uint64_t)op[1] * 950810701173480383UL) + ((((uint64_t)op[2] * 12026704580619975407UL) + ((uint64_t)op[3] * 10940822869836911902UL) + ((uint64_t)op[4] * 9175100114726909362UL) + ((uint64_t)op[5] * 5727366632989300505UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 5727366632989300505UL) + ((uint64_t)op[1] * 17336541244018628812UL) + ((uint64_t)op[2] * 950810701173480383UL) + ((((uint64_t)op[3] * 12026704580619975407UL) + ((uint64_t)op[4] * 10940822869836911902UL) + ((uint64_t)op[5] * 9175100114726909362UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 9175100114726909362UL) + ((uint64_t)op[1] * 5727366632989300505UL) + ((uint64_t)op[2] * 17336541244018628812UL) + ((uint64_t)op[3] * 950810701173480383UL) + ((((uint64_t)op[4] * 12026704580619975407UL) + ((uint64_t)op[5] * 10940822869836911902UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 10940822869836911902UL) + ((uint64_t)op[1] * 9175100114726909362UL) + ((uint64_t)op[2] * 5727366632989300505UL) + ((uint64_t)op[3] * 17336541244018628812UL) + ((uint64_t)op[4] * 950810701173480383UL) + ((uint64_t)op[5] * 12840078986179152418UL);
	tmp_q[5] = ((uint64_t)op[0] * 12026704580619975407UL) + ((uint64_t)op[1] * 10940822869836911902UL) + ((uint64_t)op[2] * 9175100114726909362UL) + ((uint64_t)op[3] * 5727366632989300505UL) + ((uint64_t)op[4] * 17336541244018628812UL) + ((uint64_t)op[5] * 950810701173480383UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1717258631L) - ((((int128)tmp_q[1] * 1055534399L) - ((int128)tmp_q[2] * 2118296529L) - ((int128)tmp_q[3] * 570831052L) - ((int128)tmp_q[4] * 687655145L) + ((int128)tmp_q[5] * 406240456L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 406240456L) + ((int128)tmp_q[1] * 1717258631L) - ((((int128)tmp_q[2] * 1055534399L) - ((int128)tmp_q[3] * 2118296529L) - ((int128)tmp_q[4] * 570831052L) - ((int128)tmp_q[5] * 687655145L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 687655145L) + ((int128)tmp_q[1] * 406240456L) + ((int128)tmp_q[2] * 1717258631L) - ((((int128)tmp_q[3] * 1055534399L) - ((int128)tmp_q[4] * 2118296529L) - ((int128)tmp_q[5] * 570831052L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 570831052L) - ((int128)tmp_q[1] * 687655145L) + ((int128)tmp_q[2] * 406240456L) + ((int128)tmp_q[3] * 1717258631L) - ((((int128)tmp_q[4] * 1055534399L) - ((int128)tmp_q[5] * 2118296529L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 2118296529L) - ((int128)tmp_q[1] * 570831052L) - ((int128)tmp_q[2] * 687655145L) + ((int128)tmp_q[3] * 406240456L) + ((int128)tmp_q[4] * 1717258631L) - ((int128)tmp_q[5] * 2111068798L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1055534399L) - ((int128)tmp_q[1] * 2118296529L) - ((int128)tmp_q[2] * 570831052L) - ((int128)tmp_q[3] * 687655145L) + ((int128)tmp_q[4] * 406240456L) + ((int128)tmp_q[5] * 1717258631L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

