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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6651969163765975942UL) + ((((uint64_t)op[1] * 3759683516736399879UL) + ((uint64_t)op[2] * 5539867077003836923UL) + ((uint64_t)op[3] * 8670760365779699297UL) + ((uint64_t)op[4] * 671701724075428873UL) + ((uint64_t)op[5] * 17839135669603212895UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 17839135669603212895UL) + ((uint64_t)op[1] * 6651969163765975942UL) + ((((uint64_t)op[2] * 3759683516736399879UL) + ((uint64_t)op[3] * 5539867077003836923UL) + ((uint64_t)op[4] * 8670760365779699297UL) + ((uint64_t)op[5] * 671701724075428873UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 671701724075428873UL) + ((uint64_t)op[1] * 17839135669603212895UL) + ((uint64_t)op[2] * 6651969163765975942UL) + ((((uint64_t)op[3] * 3759683516736399879UL) + ((uint64_t)op[4] * 5539867077003836923UL) + ((uint64_t)op[5] * 8670760365779699297UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 8670760365779699297UL) + ((uint64_t)op[1] * 671701724075428873UL) + ((uint64_t)op[2] * 17839135669603212895UL) + ((uint64_t)op[3] * 6651969163765975942UL) + ((((uint64_t)op[4] * 3759683516736399879UL) + ((uint64_t)op[5] * 5539867077003836923UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 5539867077003836923UL) + ((uint64_t)op[1] * 8670760365779699297UL) + ((uint64_t)op[2] * 671701724075428873UL) + ((uint64_t)op[3] * 17839135669603212895UL) + ((uint64_t)op[4] * 6651969163765975942UL) + ((uint64_t)op[5] * 18095070563737103837UL);
	tmp_q[5] = ((uint64_t)op[0] * 3759683516736399879UL) + ((uint64_t)op[1] * 5539867077003836923UL) + ((uint64_t)op[2] * 8670760365779699297UL) + ((uint64_t)op[3] * 671701724075428873UL) + ((uint64_t)op[4] * 17839135669603212895UL) + ((uint64_t)op[5] * 6651969163765975942UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3121696764918L) - ((((int128)tmp_q[1] * 2300954559955L) + ((int128)tmp_q[2] * 4102018223213L) - ((int128)tmp_q[3] * 2417636448027L) - ((int128)tmp_q[4] * 2181588610913L) - ((int128)tmp_q[5] * 459169466549L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 459169466549L) - ((int128)tmp_q[1] * 3121696764918L) - ((((int128)tmp_q[2] * 2300954559955L) + ((int128)tmp_q[3] * 4102018223213L) - ((int128)tmp_q[4] * 2417636448027L) - ((int128)tmp_q[5] * 2181588610913L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 2181588610913L) - ((int128)tmp_q[1] * 459169466549L) - ((int128)tmp_q[2] * 3121696764918L) - ((((int128)tmp_q[3] * 2300954559955L) + ((int128)tmp_q[4] * 4102018223213L) - ((int128)tmp_q[5] * 2417636448027L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2417636448027L) - ((int128)tmp_q[1] * 2181588610913L) - ((int128)tmp_q[2] * 459169466549L) - ((int128)tmp_q[3] * 3121696764918L) - ((((int128)tmp_q[4] * 2300954559955L) + ((int128)tmp_q[5] * 4102018223213L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 4102018223213L) - ((int128)tmp_q[1] * 2417636448027L) - ((int128)tmp_q[2] * 2181588610913L) - ((int128)tmp_q[3] * 459169466549L) - ((int128)tmp_q[4] * 3121696764918L) - ((int128)tmp_q[5] * 11504772799775L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2300954559955L) + ((int128)tmp_q[1] * 4102018223213L) - ((int128)tmp_q[2] * 2417636448027L) - ((int128)tmp_q[3] * 2181588610913L) - ((int128)tmp_q[4] * 459169466549L) - ((int128)tmp_q[5] * 3121696764918L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

