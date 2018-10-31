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
	tmp_q[0] = ((uint64_t)op[0] * 15611424517529530459UL) + ((((uint64_t)op[1] * 672908585071280556UL) + ((uint64_t)op[2] * 6613079841691256793UL) + ((uint64_t)op[3] * 15263594000514281698UL) + ((uint64_t)op[4] * 10786031619697069354UL) + ((uint64_t)op[5] * 7985775036176989518UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 7985775036176989518UL) + ((uint64_t)op[1] * 15611424517529530459UL) + ((((uint64_t)op[2] * 672908585071280556UL) + ((uint64_t)op[3] * 6613079841691256793UL) + ((uint64_t)op[4] * 15263594000514281698UL) + ((uint64_t)op[5] * 10786031619697069354UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 10786031619697069354UL) + ((uint64_t)op[1] * 7985775036176989518UL) + ((uint64_t)op[2] * 15611424517529530459UL) + ((((uint64_t)op[3] * 672908585071280556UL) + ((uint64_t)op[4] * 6613079841691256793UL) + ((uint64_t)op[5] * 15263594000514281698UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 15263594000514281698UL) + ((uint64_t)op[1] * 10786031619697069354UL) + ((uint64_t)op[2] * 7985775036176989518UL) + ((uint64_t)op[3] * 15611424517529530459UL) + ((((uint64_t)op[4] * 672908585071280556UL) + ((uint64_t)op[5] * 6613079841691256793UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 6613079841691256793UL) + ((uint64_t)op[1] * 15263594000514281698UL) + ((uint64_t)op[2] * 10786031619697069354UL) + ((uint64_t)op[3] * 7985775036176989518UL) + ((uint64_t)op[4] * 15611424517529530459UL) + ((uint64_t)op[5] * 15755109733424429392UL);
	tmp_q[5] = ((uint64_t)op[0] * 672908585071280556UL) + ((uint64_t)op[1] * 6613079841691256793UL) + ((uint64_t)op[2] * 15263594000514281698UL) + ((uint64_t)op[3] * 10786031619697069354UL) + ((uint64_t)op[4] * 7985775036176989518UL) + ((uint64_t)op[5] * 15611424517529530459UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 28809354867L) - ((-((int128)tmp_q[1] * 5714611288L) + ((int128)tmp_q[2] * 24274952157L) + ((int128)tmp_q[3] * 103596607002L) - ((int128)tmp_q[4] * 7537175798L) + ((int128)tmp_q[5] * 101850777678L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 101850777678L) - ((int128)tmp_q[1] * 28809354867L) - ((-((int128)tmp_q[2] * 5714611288L) + ((int128)tmp_q[3] * 24274952157L) + ((int128)tmp_q[4] * 103596607002L) - ((int128)tmp_q[5] * 7537175798L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 7537175798L) + ((int128)tmp_q[1] * 101850777678L) - ((int128)tmp_q[2] * 28809354867L) - ((-((int128)tmp_q[3] * 5714611288L) + ((int128)tmp_q[4] * 24274952157L) + ((int128)tmp_q[5] * 103596607002L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 103596607002L) - ((int128)tmp_q[1] * 7537175798L) + ((int128)tmp_q[2] * 101850777678L) - ((int128)tmp_q[3] * 28809354867L) - ((-((int128)tmp_q[4] * 5714611288L) + ((int128)tmp_q[5] * 24274952157L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 24274952157L) + ((int128)tmp_q[1] * 103596607002L) - ((int128)tmp_q[2] * 7537175798L) + ((int128)tmp_q[3] * 101850777678L) - ((int128)tmp_q[4] * 28809354867L) + ((int128)tmp_q[5] * 22858445152L);
	tmp_zero[5] = -((int128)tmp_q[0] * 5714611288L) + ((int128)tmp_q[1] * 24274952157L) + ((int128)tmp_q[2] * 103596607002L) - ((int128)tmp_q[3] * 7537175798L) + ((int128)tmp_q[4] * 101850777678L) - ((int128)tmp_q[5] * 28809354867L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

