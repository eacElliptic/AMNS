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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7694687945627488685UL) + ((((uint64_t)op[1] * 15036470861897899653UL) + ((uint64_t)op[2] * 2982981179193398117UL) + ((uint64_t)op[3] * 10873980292445264822UL) + ((uint64_t)op[4] * 8703530102385862299UL) + ((uint64_t)op[5] * 6226472646285889159UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 6226472646285889159UL) + ((uint64_t)op[1] * 7694687945627488685UL) + ((((uint64_t)op[2] * 15036470861897899653UL) + ((uint64_t)op[3] * 2982981179193398117UL) + ((uint64_t)op[4] * 10873980292445264822UL) + ((uint64_t)op[5] * 8703530102385862299UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 8703530102385862299UL) + ((uint64_t)op[1] * 6226472646285889159UL) + ((uint64_t)op[2] * 7694687945627488685UL) + ((((uint64_t)op[3] * 15036470861897899653UL) + ((uint64_t)op[4] * 2982981179193398117UL) + ((uint64_t)op[5] * 10873980292445264822UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 10873980292445264822UL) + ((uint64_t)op[1] * 8703530102385862299UL) + ((uint64_t)op[2] * 6226472646285889159UL) + ((uint64_t)op[3] * 7694687945627488685UL) + ((((uint64_t)op[4] * 15036470861897899653UL) + ((uint64_t)op[5] * 2982981179193398117UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 2982981179193398117UL) + ((uint64_t)op[1] * 10873980292445264822UL) + ((uint64_t)op[2] * 8703530102385862299UL) + ((uint64_t)op[3] * 6226472646285889159UL) + ((uint64_t)op[4] * 7694687945627488685UL) + ((uint64_t)op[5] * 10230819635434955889UL);
	tmp_q[5] = ((uint64_t)op[0] * 15036470861897899653UL) + ((uint64_t)op[1] * 2982981179193398117UL) + ((uint64_t)op[2] * 10873980292445264822UL) + ((uint64_t)op[3] * 8703530102385862299UL) + ((uint64_t)op[4] * 6226472646285889159UL) + ((uint64_t)op[5] * 7694687945627488685UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 905234880197L) - ((((int128)tmp_q[1] * 1508276268283L) + ((int128)tmp_q[2] * 2973307427515L) + ((int128)tmp_q[3] * 1232196310098L) + ((int128)tmp_q[4] * 2599814355593L) + ((int128)tmp_q[5] * 2531701835325L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 2531701835325L) - ((int128)tmp_q[1] * 905234880197L) - ((((int128)tmp_q[2] * 1508276268283L) + ((int128)tmp_q[3] * 2973307427515L) + ((int128)tmp_q[4] * 1232196310098L) + ((int128)tmp_q[5] * 2599814355593L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2599814355593L) + ((int128)tmp_q[1] * 2531701835325L) - ((int128)tmp_q[2] * 905234880197L) - ((((int128)tmp_q[3] * 1508276268283L) + ((int128)tmp_q[4] * 2973307427515L) + ((int128)tmp_q[5] * 1232196310098L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1232196310098L) + ((int128)tmp_q[1] * 2599814355593L) + ((int128)tmp_q[2] * 2531701835325L) - ((int128)tmp_q[3] * 905234880197L) - ((((int128)tmp_q[4] * 1508276268283L) + ((int128)tmp_q[5] * 2973307427515L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 2973307427515L) + ((int128)tmp_q[1] * 1232196310098L) + ((int128)tmp_q[2] * 2599814355593L) + ((int128)tmp_q[3] * 2531701835325L) - ((int128)tmp_q[4] * 905234880197L) - ((int128)tmp_q[5] * 4524828804849L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1508276268283L) + ((int128)tmp_q[1] * 2973307427515L) + ((int128)tmp_q[2] * 1232196310098L) + ((int128)tmp_q[3] * 2599814355593L) + ((int128)tmp_q[4] * 2531701835325L) - ((int128)tmp_q[5] * 905234880197L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

