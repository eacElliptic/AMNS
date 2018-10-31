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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8354077480609028001UL) + ((((uint64_t)op[1] * 6814122039269400776UL) + ((uint64_t)op[2] * 15994791919001558350UL) + ((uint64_t)op[3] * 13929899528598130454UL) + ((uint64_t)op[4] * 4450378556328841638UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 4450378556328841638UL) + ((uint64_t)op[1] * 8354077480609028001UL) + ((((uint64_t)op[2] * 6814122039269400776UL) + ((uint64_t)op[3] * 15994791919001558350UL) + ((uint64_t)op[4] * 13929899528598130454UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 13929899528598130454UL) + ((uint64_t)op[1] * 4450378556328841638UL) + ((uint64_t)op[2] * 8354077480609028001UL) + ((((uint64_t)op[3] * 6814122039269400776UL) + ((uint64_t)op[4] * 15994791919001558350UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 15994791919001558350UL) + ((uint64_t)op[1] * 13929899528598130454UL) + ((uint64_t)op[2] * 4450378556328841638UL) + ((uint64_t)op[3] * 8354077480609028001UL) + ((uint64_t)op[4] * 2822877951072099352UL);
	tmp_q[4] = ((uint64_t)op[0] * 6814122039269400776UL) + ((uint64_t)op[1] * 15994791919001558350UL) + ((uint64_t)op[2] * 13929899528598130454UL) + ((uint64_t)op[3] * 4450378556328841638UL) + ((uint64_t)op[4] * 8354077480609028001UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5337584234441L) - ((-((int128)tmp_q[1] * 6616729443428L) + ((int128)tmp_q[2] * 17240752967910L) + ((int128)tmp_q[3] * 595317519602L) - ((int128)tmp_q[4] * 9842373408046L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 9842373408046L) - ((int128)tmp_q[1] * 5337584234441L) - ((-((int128)tmp_q[2] * 6616729443428L) + ((int128)tmp_q[3] * 17240752967910L) + ((int128)tmp_q[4] * 595317519602L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 595317519602L) - ((int128)tmp_q[1] * 9842373408046L) - ((int128)tmp_q[2] * 5337584234441L) - ((-((int128)tmp_q[3] * 6616729443428L) + ((int128)tmp_q[4] * 17240752967910L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 17240752967910L) + ((int128)tmp_q[1] * 595317519602L) - ((int128)tmp_q[2] * 9842373408046L) - ((int128)tmp_q[3] * 5337584234441L) + ((int128)tmp_q[4] * 33083647217140L);
	tmp_zero[4] = -((int128)tmp_q[0] * 6616729443428L) + ((int128)tmp_q[1] * 17240752967910L) + ((int128)tmp_q[2] * 595317519602L) - ((int128)tmp_q[3] * 9842373408046L) - ((int128)tmp_q[4] * 5337584234441L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

