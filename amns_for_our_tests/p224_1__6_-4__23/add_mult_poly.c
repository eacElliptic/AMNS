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
	tmp_q[0] = ((uint64_t)op[0] * 1295451224129235079UL) + ((((uint64_t)op[1] * 345182176405695037UL) + ((uint64_t)op[2] * 15883458953462542526UL) + ((uint64_t)op[3] * 1854621979311610394UL) + ((uint64_t)op[4] * 44498175930170337UL) + ((uint64_t)op[5] * 7699090664147719508UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 7699090664147719508UL) + ((uint64_t)op[1] * 1295451224129235079UL) + ((((uint64_t)op[2] * 345182176405695037UL) + ((uint64_t)op[3] * 15883458953462542526UL) + ((uint64_t)op[4] * 1854621979311610394UL) + ((uint64_t)op[5] * 44498175930170337UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 44498175930170337UL) + ((uint64_t)op[1] * 7699090664147719508UL) + ((uint64_t)op[2] * 1295451224129235079UL) + ((((uint64_t)op[3] * 345182176405695037UL) + ((uint64_t)op[4] * 15883458953462542526UL) + ((uint64_t)op[5] * 1854621979311610394UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 1854621979311610394UL) + ((uint64_t)op[1] * 44498175930170337UL) + ((uint64_t)op[2] * 7699090664147719508UL) + ((uint64_t)op[3] * 1295451224129235079UL) + ((((uint64_t)op[4] * 345182176405695037UL) + ((uint64_t)op[5] * 15883458953462542526UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 15883458953462542526UL) + ((uint64_t)op[1] * 1854621979311610394UL) + ((uint64_t)op[2] * 44498175930170337UL) + ((uint64_t)op[3] * 7699090664147719508UL) + ((uint64_t)op[4] * 1295451224129235079UL) + ((uint64_t)op[5] * 17066015368086771468UL);
	tmp_q[5] = ((uint64_t)op[0] * 345182176405695037UL) + ((uint64_t)op[1] * 15883458953462542526UL) + ((uint64_t)op[2] * 1854621979311610394UL) + ((uint64_t)op[3] * 44498175930170337UL) + ((uint64_t)op[4] * 7699090664147719508UL) + ((uint64_t)op[5] * 1295451224129235079UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 497419131L) - ((-((int128)tmp_q[1] * 32648288875L) - ((int128)tmp_q[2] * 88918727937L) - ((int128)tmp_q[3] * 71617235658L) + ((int128)tmp_q[4] * 35141716309L) - ((int128)tmp_q[5] * 11361638908L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 11361638908L) - ((int128)tmp_q[1] * 497419131L) - ((-((int128)tmp_q[2] * 32648288875L) - ((int128)tmp_q[3] * 88918727937L) - ((int128)tmp_q[4] * 71617235658L) + ((int128)tmp_q[5] * 35141716309L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 35141716309L) - ((int128)tmp_q[1] * 11361638908L) - ((int128)tmp_q[2] * 497419131L) - ((-((int128)tmp_q[3] * 32648288875L) - ((int128)tmp_q[4] * 88918727937L) - ((int128)tmp_q[5] * 71617235658L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 71617235658L) + ((int128)tmp_q[1] * 35141716309L) - ((int128)tmp_q[2] * 11361638908L) - ((int128)tmp_q[3] * 497419131L) - ((-((int128)tmp_q[4] * 32648288875L) - ((int128)tmp_q[5] * 88918727937L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 88918727937L) - ((int128)tmp_q[1] * 71617235658L) + ((int128)tmp_q[2] * 35141716309L) - ((int128)tmp_q[3] * 11361638908L) - ((int128)tmp_q[4] * 497419131L) + ((int128)tmp_q[5] * 130593155500L);
	tmp_zero[5] = -((int128)tmp_q[0] * 32648288875L) - ((int128)tmp_q[1] * 88918727937L) - ((int128)tmp_q[2] * 71617235658L) + ((int128)tmp_q[3] * 35141716309L) - ((int128)tmp_q[4] * 11361638908L) - ((int128)tmp_q[5] * 497419131L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

