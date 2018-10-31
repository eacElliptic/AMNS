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
	tmp_q[0] = ((uint64_t)op[0] * 11715468295012194645UL) + ((((uint64_t)op[1] * 12345707431396833080UL) + ((uint64_t)op[2] * 5676202577114121442UL) + ((uint64_t)op[3] * 215971771419178787UL) + ((uint64_t)op[4] * 13146925077789882659UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 13146925077789882659UL) + ((uint64_t)op[1] * 11715468295012194645UL) + ((((uint64_t)op[2] * 12345707431396833080UL) + ((uint64_t)op[3] * 5676202577114121442UL) + ((uint64_t)op[4] * 215971771419178787UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 215971771419178787UL) + ((uint64_t)op[1] * 13146925077789882659UL) + ((uint64_t)op[2] * 11715468295012194645UL) + ((((uint64_t)op[3] * 12345707431396833080UL) + ((uint64_t)op[4] * 5676202577114121442UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 5676202577114121442UL) + ((uint64_t)op[1] * 215971771419178787UL) + ((uint64_t)op[2] * 13146925077789882659UL) + ((uint64_t)op[3] * 11715468295012194645UL) + ((uint64_t)op[4] * 12058439137854041064UL);
	tmp_q[4] = ((uint64_t)op[0] * 12345707431396833080UL) + ((uint64_t)op[1] * 5676202577114121442UL) + ((uint64_t)op[2] * 215971771419178787UL) + ((uint64_t)op[3] * 13146925077789882659UL) + ((uint64_t)op[4] * 11715468295012194645UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 292533345974L) - ((((int128)tmp_q[1] * 227826381805L) + ((int128)tmp_q[2] * 42981250982L) + ((int128)tmp_q[3] * 232908254553L) + ((int128)tmp_q[4] * 231045908321L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 231045908321L) + ((int128)tmp_q[1] * 292533345974L) - ((((int128)tmp_q[2] * 227826381805L) + ((int128)tmp_q[3] * 42981250982L) + ((int128)tmp_q[4] * 232908254553L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 232908254553L) + ((int128)tmp_q[1] * 231045908321L) + ((int128)tmp_q[2] * 292533345974L) - ((((int128)tmp_q[3] * 227826381805L) + ((int128)tmp_q[4] * 42981250982L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 42981250982L) + ((int128)tmp_q[1] * 232908254553L) + ((int128)tmp_q[2] * 231045908321L) + ((int128)tmp_q[3] * 292533345974L) - ((int128)tmp_q[4] * 1139131909025L);
	tmp_zero[4] = ((int128)tmp_q[0] * 227826381805L) + ((int128)tmp_q[1] * 42981250982L) + ((int128)tmp_q[2] * 232908254553L) + ((int128)tmp_q[3] * 231045908321L) + ((int128)tmp_q[4] * 292533345974L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

