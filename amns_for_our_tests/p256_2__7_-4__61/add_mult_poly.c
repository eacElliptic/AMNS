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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9030525696376491695UL) + ((((uint64_t)op[1] * 8340024073186895658UL) + ((uint64_t)op[2] * 5955872375714650295UL) + ((uint64_t)op[3] * 14449361313192027176UL) + ((uint64_t)op[4] * 11577938930941387980UL) + ((uint64_t)op[5] * 5298209904972153097UL) + ((uint64_t)op[6] * 13301101424855230261UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 13301101424855230261UL) + ((uint64_t)op[1] * 9030525696376491695UL) + ((((uint64_t)op[2] * 8340024073186895658UL) + ((uint64_t)op[3] * 5955872375714650295UL) + ((uint64_t)op[4] * 14449361313192027176UL) + ((uint64_t)op[5] * 11577938930941387980UL) + ((uint64_t)op[6] * 5298209904972153097UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 5298209904972153097UL) + ((uint64_t)op[1] * 13301101424855230261UL) + ((uint64_t)op[2] * 9030525696376491695UL) + ((((uint64_t)op[3] * 8340024073186895658UL) + ((uint64_t)op[4] * 5955872375714650295UL) + ((uint64_t)op[5] * 14449361313192027176UL) + ((uint64_t)op[6] * 11577938930941387980UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 11577938930941387980UL) + ((uint64_t)op[1] * 5298209904972153097UL) + ((uint64_t)op[2] * 13301101424855230261UL) + ((uint64_t)op[3] * 9030525696376491695UL) + ((((uint64_t)op[4] * 8340024073186895658UL) + ((uint64_t)op[5] * 5955872375714650295UL) + ((uint64_t)op[6] * 14449361313192027176UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 14449361313192027176UL) + ((uint64_t)op[1] * 11577938930941387980UL) + ((uint64_t)op[2] * 5298209904972153097UL) + ((uint64_t)op[3] * 13301101424855230261UL) + ((uint64_t)op[4] * 9030525696376491695UL) + ((((uint64_t)op[5] * 8340024073186895658UL) + ((uint64_t)op[6] * 5955872375714650295UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 5955872375714650295UL) + ((uint64_t)op[1] * 14449361313192027176UL) + ((uint64_t)op[2] * 11577938930941387980UL) + ((uint64_t)op[3] * 5298209904972153097UL) + ((uint64_t)op[4] * 13301101424855230261UL) + ((uint64_t)op[5] * 9030525696376491695UL) + ((uint64_t)op[6] * 3533391854671520600UL);
	tmp_q[6] = ((uint64_t)op[0] * 8340024073186895658UL) + ((uint64_t)op[1] * 5955872375714650295UL) + ((uint64_t)op[2] * 14449361313192027176UL) + ((uint64_t)op[3] * 11577938930941387980UL) + ((uint64_t)op[4] * 5298209904972153097UL) + ((uint64_t)op[5] * 13301101424855230261UL) + ((uint64_t)op[6] * 9030525696376491695UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 70355632673L) - ((((int128)tmp_q[1] * 39154028245L) - ((int128)tmp_q[2] * 19843409821L) - ((int128)tmp_q[3] * 24237525435L) + ((int128)tmp_q[4] * 31793952531L) - ((int128)tmp_q[5] * 26998175242L) - ((int128)tmp_q[6] * 14984060579L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 14984060579L) + ((int128)tmp_q[1] * 70355632673L) - ((((int128)tmp_q[2] * 39154028245L) - ((int128)tmp_q[3] * 19843409821L) - ((int128)tmp_q[4] * 24237525435L) + ((int128)tmp_q[5] * 31793952531L) - ((int128)tmp_q[6] * 26998175242L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 26998175242L) - ((int128)tmp_q[1] * 14984060579L) + ((int128)tmp_q[2] * 70355632673L) - ((((int128)tmp_q[3] * 39154028245L) - ((int128)tmp_q[4] * 19843409821L) - ((int128)tmp_q[5] * 24237525435L) + ((int128)tmp_q[6] * 31793952531L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 31793952531L) - ((int128)tmp_q[1] * 26998175242L) - ((int128)tmp_q[2] * 14984060579L) + ((int128)tmp_q[3] * 70355632673L) - ((((int128)tmp_q[4] * 39154028245L) - ((int128)tmp_q[5] * 19843409821L) - ((int128)tmp_q[6] * 24237525435L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 24237525435L) + ((int128)tmp_q[1] * 31793952531L) - ((int128)tmp_q[2] * 26998175242L) - ((int128)tmp_q[3] * 14984060579L) + ((int128)tmp_q[4] * 70355632673L) - ((((int128)tmp_q[5] * 39154028245L) - ((int128)tmp_q[6] * 19843409821L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 19843409821L) - ((int128)tmp_q[1] * 24237525435L) + ((int128)tmp_q[2] * 31793952531L) - ((int128)tmp_q[3] * 26998175242L) - ((int128)tmp_q[4] * 14984060579L) + ((int128)tmp_q[5] * 70355632673L) - ((int128)tmp_q[6] * 156616112980L);
	tmp_zero[6] = ((int128)tmp_q[0] * 39154028245L) - ((int128)tmp_q[1] * 19843409821L) - ((int128)tmp_q[2] * 24237525435L) + ((int128)tmp_q[3] * 31793952531L) - ((int128)tmp_q[4] * 26998175242L) - ((int128)tmp_q[5] * 14984060579L) + ((int128)tmp_q[6] * 70355632673L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

