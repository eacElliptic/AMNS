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
	tmp_q[0] = ((uint64_t)op[0] * 5880306835568429345UL) + ((((uint64_t)op[1] * 13335133088743606744UL) + ((uint64_t)op[2] * 6387098314028184975UL) + ((uint64_t)op[3] * 16698693018273869621UL) + ((uint64_t)op[4] * 16723948618958320586UL) + ((uint64_t)op[5] * 15636737779396279428UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 15636737779396279428UL) + ((uint64_t)op[1] * 5880306835568429345UL) + ((((uint64_t)op[2] * 13335133088743606744UL) + ((uint64_t)op[3] * 6387098314028184975UL) + ((uint64_t)op[4] * 16698693018273869621UL) + ((uint64_t)op[5] * 16723948618958320586UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 16723948618958320586UL) + ((uint64_t)op[1] * 15636737779396279428UL) + ((uint64_t)op[2] * 5880306835568429345UL) + ((((uint64_t)op[3] * 13335133088743606744UL) + ((uint64_t)op[4] * 6387098314028184975UL) + ((uint64_t)op[5] * 16698693018273869621UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 16698693018273869621UL) + ((uint64_t)op[1] * 16723948618958320586UL) + ((uint64_t)op[2] * 15636737779396279428UL) + ((uint64_t)op[3] * 5880306835568429345UL) + ((((uint64_t)op[4] * 13335133088743606744UL) + ((uint64_t)op[5] * 6387098314028184975UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 6387098314028184975UL) + ((uint64_t)op[1] * 16698693018273869621UL) + ((uint64_t)op[2] * 16723948618958320586UL) + ((uint64_t)op[3] * 15636737779396279428UL) + ((uint64_t)op[4] * 5880306835568429345UL) + ((uint64_t)op[5] * 11335433222589378872UL);
	tmp_q[5] = ((uint64_t)op[0] * 13335133088743606744UL) + ((uint64_t)op[1] * 6387098314028184975UL) + ((uint64_t)op[2] * 16698693018273869621UL) + ((uint64_t)op[3] * 16723948618958320586UL) + ((uint64_t)op[4] * 15636737779396279428UL) + ((uint64_t)op[5] * 5880306835568429345UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 503742840876L) + ((((int128)tmp_q[1] * 1971007731576L) - ((int128)tmp_q[2] * 427474411999L) - ((int128)tmp_q[3] * 677518238282L) + ((int128)tmp_q[4] * 2261408742193L) + ((int128)tmp_q[5] * 147788887149L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 147788887149L) - ((int128)tmp_q[1] * 503742840876L) + ((((int128)tmp_q[2] * 1971007731576L) - ((int128)tmp_q[3] * 427474411999L) - ((int128)tmp_q[4] * 677518238282L) + ((int128)tmp_q[5] * 2261408742193L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 2261408742193L) + ((int128)tmp_q[1] * 147788887149L) - ((int128)tmp_q[2] * 503742840876L) + ((((int128)tmp_q[3] * 1971007731576L) - ((int128)tmp_q[4] * 427474411999L) - ((int128)tmp_q[5] * 677518238282L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 677518238282L) + ((int128)tmp_q[1] * 2261408742193L) + ((int128)tmp_q[2] * 147788887149L) - ((int128)tmp_q[3] * 503742840876L) + ((((int128)tmp_q[4] * 1971007731576L) - ((int128)tmp_q[5] * 427474411999L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 427474411999L) - ((int128)tmp_q[1] * 677518238282L) + ((int128)tmp_q[2] * 2261408742193L) + ((int128)tmp_q[3] * 147788887149L) - ((int128)tmp_q[4] * 503742840876L) + ((int128)tmp_q[5] * 9855038657880L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1971007731576L) - ((int128)tmp_q[1] * 427474411999L) - ((int128)tmp_q[2] * 677518238282L) + ((int128)tmp_q[3] * 2261408742193L) + ((int128)tmp_q[4] * 147788887149L) - ((int128)tmp_q[5] * 503742840876L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

