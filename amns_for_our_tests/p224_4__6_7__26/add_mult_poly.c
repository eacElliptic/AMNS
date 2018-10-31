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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13528445699159515380UL) + ((((uint64_t)op[1] * 16797350767652505712UL) + ((uint64_t)op[2] * 7190411189035765UL) + ((uint64_t)op[3] * 16935525157327162425UL) + ((uint64_t)op[4] * 3365919721454478312UL) + ((uint64_t)op[5] * 2387172697875651221UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 2387172697875651221UL) + ((uint64_t)op[1] * 13528445699159515380UL) + ((((uint64_t)op[2] * 16797350767652505712UL) + ((uint64_t)op[3] * 7190411189035765UL) + ((uint64_t)op[4] * 16935525157327162425UL) + ((uint64_t)op[5] * 3365919721454478312UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 3365919721454478312UL) + ((uint64_t)op[1] * 2387172697875651221UL) + ((uint64_t)op[2] * 13528445699159515380UL) + ((((uint64_t)op[3] * 16797350767652505712UL) + ((uint64_t)op[4] * 7190411189035765UL) + ((uint64_t)op[5] * 16935525157327162425UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 16935525157327162425UL) + ((uint64_t)op[1] * 3365919721454478312UL) + ((uint64_t)op[2] * 2387172697875651221UL) + ((uint64_t)op[3] * 13528445699159515380UL) + ((((uint64_t)op[4] * 16797350767652505712UL) + ((uint64_t)op[5] * 7190411189035765UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 7190411189035765UL) + ((uint64_t)op[1] * 16935525157327162425UL) + ((uint64_t)op[2] * 3365919721454478312UL) + ((uint64_t)op[3] * 2387172697875651221UL) + ((uint64_t)op[4] * 13528445699159515380UL) + ((uint64_t)op[5] * 6900990931310230288UL);
	tmp_q[5] = ((uint64_t)op[0] * 16797350767652505712UL) + ((uint64_t)op[1] * 7190411189035765UL) + ((uint64_t)op[2] * 16935525157327162425UL) + ((uint64_t)op[3] * 3365919721454478312UL) + ((uint64_t)op[4] * 2387172697875651221UL) + ((uint64_t)op[5] * 13528445699159515380UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 47536075906L) + ((-((int128)tmp_q[1] * 73524418268L) - ((int128)tmp_q[2] * 96862748301L) - ((int128)tmp_q[3] * 71586619937L) + ((int128)tmp_q[4] * 42974831254L) + ((int128)tmp_q[5] * 63422537463L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 63422537463L) + ((int128)tmp_q[1] * 47536075906L) + ((-((int128)tmp_q[2] * 73524418268L) - ((int128)tmp_q[3] * 96862748301L) - ((int128)tmp_q[4] * 71586619937L) + ((int128)tmp_q[5] * 42974831254L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 42974831254L) + ((int128)tmp_q[1] * 63422537463L) + ((int128)tmp_q[2] * 47536075906L) + ((-((int128)tmp_q[3] * 73524418268L) - ((int128)tmp_q[4] * 96862748301L) - ((int128)tmp_q[5] * 71586619937L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 71586619937L) + ((int128)tmp_q[1] * 42974831254L) + ((int128)tmp_q[2] * 63422537463L) + ((int128)tmp_q[3] * 47536075906L) + ((-((int128)tmp_q[4] * 73524418268L) - ((int128)tmp_q[5] * 96862748301L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 96862748301L) - ((int128)tmp_q[1] * 71586619937L) + ((int128)tmp_q[2] * 42974831254L) + ((int128)tmp_q[3] * 63422537463L) + ((int128)tmp_q[4] * 47536075906L) - ((int128)tmp_q[5] * 514670927876L);
	tmp_zero[5] = -((int128)tmp_q[0] * 73524418268L) - ((int128)tmp_q[1] * 96862748301L) - ((int128)tmp_q[2] * 71586619937L) + ((int128)tmp_q[3] * 42974831254L) + ((int128)tmp_q[4] * 63422537463L) + ((int128)tmp_q[5] * 47536075906L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

