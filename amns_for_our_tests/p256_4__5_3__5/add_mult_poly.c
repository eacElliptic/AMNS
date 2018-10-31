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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8255914983287176538UL) + ((((uint64_t)op[1] * 16153039029075854701UL) + ((uint64_t)op[2] * 3577380165238963717UL) + ((uint64_t)op[3] * 13495628561828852981UL) + ((uint64_t)op[4] * 10935977291661429982UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 10935977291661429982UL) + ((uint64_t)op[1] * 8255914983287176538UL) + ((((uint64_t)op[2] * 16153039029075854701UL) + ((uint64_t)op[3] * 3577380165238963717UL) + ((uint64_t)op[4] * 13495628561828852981UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 13495628561828852981UL) + ((uint64_t)op[1] * 10935977291661429982UL) + ((uint64_t)op[2] * 8255914983287176538UL) + ((((uint64_t)op[3] * 16153039029075854701UL) + ((uint64_t)op[4] * 3577380165238963717UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 3577380165238963717UL) + ((uint64_t)op[1] * 13495628561828852981UL) + ((uint64_t)op[2] * 10935977291661429982UL) + ((uint64_t)op[3] * 8255914983287176538UL) + ((uint64_t)op[4] * 11565628939808460871UL);
	tmp_q[4] = ((uint64_t)op[0] * 16153039029075854701UL) + ((uint64_t)op[1] * 3577380165238963717UL) + ((uint64_t)op[2] * 13495628561828852981UL) + ((uint64_t)op[3] * 10935977291661429982UL) + ((uint64_t)op[4] * 8255914983287176538UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1389634421480411L) + ((((int128)tmp_q[1] * 1191763858762411L) + ((int128)tmp_q[2] * 11631372666706L) + ((int128)tmp_q[3] * 1428749020298217L) + ((int128)tmp_q[4] * 522524942087298L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 522524942087298L) + ((int128)tmp_q[1] * 1389634421480411L) + ((((int128)tmp_q[2] * 1191763858762411L) + ((int128)tmp_q[3] * 11631372666706L) + ((int128)tmp_q[4] * 1428749020298217L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 1428749020298217L) + ((int128)tmp_q[1] * 522524942087298L) + ((int128)tmp_q[2] * 1389634421480411L) + ((((int128)tmp_q[3] * 1191763858762411L) + ((int128)tmp_q[4] * 11631372666706L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 11631372666706L) + ((int128)tmp_q[1] * 1428749020298217L) + ((int128)tmp_q[2] * 522524942087298L) + ((int128)tmp_q[3] * 1389634421480411L) + ((int128)tmp_q[4] * 3575291576287233L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1191763858762411L) + ((int128)tmp_q[1] * 11631372666706L) + ((int128)tmp_q[2] * 1428749020298217L) + ((int128)tmp_q[3] * 522524942087298L) + ((int128)tmp_q[4] * 1389634421480411L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

