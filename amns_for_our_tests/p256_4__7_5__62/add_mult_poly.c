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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6203048794235726393UL) + ((((uint64_t)op[1] * 10493488509697073148UL) + ((uint64_t)op[2] * 10458323379115823329UL) + ((uint64_t)op[3] * 17924663376254293634UL) + ((uint64_t)op[4] * 14401118761468768990UL) + ((uint64_t)op[5] * 7865419579248818977UL) + ((uint64_t)op[6] * 6024154042657686016UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 6024154042657686016UL) + ((uint64_t)op[1] * 6203048794235726393UL) + ((((uint64_t)op[2] * 10493488509697073148UL) + ((uint64_t)op[3] * 10458323379115823329UL) + ((uint64_t)op[4] * 17924663376254293634UL) + ((uint64_t)op[5] * 14401118761468768990UL) + ((uint64_t)op[6] * 7865419579248818977UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 7865419579248818977UL) + ((uint64_t)op[1] * 6024154042657686016UL) + ((uint64_t)op[2] * 6203048794235726393UL) + ((((uint64_t)op[3] * 10493488509697073148UL) + ((uint64_t)op[4] * 10458323379115823329UL) + ((uint64_t)op[5] * 17924663376254293634UL) + ((uint64_t)op[6] * 14401118761468768990UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 14401118761468768990UL) + ((uint64_t)op[1] * 7865419579248818977UL) + ((uint64_t)op[2] * 6024154042657686016UL) + ((uint64_t)op[3] * 6203048794235726393UL) + ((((uint64_t)op[4] * 10493488509697073148UL) + ((uint64_t)op[5] * 10458323379115823329UL) + ((uint64_t)op[6] * 17924663376254293634UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 17924663376254293634UL) + ((uint64_t)op[1] * 14401118761468768990UL) + ((uint64_t)op[2] * 7865419579248818977UL) + ((uint64_t)op[3] * 6024154042657686016UL) + ((uint64_t)op[4] * 6203048794235726393UL) + ((((uint64_t)op[5] * 10493488509697073148UL) + ((uint64_t)op[6] * 10458323379115823329UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 10458323379115823329UL) + ((uint64_t)op[1] * 17924663376254293634UL) + ((uint64_t)op[2] * 14401118761468768990UL) + ((uint64_t)op[3] * 7865419579248818977UL) + ((uint64_t)op[4] * 6024154042657686016UL) + ((uint64_t)op[5] * 6203048794235726393UL) + ((uint64_t)op[6] * 15573954401066262508UL);
	tmp_q[6] = ((uint64_t)op[0] * 10493488509697073148UL) + ((uint64_t)op[1] * 10458323379115823329UL) + ((uint64_t)op[2] * 17924663376254293634UL) + ((uint64_t)op[3] * 14401118761468768990UL) + ((uint64_t)op[4] * 7865419579248818977UL) + ((uint64_t)op[5] * 6024154042657686016UL) + ((uint64_t)op[6] * 6203048794235726393UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 23813745937L) + ((((int128)tmp_q[1] * 42533645429L) + ((int128)tmp_q[2] * 17500256115L) + ((int128)tmp_q[3] * 14954395452L) + ((int128)tmp_q[4] * 51702210316L) + ((int128)tmp_q[5] * 49708660479L) + ((int128)tmp_q[6] * 69303167285L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 69303167285L) + ((int128)tmp_q[1] * 23813745937L) + ((((int128)tmp_q[2] * 42533645429L) + ((int128)tmp_q[3] * 17500256115L) + ((int128)tmp_q[4] * 14954395452L) + ((int128)tmp_q[5] * 51702210316L) + ((int128)tmp_q[6] * 49708660479L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 49708660479L) + ((int128)tmp_q[1] * 69303167285L) + ((int128)tmp_q[2] * 23813745937L) + ((((int128)tmp_q[3] * 42533645429L) + ((int128)tmp_q[4] * 17500256115L) + ((int128)tmp_q[5] * 14954395452L) + ((int128)tmp_q[6] * 51702210316L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 51702210316L) + ((int128)tmp_q[1] * 49708660479L) + ((int128)tmp_q[2] * 69303167285L) + ((int128)tmp_q[3] * 23813745937L) + ((((int128)tmp_q[4] * 42533645429L) + ((int128)tmp_q[5] * 17500256115L) + ((int128)tmp_q[6] * 14954395452L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 14954395452L) + ((int128)tmp_q[1] * 51702210316L) + ((int128)tmp_q[2] * 49708660479L) + ((int128)tmp_q[3] * 69303167285L) + ((int128)tmp_q[4] * 23813745937L) + ((((int128)tmp_q[5] * 42533645429L) + ((int128)tmp_q[6] * 17500256115L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 17500256115L) + ((int128)tmp_q[1] * 14954395452L) + ((int128)tmp_q[2] * 51702210316L) + ((int128)tmp_q[3] * 49708660479L) + ((int128)tmp_q[4] * 69303167285L) + ((int128)tmp_q[5] * 23813745937L) + ((int128)tmp_q[6] * 212668227145L);
	tmp_zero[6] = ((int128)tmp_q[0] * 42533645429L) + ((int128)tmp_q[1] * 17500256115L) + ((int128)tmp_q[2] * 14954395452L) + ((int128)tmp_q[3] * 51702210316L) + ((int128)tmp_q[4] * 49708660479L) + ((int128)tmp_q[5] * 69303167285L) + ((int128)tmp_q[6] * 23813745937L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

