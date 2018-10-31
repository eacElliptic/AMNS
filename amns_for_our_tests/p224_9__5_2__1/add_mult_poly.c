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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15997305282925435235UL) + ((((uint64_t)op[1] * 8644728306076232268UL) + ((uint64_t)op[2] * 10748588077477212346UL) + ((uint64_t)op[3] * 9217331080336101187UL) + ((uint64_t)op[4] * 3107332011317222984UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 3107332011317222984UL) + ((uint64_t)op[1] * 15997305282925435235UL) + ((((uint64_t)op[2] * 8644728306076232268UL) + ((uint64_t)op[3] * 10748588077477212346UL) + ((uint64_t)op[4] * 9217331080336101187UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 9217331080336101187UL) + ((uint64_t)op[1] * 3107332011317222984UL) + ((uint64_t)op[2] * 15997305282925435235UL) + ((((uint64_t)op[3] * 8644728306076232268UL) + ((uint64_t)op[4] * 10748588077477212346UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 10748588077477212346UL) + ((uint64_t)op[1] * 9217331080336101187UL) + ((uint64_t)op[2] * 3107332011317222984UL) + ((uint64_t)op[3] * 15997305282925435235UL) + ((uint64_t)op[4] * 17289456612152464536UL);
	tmp_q[4] = ((uint64_t)op[0] * 8644728306076232268UL) + ((uint64_t)op[1] * 10748588077477212346UL) + ((uint64_t)op[2] * 9217331080336101187UL) + ((uint64_t)op[3] * 3107332011317222984UL) + ((uint64_t)op[4] * 15997305282925435235UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 7690609000823L) + ((((int128)tmp_q[1] * 14502883148621L) - ((int128)tmp_q[2] * 6997426790020L) - ((int128)tmp_q[3] * 1941216538981L) + ((int128)tmp_q[4] * 19888026055910L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 19888026055910L) - ((int128)tmp_q[1] * 7690609000823L) + ((((int128)tmp_q[2] * 14502883148621L) - ((int128)tmp_q[3] * 6997426790020L) - ((int128)tmp_q[4] * 1941216538981L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 1941216538981L) + ((int128)tmp_q[1] * 19888026055910L) - ((int128)tmp_q[2] * 7690609000823L) + ((((int128)tmp_q[3] * 14502883148621L) - ((int128)tmp_q[4] * 6997426790020L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 6997426790020L) - ((int128)tmp_q[1] * 1941216538981L) + ((int128)tmp_q[2] * 19888026055910L) - ((int128)tmp_q[3] * 7690609000823L) + ((int128)tmp_q[4] * 29005766297242L);
	tmp_zero[4] = ((int128)tmp_q[0] * 14502883148621L) - ((int128)tmp_q[1] * 6997426790020L) - ((int128)tmp_q[2] * 1941216538981L) + ((int128)tmp_q[3] * 19888026055910L) - ((int128)tmp_q[4] * 7690609000823L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

