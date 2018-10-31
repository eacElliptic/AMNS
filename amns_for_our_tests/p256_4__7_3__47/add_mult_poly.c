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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2007803222218272583UL) + ((((uint64_t)op[1] * 7730191527778639246UL) + ((uint64_t)op[2] * 10435629194004976732UL) + ((uint64_t)op[3] * 3886656563343595630UL) + ((uint64_t)op[4] * 2810704570626721312UL) + ((uint64_t)op[5] * 5931462719672619731UL) + ((uint64_t)op[6] * 6519101444511327421UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 6519101444511327421UL) + ((uint64_t)op[1] * 2007803222218272583UL) + ((((uint64_t)op[2] * 7730191527778639246UL) + ((uint64_t)op[3] * 10435629194004976732UL) + ((uint64_t)op[4] * 3886656563343595630UL) + ((uint64_t)op[5] * 2810704570626721312UL) + ((uint64_t)op[6] * 5931462719672619731UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 5931462719672619731UL) + ((uint64_t)op[1] * 6519101444511327421UL) + ((uint64_t)op[2] * 2007803222218272583UL) + ((((uint64_t)op[3] * 7730191527778639246UL) + ((uint64_t)op[4] * 10435629194004976732UL) + ((uint64_t)op[5] * 3886656563343595630UL) + ((uint64_t)op[6] * 2810704570626721312UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 2810704570626721312UL) + ((uint64_t)op[1] * 5931462719672619731UL) + ((uint64_t)op[2] * 6519101444511327421UL) + ((uint64_t)op[3] * 2007803222218272583UL) + ((((uint64_t)op[4] * 7730191527778639246UL) + ((uint64_t)op[5] * 10435629194004976732UL) + ((uint64_t)op[6] * 3886656563343595630UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 3886656563343595630UL) + ((uint64_t)op[1] * 2810704570626721312UL) + ((uint64_t)op[2] * 5931462719672619731UL) + ((uint64_t)op[3] * 6519101444511327421UL) + ((uint64_t)op[4] * 2007803222218272583UL) + ((((uint64_t)op[5] * 7730191527778639246UL) + ((uint64_t)op[6] * 10435629194004976732UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 10435629194004976732UL) + ((uint64_t)op[1] * 3886656563343595630UL) + ((uint64_t)op[2] * 2810704570626721312UL) + ((uint64_t)op[3] * 5931462719672619731UL) + ((uint64_t)op[4] * 6519101444511327421UL) + ((uint64_t)op[5] * 2007803222218272583UL) + ((uint64_t)op[6] * 4743830509626366122UL);
	tmp_q[6] = ((uint64_t)op[0] * 7730191527778639246UL) + ((uint64_t)op[1] * 10435629194004976732UL) + ((uint64_t)op[2] * 3886656563343595630UL) + ((uint64_t)op[3] * 2810704570626721312UL) + ((uint64_t)op[4] * 5931462719672619731UL) + ((uint64_t)op[5] * 6519101444511327421UL) + ((uint64_t)op[6] * 2007803222218272583UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 23253353873L) + ((-((int128)tmp_q[1] * 14031500763L) + ((int128)tmp_q[2] * 53112075975L) - ((int128)tmp_q[3] * 47885737810L) - ((int128)tmp_q[4] * 29741883015L) + ((int128)tmp_q[5] * 5992630027L) - ((int128)tmp_q[6] * 1037469192L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 1037469192L) + ((int128)tmp_q[1] * 23253353873L) + ((-((int128)tmp_q[2] * 14031500763L) + ((int128)tmp_q[3] * 53112075975L) - ((int128)tmp_q[4] * 47885737810L) - ((int128)tmp_q[5] * 29741883015L) + ((int128)tmp_q[6] * 5992630027L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 5992630027L) - ((int128)tmp_q[1] * 1037469192L) + ((int128)tmp_q[2] * 23253353873L) + ((-((int128)tmp_q[3] * 14031500763L) + ((int128)tmp_q[4] * 53112075975L) - ((int128)tmp_q[5] * 47885737810L) - ((int128)tmp_q[6] * 29741883015L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 29741883015L) + ((int128)tmp_q[1] * 5992630027L) - ((int128)tmp_q[2] * 1037469192L) + ((int128)tmp_q[3] * 23253353873L) + ((-((int128)tmp_q[4] * 14031500763L) + ((int128)tmp_q[5] * 53112075975L) - ((int128)tmp_q[6] * 47885737810L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 47885737810L) - ((int128)tmp_q[1] * 29741883015L) + ((int128)tmp_q[2] * 5992630027L) - ((int128)tmp_q[3] * 1037469192L) + ((int128)tmp_q[4] * 23253353873L) + ((-((int128)tmp_q[5] * 14031500763L) + ((int128)tmp_q[6] * 53112075975L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 53112075975L) - ((int128)tmp_q[1] * 47885737810L) - ((int128)tmp_q[2] * 29741883015L) + ((int128)tmp_q[3] * 5992630027L) - ((int128)tmp_q[4] * 1037469192L) + ((int128)tmp_q[5] * 23253353873L) - ((int128)tmp_q[6] * 42094502289L);
	tmp_zero[6] = -((int128)tmp_q[0] * 14031500763L) + ((int128)tmp_q[1] * 53112075975L) - ((int128)tmp_q[2] * 47885737810L) - ((int128)tmp_q[3] * 29741883015L) + ((int128)tmp_q[4] * 5992630027L) - ((int128)tmp_q[5] * 1037469192L) + ((int128)tmp_q[6] * 23253353873L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

