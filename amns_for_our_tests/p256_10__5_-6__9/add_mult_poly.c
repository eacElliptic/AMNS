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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12469694602983317105UL) + ((((uint64_t)op[1] * 951607894289758993UL) + ((uint64_t)op[2] * 16723032593347104242UL) + ((uint64_t)op[3] * 15895584043057354906UL) + ((uint64_t)op[4] * 1712864738748613933UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 1712864738748613933UL) + ((uint64_t)op[1] * 12469694602983317105UL) + ((((uint64_t)op[2] * 951607894289758993UL) + ((uint64_t)op[3] * 16723032593347104242UL) + ((uint64_t)op[4] * 15895584043057354906UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 15895584043057354906UL) + ((uint64_t)op[1] * 1712864738748613933UL) + ((uint64_t)op[2] * 12469694602983317105UL) + ((((uint64_t)op[3] * 951607894289758993UL) + ((uint64_t)op[4] * 16723032593347104242UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 16723032593347104242UL) + ((uint64_t)op[1] * 15895584043057354906UL) + ((uint64_t)op[2] * 1712864738748613933UL) + ((uint64_t)op[3] * 12469694602983317105UL) + ((uint64_t)op[4] * 12737096707970997658UL);
	tmp_q[4] = ((uint64_t)op[0] * 951607894289758993UL) + ((uint64_t)op[1] * 16723032593347104242UL) + ((uint64_t)op[2] * 15895584043057354906UL) + ((uint64_t)op[3] * 1712864738748613933UL) + ((uint64_t)op[4] * 12469694602983317105UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 446582016045639L) - ((((int128)tmp_q[1] * 2364280385481318L) - ((int128)tmp_q[2] * 537965064449347L) + ((int128)tmp_q[3] * 749770862652887L) + ((int128)tmp_q[4] * 2490853758843993L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 2490853758843993L) - ((int128)tmp_q[1] * 446582016045639L) - ((((int128)tmp_q[2] * 2364280385481318L) - ((int128)tmp_q[3] * 537965064449347L) + ((int128)tmp_q[4] * 749770862652887L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 749770862652887L) + ((int128)tmp_q[1] * 2490853758843993L) - ((int128)tmp_q[2] * 446582016045639L) - ((((int128)tmp_q[3] * 2364280385481318L) - ((int128)tmp_q[4] * 537965064449347L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 537965064449347L) + ((int128)tmp_q[1] * 749770862652887L) + ((int128)tmp_q[2] * 2490853758843993L) - ((int128)tmp_q[3] * 446582016045639L) - ((int128)tmp_q[4] * 14185682312887908L);
	tmp_zero[4] = ((int128)tmp_q[0] * 2364280385481318L) - ((int128)tmp_q[1] * 537965064449347L) + ((int128)tmp_q[2] * 749770862652887L) + ((int128)tmp_q[3] * 2490853758843993L) - ((int128)tmp_q[4] * 446582016045639L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

