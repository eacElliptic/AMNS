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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9895543390708696005UL) + ((((uint64_t)op[1] * 15194578426280676420UL) + ((uint64_t)op[2] * 5614683772940367772UL) + ((uint64_t)op[3] * 10540415567984782927UL) + ((uint64_t)op[4] * 17308371798148554186UL) + ((uint64_t)op[5] * 7221600412973833048UL) + ((uint64_t)op[6] * 7288126539530055614UL) + ((uint64_t)op[7] * 7971264918158028471UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 7971264918158028471UL) + ((uint64_t)op[1] * 9895543390708696005UL) + ((((uint64_t)op[2] * 15194578426280676420UL) + ((uint64_t)op[3] * 5614683772940367772UL) + ((uint64_t)op[4] * 10540415567984782927UL) + ((uint64_t)op[5] * 17308371798148554186UL) + ((uint64_t)op[6] * 7221600412973833048UL) + ((uint64_t)op[7] * 7288126539530055614UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 7288126539530055614UL) + ((uint64_t)op[1] * 7971264918158028471UL) + ((uint64_t)op[2] * 9895543390708696005UL) + ((((uint64_t)op[3] * 15194578426280676420UL) + ((uint64_t)op[4] * 5614683772940367772UL) + ((uint64_t)op[5] * 10540415567984782927UL) + ((uint64_t)op[6] * 17308371798148554186UL) + ((uint64_t)op[7] * 7221600412973833048UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 7221600412973833048UL) + ((uint64_t)op[1] * 7288126539530055614UL) + ((uint64_t)op[2] * 7971264918158028471UL) + ((uint64_t)op[3] * 9895543390708696005UL) + ((((uint64_t)op[4] * 15194578426280676420UL) + ((uint64_t)op[5] * 5614683772940367772UL) + ((uint64_t)op[6] * 10540415567984782927UL) + ((uint64_t)op[7] * 17308371798148554186UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 17308371798148554186UL) + ((uint64_t)op[1] * 7221600412973833048UL) + ((uint64_t)op[2] * 7288126539530055614UL) + ((uint64_t)op[3] * 7971264918158028471UL) + ((uint64_t)op[4] * 9895543390708696005UL) + ((((uint64_t)op[5] * 15194578426280676420UL) + ((uint64_t)op[6] * 5614683772940367772UL) + ((uint64_t)op[7] * 10540415567984782927UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 10540415567984782927UL) + ((uint64_t)op[1] * 17308371798148554186UL) + ((uint64_t)op[2] * 7221600412973833048UL) + ((uint64_t)op[3] * 7288126539530055614UL) + ((uint64_t)op[4] * 7971264918158028471UL) + ((uint64_t)op[5] * 9895543390708696005UL) + ((((uint64_t)op[6] * 15194578426280676420UL) + ((uint64_t)op[7] * 5614683772940367772UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 5614683772940367772UL) + ((uint64_t)op[1] * 10540415567984782927UL) + ((uint64_t)op[2] * 17308371798148554186UL) + ((uint64_t)op[3] * 7221600412973833048UL) + ((uint64_t)op[4] * 7288126539530055614UL) + ((uint64_t)op[5] * 7971264918158028471UL) + ((uint64_t)op[6] * 9895543390708696005UL) + ((uint64_t)op[7] * 6504331294857750392UL);
	tmp_q[7] = ((uint64_t)op[0] * 15194578426280676420UL) + ((uint64_t)op[1] * 5614683772940367772UL) + ((uint64_t)op[2] * 10540415567984782927UL) + ((uint64_t)op[3] * 17308371798148554186UL) + ((uint64_t)op[4] * 7221600412973833048UL) + ((uint64_t)op[5] * 7288126539530055614UL) + ((uint64_t)op[6] * 7971264918158028471UL) + ((uint64_t)op[7] * 9895543390708696005UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 17868999550567L) - ((-((int128)tmp_q[1] * 106641464582520L) - ((int128)tmp_q[2] * 11241434641171L) - ((int128)tmp_q[3] * 39331931954430L) + ((int128)tmp_q[4] * 77221694373685L) + ((int128)tmp_q[5] * 8853950268889L) - ((int128)tmp_q[6] * 113960445753783L) + ((int128)tmp_q[7] * 108797066620739L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 108797066620739L) - ((int128)tmp_q[1] * 17868999550567L) - ((-((int128)tmp_q[2] * 106641464582520L) - ((int128)tmp_q[3] * 11241434641171L) - ((int128)tmp_q[4] * 39331931954430L) + ((int128)tmp_q[5] * 77221694373685L) + ((int128)tmp_q[6] * 8853950268889L) - ((int128)tmp_q[7] * 113960445753783L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 113960445753783L) + ((int128)tmp_q[1] * 108797066620739L) - ((int128)tmp_q[2] * 17868999550567L) - ((-((int128)tmp_q[3] * 106641464582520L) - ((int128)tmp_q[4] * 11241434641171L) - ((int128)tmp_q[5] * 39331931954430L) + ((int128)tmp_q[6] * 77221694373685L) + ((int128)tmp_q[7] * 8853950268889L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 8853950268889L) - ((int128)tmp_q[1] * 113960445753783L) + ((int128)tmp_q[2] * 108797066620739L) - ((int128)tmp_q[3] * 17868999550567L) - ((-((int128)tmp_q[4] * 106641464582520L) - ((int128)tmp_q[5] * 11241434641171L) - ((int128)tmp_q[6] * 39331931954430L) + ((int128)tmp_q[7] * 77221694373685L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 77221694373685L) + ((int128)tmp_q[1] * 8853950268889L) - ((int128)tmp_q[2] * 113960445753783L) + ((int128)tmp_q[3] * 108797066620739L) - ((int128)tmp_q[4] * 17868999550567L) - ((-((int128)tmp_q[5] * 106641464582520L) - ((int128)tmp_q[6] * 11241434641171L) - ((int128)tmp_q[7] * 39331931954430L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 39331931954430L) + ((int128)tmp_q[1] * 77221694373685L) + ((int128)tmp_q[2] * 8853950268889L) - ((int128)tmp_q[3] * 113960445753783L) + ((int128)tmp_q[4] * 108797066620739L) - ((int128)tmp_q[5] * 17868999550567L) - ((-((int128)tmp_q[6] * 106641464582520L) - ((int128)tmp_q[7] * 11241434641171L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 11241434641171L) - ((int128)tmp_q[1] * 39331931954430L) + ((int128)tmp_q[2] * 77221694373685L) + ((int128)tmp_q[3] * 8853950268889L) - ((int128)tmp_q[4] * 113960445753783L) + ((int128)tmp_q[5] * 108797066620739L) - ((int128)tmp_q[6] * 17868999550567L) + ((int128)tmp_q[7] * 213282929165040L);
	tmp_zero[7] = -((int128)tmp_q[0] * 106641464582520L) - ((int128)tmp_q[1] * 11241434641171L) - ((int128)tmp_q[2] * 39331931954430L) + ((int128)tmp_q[3] * 77221694373685L) + ((int128)tmp_q[4] * 8853950268889L) - ((int128)tmp_q[5] * 113960445753783L) + ((int128)tmp_q[6] * 108797066620739L) - ((int128)tmp_q[7] * 17868999550567L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

