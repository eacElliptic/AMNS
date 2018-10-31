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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18180271275750473243UL) + ((((uint64_t)op[1] * 6678436866350650078UL) + ((uint64_t)op[2] * 14112840067739541616UL) + ((uint64_t)op[3] * 970040583392188313UL) + ((uint64_t)op[4] * 7951402223663005635UL) + ((uint64_t)op[5] * 17369630280631895697UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 17369630280631895697UL) + ((uint64_t)op[1] * 18180271275750473243UL) + ((((uint64_t)op[2] * 6678436866350650078UL) + ((uint64_t)op[3] * 14112840067739541616UL) + ((uint64_t)op[4] * 970040583392188313UL) + ((uint64_t)op[5] * 7951402223663005635UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 7951402223663005635UL) + ((uint64_t)op[1] * 17369630280631895697UL) + ((uint64_t)op[2] * 18180271275750473243UL) + ((((uint64_t)op[3] * 6678436866350650078UL) + ((uint64_t)op[4] * 14112840067739541616UL) + ((uint64_t)op[5] * 970040583392188313UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 970040583392188313UL) + ((uint64_t)op[1] * 7951402223663005635UL) + ((uint64_t)op[2] * 17369630280631895697UL) + ((uint64_t)op[3] * 18180271275750473243UL) + ((((uint64_t)op[4] * 6678436866350650078UL) + ((uint64_t)op[5] * 14112840067739541616UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 14112840067739541616UL) + ((uint64_t)op[1] * 970040583392188313UL) + ((uint64_t)op[2] * 7951402223663005635UL) + ((uint64_t)op[3] * 17369630280631895697UL) + ((uint64_t)op[4] * 18180271275750473243UL) + ((uint64_t)op[5] * 15269611023024754380UL);
	tmp_q[5] = ((uint64_t)op[0] * 6678436866350650078UL) + ((uint64_t)op[1] * 14112840067739541616UL) + ((uint64_t)op[2] * 970040583392188313UL) + ((uint64_t)op[3] * 7951402223663005635UL) + ((uint64_t)op[4] * 17369630280631895697UL) + ((uint64_t)op[5] * 18180271275750473243UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 534164946215L) - ((((int128)tmp_q[1] * 2695695255163L) - ((int128)tmp_q[2] * 2440183712663L) + ((int128)tmp_q[3] * 1785766498354L) + ((int128)tmp_q[4] * 1635890069226L) + ((int128)tmp_q[5] * 4236864431697L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 4236864431697L) - ((int128)tmp_q[1] * 534164946215L) - ((((int128)tmp_q[2] * 2695695255163L) - ((int128)tmp_q[3] * 2440183712663L) + ((int128)tmp_q[4] * 1785766498354L) + ((int128)tmp_q[5] * 1635890069226L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 1635890069226L) + ((int128)tmp_q[1] * 4236864431697L) - ((int128)tmp_q[2] * 534164946215L) - ((((int128)tmp_q[3] * 2695695255163L) - ((int128)tmp_q[4] * 2440183712663L) + ((int128)tmp_q[5] * 1785766498354L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1785766498354L) + ((int128)tmp_q[1] * 1635890069226L) + ((int128)tmp_q[2] * 4236864431697L) - ((int128)tmp_q[3] * 534164946215L) - ((((int128)tmp_q[4] * 2695695255163L) - ((int128)tmp_q[5] * 2440183712663L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 2440183712663L) + ((int128)tmp_q[1] * 1785766498354L) + ((int128)tmp_q[2] * 1635890069226L) + ((int128)tmp_q[3] * 4236864431697L) - ((int128)tmp_q[4] * 534164946215L) - ((int128)tmp_q[5] * 16174171530978L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2695695255163L) - ((int128)tmp_q[1] * 2440183712663L) + ((int128)tmp_q[2] * 1785766498354L) + ((int128)tmp_q[3] * 1635890069226L) + ((int128)tmp_q[4] * 4236864431697L) - ((int128)tmp_q[5] * 534164946215L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

