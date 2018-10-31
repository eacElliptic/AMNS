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
	tmp_q[0] = ((uint64_t)op[0] * 9575239449420039377UL) + ((((uint64_t)op[1] * 18200779438629337414UL) + ((uint64_t)op[2] * 13488360666086598926UL) + ((uint64_t)op[3] * 9789571939015084809UL) + ((uint64_t)op[4] * 16732811939981982179UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16732811939981982179UL) + ((uint64_t)op[1] * 9575239449420039377UL) + ((((uint64_t)op[2] * 18200779438629337414UL) + ((uint64_t)op[3] * 13488360666086598926UL) + ((uint64_t)op[4] * 9789571939015084809UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 9789571939015084809UL) + ((uint64_t)op[1] * 16732811939981982179UL) + ((uint64_t)op[2] * 9575239449420039377UL) + ((((uint64_t)op[3] * 18200779438629337414UL) + ((uint64_t)op[4] * 13488360666086598926UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 13488360666086598926UL) + ((uint64_t)op[1] * 9789571939015084809UL) + ((uint64_t)op[2] * 16732811939981982179UL) + ((uint64_t)op[3] * 9575239449420039377UL) + ((uint64_t)op[4] * 1475787810481285212UL);
	tmp_q[4] = ((uint64_t)op[0] * 18200779438629337414UL) + ((uint64_t)op[1] * 13488360666086598926UL) + ((uint64_t)op[2] * 9789571939015084809UL) + ((uint64_t)op[3] * 16732811939981982179UL) + ((uint64_t)op[4] * 9575239449420039377UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1671048223693017L) - ((-((int128)tmp_q[1] * 1127539810489199L) + ((int128)tmp_q[2] * 1014173957885599L) + ((int128)tmp_q[3] * 295920153948514L) + ((int128)tmp_q[4] * 1363850448796389L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 1363850448796389L) - ((int128)tmp_q[1] * 1671048223693017L) - ((-((int128)tmp_q[2] * 1127539810489199L) + ((int128)tmp_q[3] * 1014173957885599L) + ((int128)tmp_q[4] * 295920153948514L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 295920153948514L) + ((int128)tmp_q[1] * 1363850448796389L) - ((int128)tmp_q[2] * 1671048223693017L) - ((-((int128)tmp_q[3] * 1127539810489199L) + ((int128)tmp_q[4] * 1014173957885599L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1014173957885599L) + ((int128)tmp_q[1] * 295920153948514L) + ((int128)tmp_q[2] * 1363850448796389L) - ((int128)tmp_q[3] * 1671048223693017L) + ((int128)tmp_q[4] * 6765238862935194L);
	tmp_zero[4] = -((int128)tmp_q[0] * 1127539810489199L) + ((int128)tmp_q[1] * 1014173957885599L) + ((int128)tmp_q[2] * 295920153948514L) + ((int128)tmp_q[3] * 1363850448796389L) - ((int128)tmp_q[4] * 1671048223693017L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

