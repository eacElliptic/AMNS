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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1055718445659436695UL) + ((((uint64_t)op[1] * 3531352961957945671UL) + ((uint64_t)op[2] * 1486631529826560506UL) + ((uint64_t)op[3] * 17979628685967020513UL) + ((uint64_t)op[4] * 4084963138601363072UL) + ((uint64_t)op[5] * 922431341792761191UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 922431341792761191UL) + ((uint64_t)op[1] * 1055718445659436695UL) + ((((uint64_t)op[2] * 3531352961957945671UL) + ((uint64_t)op[3] * 1486631529826560506UL) + ((uint64_t)op[4] * 17979628685967020513UL) + ((uint64_t)op[5] * 4084963138601363072UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 4084963138601363072UL) + ((uint64_t)op[1] * 922431341792761191UL) + ((uint64_t)op[2] * 1055718445659436695UL) + ((((uint64_t)op[3] * 3531352961957945671UL) + ((uint64_t)op[4] * 1486631529826560506UL) + ((uint64_t)op[5] * 17979628685967020513UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 17979628685967020513UL) + ((uint64_t)op[1] * 4084963138601363072UL) + ((uint64_t)op[2] * 922431341792761191UL) + ((uint64_t)op[3] * 1055718445659436695UL) + ((((uint64_t)op[4] * 3531352961957945671UL) + ((uint64_t)op[5] * 1486631529826560506UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 1486631529826560506UL) + ((uint64_t)op[1] * 17979628685967020513UL) + ((uint64_t)op[2] * 4084963138601363072UL) + ((uint64_t)op[3] * 922431341792761191UL) + ((uint64_t)op[4] * 1055718445659436695UL) + ((uint64_t)op[5] * 8642664451755537864UL);
	tmp_q[5] = ((uint64_t)op[0] * 3531352961957945671UL) + ((uint64_t)op[1] * 1486631529826560506UL) + ((uint64_t)op[2] * 17979628685967020513UL) + ((uint64_t)op[3] * 4084963138601363072UL) + ((uint64_t)op[4] * 922431341792761191UL) + ((uint64_t)op[5] * 1055718445659436695UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 27707923652201L) - ((((int128)tmp_q[1] * 37146878801421L) + ((int128)tmp_q[2] * 42676648121465L) + ((int128)tmp_q[3] * 33059442598872L) + ((int128)tmp_q[4] * 76754171365137L) + ((int128)tmp_q[5] * 7243001599079L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 7243001599079L) + ((int128)tmp_q[1] * 27707923652201L) - ((((int128)tmp_q[2] * 37146878801421L) + ((int128)tmp_q[3] * 42676648121465L) + ((int128)tmp_q[4] * 33059442598872L) + ((int128)tmp_q[5] * 76754171365137L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 76754171365137L) + ((int128)tmp_q[1] * 7243001599079L) + ((int128)tmp_q[2] * 27707923652201L) - ((((int128)tmp_q[3] * 37146878801421L) + ((int128)tmp_q[4] * 42676648121465L) + ((int128)tmp_q[5] * 33059442598872L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 33059442598872L) + ((int128)tmp_q[1] * 76754171365137L) + ((int128)tmp_q[2] * 7243001599079L) + ((int128)tmp_q[3] * 27707923652201L) - ((((int128)tmp_q[4] * 37146878801421L) + ((int128)tmp_q[5] * 42676648121465L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 42676648121465L) + ((int128)tmp_q[1] * 33059442598872L) + ((int128)tmp_q[2] * 76754171365137L) + ((int128)tmp_q[3] * 7243001599079L) + ((int128)tmp_q[4] * 27707923652201L) - ((int128)tmp_q[5] * 297175030411368L);
	tmp_zero[5] = ((int128)tmp_q[0] * 37146878801421L) + ((int128)tmp_q[1] * 42676648121465L) + ((int128)tmp_q[2] * 33059442598872L) + ((int128)tmp_q[3] * 76754171365137L) + ((int128)tmp_q[4] * 7243001599079L) + ((int128)tmp_q[5] * 27707923652201L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

