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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2085257375062566816UL) + ((((uint64_t)op[1] * 11971594155220829825UL) + ((uint64_t)op[2] * 9295415239576128668UL) + ((uint64_t)op[3] * 768659598809332199UL) + ((uint64_t)op[4] * 3018753352623969973UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 3018753352623969973UL) + ((uint64_t)op[1] * 2085257375062566816UL) + ((((uint64_t)op[2] * 11971594155220829825UL) + ((uint64_t)op[3] * 9295415239576128668UL) + ((uint64_t)op[4] * 768659598809332199UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 768659598809332199UL) + ((uint64_t)op[1] * 3018753352623969973UL) + ((uint64_t)op[2] * 2085257375062566816UL) + ((((uint64_t)op[3] * 11971594155220829825UL) + ((uint64_t)op[4] * 9295415239576128668UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 9295415239576128668UL) + ((uint64_t)op[1] * 768659598809332199UL) + ((uint64_t)op[2] * 3018753352623969973UL) + ((uint64_t)op[3] * 2085257375062566816UL) + ((uint64_t)op[4] * 4517738554975494277UL);
	tmp_q[4] = ((uint64_t)op[0] * 11971594155220829825UL) + ((uint64_t)op[1] * 9295415239576128668UL) + ((uint64_t)op[2] * 768659598809332199UL) + ((uint64_t)op[3] * 3018753352623969973UL) + ((uint64_t)op[4] * 2085257375062566816UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 358671532963L) + ((((int128)tmp_q[1] * 220136546860L) - ((int128)tmp_q[2] * 146148960478L) - ((int128)tmp_q[3] * 186001793023L) - ((int128)tmp_q[4] * 189798319635L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 189798319635L) + ((int128)tmp_q[1] * 358671532963L) + ((((int128)tmp_q[2] * 220136546860L) - ((int128)tmp_q[3] * 146148960478L) - ((int128)tmp_q[4] * 186001793023L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 186001793023L) - ((int128)tmp_q[1] * 189798319635L) + ((int128)tmp_q[2] * 358671532963L) + ((((int128)tmp_q[3] * 220136546860L) - ((int128)tmp_q[4] * 146148960478L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 146148960478L) - ((int128)tmp_q[1] * 186001793023L) - ((int128)tmp_q[2] * 189798319635L) + ((int128)tmp_q[3] * 358671532963L) + ((int128)tmp_q[4] * 1100682734300L);
	tmp_zero[4] = ((int128)tmp_q[0] * 220136546860L) - ((int128)tmp_q[1] * 146148960478L) - ((int128)tmp_q[2] * 186001793023L) - ((int128)tmp_q[3] * 189798319635L) + ((int128)tmp_q[4] * 358671532963L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

