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
	tmp_q[0] = ((uint64_t)op[0] * 14111439465940304395UL) + ((((uint64_t)op[1] * 2359354419849899615UL) + ((uint64_t)op[2] * 11247921608585129018UL) + ((uint64_t)op[3] * 16185105432950119797UL) + ((uint64_t)op[4] * 1769361421402385710UL) + ((uint64_t)op[5] * 13959764442608074216UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 13959764442608074216UL) + ((uint64_t)op[1] * 14111439465940304395UL) + ((((uint64_t)op[2] * 2359354419849899615UL) + ((uint64_t)op[3] * 11247921608585129018UL) + ((uint64_t)op[4] * 16185105432950119797UL) + ((uint64_t)op[5] * 1769361421402385710UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 1769361421402385710UL) + ((uint64_t)op[1] * 13959764442608074216UL) + ((uint64_t)op[2] * 14111439465940304395UL) + ((((uint64_t)op[3] * 2359354419849899615UL) + ((uint64_t)op[4] * 11247921608585129018UL) + ((uint64_t)op[5] * 16185105432950119797UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 16185105432950119797UL) + ((uint64_t)op[1] * 1769361421402385710UL) + ((uint64_t)op[2] * 13959764442608074216UL) + ((uint64_t)op[3] * 14111439465940304395UL) + ((((uint64_t)op[4] * 2359354419849899615UL) + ((uint64_t)op[5] * 11247921608585129018UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 11247921608585129018UL) + ((uint64_t)op[1] * 16185105432950119797UL) + ((uint64_t)op[2] * 1769361421402385710UL) + ((uint64_t)op[3] * 13959764442608074216UL) + ((uint64_t)op[4] * 14111439465940304395UL) + ((uint64_t)op[5] * 4290617554610153926UL);
	tmp_q[5] = ((uint64_t)op[0] * 2359354419849899615UL) + ((uint64_t)op[1] * 11247921608585129018UL) + ((uint64_t)op[2] * 16185105432950119797UL) + ((uint64_t)op[3] * 1769361421402385710UL) + ((uint64_t)op[4] * 13959764442608074216UL) + ((uint64_t)op[5] * 14111439465940304395UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2940318046925L) - ((-((int128)tmp_q[1] * 1156744925367L) + ((int128)tmp_q[2] * 4158934301500L) - ((int128)tmp_q[3] * 886458088121L) + ((int128)tmp_q[4] * 2146296991262L) + ((int128)tmp_q[5] * 3867803861116L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 3867803861116L) - ((int128)tmp_q[1] * 2940318046925L) - ((-((int128)tmp_q[2] * 1156744925367L) + ((int128)tmp_q[3] * 4158934301500L) - ((int128)tmp_q[4] * 886458088121L) + ((int128)tmp_q[5] * 2146296991262L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 2146296991262L) + ((int128)tmp_q[1] * 3867803861116L) - ((int128)tmp_q[2] * 2940318046925L) - ((-((int128)tmp_q[3] * 1156744925367L) + ((int128)tmp_q[4] * 4158934301500L) - ((int128)tmp_q[5] * 886458088121L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 886458088121L) + ((int128)tmp_q[1] * 2146296991262L) + ((int128)tmp_q[2] * 3867803861116L) - ((int128)tmp_q[3] * 2940318046925L) - ((-((int128)tmp_q[4] * 1156744925367L) + ((int128)tmp_q[5] * 4158934301500L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 4158934301500L) - ((int128)tmp_q[1] * 886458088121L) + ((int128)tmp_q[2] * 2146296991262L) + ((int128)tmp_q[3] * 3867803861116L) - ((int128)tmp_q[4] * 2940318046925L) + ((int128)tmp_q[5] * 6940469552202L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1156744925367L) + ((int128)tmp_q[1] * 4158934301500L) - ((int128)tmp_q[2] * 886458088121L) + ((int128)tmp_q[3] * 2146296991262L) + ((int128)tmp_q[4] * 3867803861116L) - ((int128)tmp_q[5] * 2940318046925L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

