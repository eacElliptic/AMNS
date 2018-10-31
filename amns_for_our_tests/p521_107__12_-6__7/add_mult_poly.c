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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 6);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 12);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 12);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 6);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9600161852058675165UL) + ((((uint64_t)op[1] * 13634178712876501333UL) + ((uint64_t)op[2] * 4665101481282005532UL) + ((uint64_t)op[3] * 8309209697472775201UL) + ((uint64_t)op[4] * 1569882885229354159UL) + ((uint64_t)op[5] * 17521546868388607043UL) + ((uint64_t)op[6] * 14964450715530679832UL) + ((uint64_t)op[7] * 8214389023539642784UL) + ((uint64_t)op[8] * 12202411481317448553UL) + ((uint64_t)op[9] * 11571515075593750888UL) + ((uint64_t)op[10] * 16867362105146066712UL) + ((uint64_t)op[11] * 2685244346811222983UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 2685244346811222983UL) + ((uint64_t)op[1] * 9600161852058675165UL) + ((((uint64_t)op[2] * 13634178712876501333UL) + ((uint64_t)op[3] * 4665101481282005532UL) + ((uint64_t)op[4] * 8309209697472775201UL) + ((uint64_t)op[5] * 1569882885229354159UL) + ((uint64_t)op[6] * 17521546868388607043UL) + ((uint64_t)op[7] * 14964450715530679832UL) + ((uint64_t)op[8] * 8214389023539642784UL) + ((uint64_t)op[9] * 12202411481317448553UL) + ((uint64_t)op[10] * 11571515075593750888UL) + ((uint64_t)op[11] * 16867362105146066712UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 16867362105146066712UL) + ((uint64_t)op[1] * 2685244346811222983UL) + ((uint64_t)op[2] * 9600161852058675165UL) + ((((uint64_t)op[3] * 13634178712876501333UL) + ((uint64_t)op[4] * 4665101481282005532UL) + ((uint64_t)op[5] * 8309209697472775201UL) + ((uint64_t)op[6] * 1569882885229354159UL) + ((uint64_t)op[7] * 17521546868388607043UL) + ((uint64_t)op[8] * 14964450715530679832UL) + ((uint64_t)op[9] * 8214389023539642784UL) + ((uint64_t)op[10] * 12202411481317448553UL) + ((uint64_t)op[11] * 11571515075593750888UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 11571515075593750888UL) + ((uint64_t)op[1] * 16867362105146066712UL) + ((uint64_t)op[2] * 2685244346811222983UL) + ((uint64_t)op[3] * 9600161852058675165UL) + ((((uint64_t)op[4] * 13634178712876501333UL) + ((uint64_t)op[5] * 4665101481282005532UL) + ((uint64_t)op[6] * 8309209697472775201UL) + ((uint64_t)op[7] * 1569882885229354159UL) + ((uint64_t)op[8] * 17521546868388607043UL) + ((uint64_t)op[9] * 14964450715530679832UL) + ((uint64_t)op[10] * 8214389023539642784UL) + ((uint64_t)op[11] * 12202411481317448553UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 12202411481317448553UL) + ((uint64_t)op[1] * 11571515075593750888UL) + ((uint64_t)op[2] * 16867362105146066712UL) + ((uint64_t)op[3] * 2685244346811222983UL) + ((uint64_t)op[4] * 9600161852058675165UL) + ((((uint64_t)op[5] * 13634178712876501333UL) + ((uint64_t)op[6] * 4665101481282005532UL) + ((uint64_t)op[7] * 8309209697472775201UL) + ((uint64_t)op[8] * 1569882885229354159UL) + ((uint64_t)op[9] * 17521546868388607043UL) + ((uint64_t)op[10] * 14964450715530679832UL) + ((uint64_t)op[11] * 8214389023539642784UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 8214389023539642784UL) + ((uint64_t)op[1] * 12202411481317448553UL) + ((uint64_t)op[2] * 11571515075593750888UL) + ((uint64_t)op[3] * 16867362105146066712UL) + ((uint64_t)op[4] * 2685244346811222983UL) + ((uint64_t)op[5] * 9600161852058675165UL) + ((((uint64_t)op[6] * 13634178712876501333UL) + ((uint64_t)op[7] * 4665101481282005532UL) + ((uint64_t)op[8] * 8309209697472775201UL) + ((uint64_t)op[9] * 1569882885229354159UL) + ((uint64_t)op[10] * 17521546868388607043UL) + ((uint64_t)op[11] * 14964450715530679832UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 14964450715530679832UL) + ((uint64_t)op[1] * 8214389023539642784UL) + ((uint64_t)op[2] * 12202411481317448553UL) + ((uint64_t)op[3] * 11571515075593750888UL) + ((uint64_t)op[4] * 16867362105146066712UL) + ((uint64_t)op[5] * 2685244346811222983UL) + ((uint64_t)op[6] * 9600161852058675165UL) + ((((uint64_t)op[7] * 13634178712876501333UL) + ((uint64_t)op[8] * 4665101481282005532UL) + ((uint64_t)op[9] * 8309209697472775201UL) + ((uint64_t)op[10] * 1569882885229354159UL) + ((uint64_t)op[11] * 17521546868388607043UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 17521546868388607043UL) + ((uint64_t)op[1] * 14964450715530679832UL) + ((uint64_t)op[2] * 8214389023539642784UL) + ((uint64_t)op[3] * 12202411481317448553UL) + ((uint64_t)op[4] * 11571515075593750888UL) + ((uint64_t)op[5] * 16867362105146066712UL) + ((uint64_t)op[6] * 2685244346811222983UL) + ((uint64_t)op[7] * 9600161852058675165UL) + ((((uint64_t)op[8] * 13634178712876501333UL) + ((uint64_t)op[9] * 4665101481282005532UL) + ((uint64_t)op[10] * 8309209697472775201UL) + ((uint64_t)op[11] * 1569882885229354159UL)) * 18446744073709551610);
	tmp_q[8] = ((uint64_t)op[0] * 1569882885229354159UL) + ((uint64_t)op[1] * 17521546868388607043UL) + ((uint64_t)op[2] * 14964450715530679832UL) + ((uint64_t)op[3] * 8214389023539642784UL) + ((uint64_t)op[4] * 12202411481317448553UL) + ((uint64_t)op[5] * 11571515075593750888UL) + ((uint64_t)op[6] * 16867362105146066712UL) + ((uint64_t)op[7] * 2685244346811222983UL) + ((uint64_t)op[8] * 9600161852058675165UL) + ((((uint64_t)op[9] * 13634178712876501333UL) + ((uint64_t)op[10] * 4665101481282005532UL) + ((uint64_t)op[11] * 8309209697472775201UL)) * 18446744073709551610);
	tmp_q[9] = ((uint64_t)op[0] * 8309209697472775201UL) + ((uint64_t)op[1] * 1569882885229354159UL) + ((uint64_t)op[2] * 17521546868388607043UL) + ((uint64_t)op[3] * 14964450715530679832UL) + ((uint64_t)op[4] * 8214389023539642784UL) + ((uint64_t)op[5] * 12202411481317448553UL) + ((uint64_t)op[6] * 11571515075593750888UL) + ((uint64_t)op[7] * 16867362105146066712UL) + ((uint64_t)op[8] * 2685244346811222983UL) + ((uint64_t)op[9] * 9600161852058675165UL) + ((((uint64_t)op[10] * 13634178712876501333UL) + ((uint64_t)op[11] * 4665101481282005532UL)) * 18446744073709551610);
	tmp_q[10] = ((uint64_t)op[0] * 4665101481282005532UL) + ((uint64_t)op[1] * 8309209697472775201UL) + ((uint64_t)op[2] * 1569882885229354159UL) + ((uint64_t)op[3] * 17521546868388607043UL) + ((uint64_t)op[4] * 14964450715530679832UL) + ((uint64_t)op[5] * 8214389023539642784UL) + ((uint64_t)op[6] * 12202411481317448553UL) + ((uint64_t)op[7] * 11571515075593750888UL) + ((uint64_t)op[8] * 16867362105146066712UL) + ((uint64_t)op[9] * 2685244346811222983UL) + ((uint64_t)op[10] * 9600161852058675165UL) + ((uint64_t)op[11] * 10428648091288750082UL);
	tmp_q[11] = ((uint64_t)op[0] * 13634178712876501333UL) + ((uint64_t)op[1] * 4665101481282005532UL) + ((uint64_t)op[2] * 8309209697472775201UL) + ((uint64_t)op[3] * 1569882885229354159UL) + ((uint64_t)op[4] * 17521546868388607043UL) + ((uint64_t)op[5] * 14964450715530679832UL) + ((uint64_t)op[6] * 8214389023539642784UL) + ((uint64_t)op[7] * 12202411481317448553UL) + ((uint64_t)op[8] * 11571515075593750888UL) + ((uint64_t)op[9] * 16867362105146066712UL) + ((uint64_t)op[10] * 2685244346811222983UL) + ((uint64_t)op[11] * 9600161852058675165UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2019131915395L) - ((((int128)tmp_q[1] * 5246195730242L) + ((int128)tmp_q[2] * 4540522082875L) - ((int128)tmp_q[3] * 252367658450L) - ((int128)tmp_q[4] * 4634026393462L) - ((int128)tmp_q[5] * 6913985487204L) - ((int128)tmp_q[6] * 2891002178804L) - ((int128)tmp_q[7] * 1077240298019L) - ((int128)tmp_q[8] * 4753277710636L) + ((int128)tmp_q[9] * 3035936621427L) - ((int128)tmp_q[10] * 6515335234081L) - ((int128)tmp_q[11] * 627495443143L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 627495443143L) - ((int128)tmp_q[1] * 2019131915395L) - ((((int128)tmp_q[2] * 5246195730242L) + ((int128)tmp_q[3] * 4540522082875L) - ((int128)tmp_q[4] * 252367658450L) - ((int128)tmp_q[5] * 4634026393462L) - ((int128)tmp_q[6] * 6913985487204L) - ((int128)tmp_q[7] * 2891002178804L) - ((int128)tmp_q[8] * 1077240298019L) - ((int128)tmp_q[9] * 4753277710636L) + ((int128)tmp_q[10] * 3035936621427L) - ((int128)tmp_q[11] * 6515335234081L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 6515335234081L) - ((int128)tmp_q[1] * 627495443143L) - ((int128)tmp_q[2] * 2019131915395L) - ((((int128)tmp_q[3] * 5246195730242L) + ((int128)tmp_q[4] * 4540522082875L) - ((int128)tmp_q[5] * 252367658450L) - ((int128)tmp_q[6] * 4634026393462L) - ((int128)tmp_q[7] * 6913985487204L) - ((int128)tmp_q[8] * 2891002178804L) - ((int128)tmp_q[9] * 1077240298019L) - ((int128)tmp_q[10] * 4753277710636L) + ((int128)tmp_q[11] * 3035936621427L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 3035936621427L) - ((int128)tmp_q[1] * 6515335234081L) - ((int128)tmp_q[2] * 627495443143L) - ((int128)tmp_q[3] * 2019131915395L) - ((((int128)tmp_q[4] * 5246195730242L) + ((int128)tmp_q[5] * 4540522082875L) - ((int128)tmp_q[6] * 252367658450L) - ((int128)tmp_q[7] * 4634026393462L) - ((int128)tmp_q[8] * 6913985487204L) - ((int128)tmp_q[9] * 2891002178804L) - ((int128)tmp_q[10] * 1077240298019L) - ((int128)tmp_q[11] * 4753277710636L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 4753277710636L) + ((int128)tmp_q[1] * 3035936621427L) - ((int128)tmp_q[2] * 6515335234081L) - ((int128)tmp_q[3] * 627495443143L) - ((int128)tmp_q[4] * 2019131915395L) - ((((int128)tmp_q[5] * 5246195730242L) + ((int128)tmp_q[6] * 4540522082875L) - ((int128)tmp_q[7] * 252367658450L) - ((int128)tmp_q[8] * 4634026393462L) - ((int128)tmp_q[9] * 6913985487204L) - ((int128)tmp_q[10] * 2891002178804L) - ((int128)tmp_q[11] * 1077240298019L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 1077240298019L) - ((int128)tmp_q[1] * 4753277710636L) + ((int128)tmp_q[2] * 3035936621427L) - ((int128)tmp_q[3] * 6515335234081L) - ((int128)tmp_q[4] * 627495443143L) - ((int128)tmp_q[5] * 2019131915395L) - ((((int128)tmp_q[6] * 5246195730242L) + ((int128)tmp_q[7] * 4540522082875L) - ((int128)tmp_q[8] * 252367658450L) - ((int128)tmp_q[9] * 4634026393462L) - ((int128)tmp_q[10] * 6913985487204L) - ((int128)tmp_q[11] * 2891002178804L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 2891002178804L) - ((int128)tmp_q[1] * 1077240298019L) - ((int128)tmp_q[2] * 4753277710636L) + ((int128)tmp_q[3] * 3035936621427L) - ((int128)tmp_q[4] * 6515335234081L) - ((int128)tmp_q[5] * 627495443143L) - ((int128)tmp_q[6] * 2019131915395L) - ((((int128)tmp_q[7] * 5246195730242L) + ((int128)tmp_q[8] * 4540522082875L) - ((int128)tmp_q[9] * 252367658450L) - ((int128)tmp_q[10] * 4634026393462L) - ((int128)tmp_q[11] * 6913985487204L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 6913985487204L) - ((int128)tmp_q[1] * 2891002178804L) - ((int128)tmp_q[2] * 1077240298019L) - ((int128)tmp_q[3] * 4753277710636L) + ((int128)tmp_q[4] * 3035936621427L) - ((int128)tmp_q[5] * 6515335234081L) - ((int128)tmp_q[6] * 627495443143L) - ((int128)tmp_q[7] * 2019131915395L) - ((((int128)tmp_q[8] * 5246195730242L) + ((int128)tmp_q[9] * 4540522082875L) - ((int128)tmp_q[10] * 252367658450L) - ((int128)tmp_q[11] * 4634026393462L)) * 6);
	tmp_zero[8] = -((int128)tmp_q[0] * 4634026393462L) - ((int128)tmp_q[1] * 6913985487204L) - ((int128)tmp_q[2] * 2891002178804L) - ((int128)tmp_q[3] * 1077240298019L) - ((int128)tmp_q[4] * 4753277710636L) + ((int128)tmp_q[5] * 3035936621427L) - ((int128)tmp_q[6] * 6515335234081L) - ((int128)tmp_q[7] * 627495443143L) - ((int128)tmp_q[8] * 2019131915395L) - ((((int128)tmp_q[9] * 5246195730242L) + ((int128)tmp_q[10] * 4540522082875L) - ((int128)tmp_q[11] * 252367658450L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 252367658450L) - ((int128)tmp_q[1] * 4634026393462L) - ((int128)tmp_q[2] * 6913985487204L) - ((int128)tmp_q[3] * 2891002178804L) - ((int128)tmp_q[4] * 1077240298019L) - ((int128)tmp_q[5] * 4753277710636L) + ((int128)tmp_q[6] * 3035936621427L) - ((int128)tmp_q[7] * 6515335234081L) - ((int128)tmp_q[8] * 627495443143L) - ((int128)tmp_q[9] * 2019131915395L) - ((((int128)tmp_q[10] * 5246195730242L) + ((int128)tmp_q[11] * 4540522082875L)) * 6);
	tmp_zero[10] = ((int128)tmp_q[0] * 4540522082875L) - ((int128)tmp_q[1] * 252367658450L) - ((int128)tmp_q[2] * 4634026393462L) - ((int128)tmp_q[3] * 6913985487204L) - ((int128)tmp_q[4] * 2891002178804L) - ((int128)tmp_q[5] * 1077240298019L) - ((int128)tmp_q[6] * 4753277710636L) + ((int128)tmp_q[7] * 3035936621427L) - ((int128)tmp_q[8] * 6515335234081L) - ((int128)tmp_q[9] * 627495443143L) - ((int128)tmp_q[10] * 2019131915395L) - ((int128)tmp_q[11] * 31477174381452L);
	tmp_zero[11] = ((int128)tmp_q[0] * 5246195730242L) + ((int128)tmp_q[1] * 4540522082875L) - ((int128)tmp_q[2] * 252367658450L) - ((int128)tmp_q[3] * 4634026393462L) - ((int128)tmp_q[4] * 6913985487204L) - ((int128)tmp_q[5] * 2891002178804L) - ((int128)tmp_q[6] * 1077240298019L) - ((int128)tmp_q[7] * 4753277710636L) + ((int128)tmp_q[8] * 3035936621427L) - ((int128)tmp_q[9] * 6515335234081L) - ((int128)tmp_q[10] * 627495443143L) - ((int128)tmp_q[11] * 2019131915395L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

