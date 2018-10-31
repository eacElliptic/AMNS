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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4920811176427688734UL) + ((((uint64_t)op[1] * 9741493822036522217UL) + ((uint64_t)op[2] * 8411878198808107733UL) + ((uint64_t)op[3] * 5377654014161335158UL) + ((uint64_t)op[4] * 681794268024127868UL) + ((uint64_t)op[5] * 6060957619124398383UL) + ((uint64_t)op[6] * 10538162060124678713UL) + ((uint64_t)op[7] * 9064121168597557845UL) + ((uint64_t)op[8] * 10072093646113592444UL) + ((uint64_t)op[9] * 16374808041100197840UL) + ((uint64_t)op[10] * 8300450325308716808UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 8300450325308716808UL) + ((uint64_t)op[1] * 4920811176427688734UL) + ((((uint64_t)op[2] * 9741493822036522217UL) + ((uint64_t)op[3] * 8411878198808107733UL) + ((uint64_t)op[4] * 5377654014161335158UL) + ((uint64_t)op[5] * 681794268024127868UL) + ((uint64_t)op[6] * 6060957619124398383UL) + ((uint64_t)op[7] * 10538162060124678713UL) + ((uint64_t)op[8] * 9064121168597557845UL) + ((uint64_t)op[9] * 10072093646113592444UL) + ((uint64_t)op[10] * 16374808041100197840UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 16374808041100197840UL) + ((uint64_t)op[1] * 8300450325308716808UL) + ((uint64_t)op[2] * 4920811176427688734UL) + ((((uint64_t)op[3] * 9741493822036522217UL) + ((uint64_t)op[4] * 8411878198808107733UL) + ((uint64_t)op[5] * 5377654014161335158UL) + ((uint64_t)op[6] * 681794268024127868UL) + ((uint64_t)op[7] * 6060957619124398383UL) + ((uint64_t)op[8] * 10538162060124678713UL) + ((uint64_t)op[9] * 9064121168597557845UL) + ((uint64_t)op[10] * 10072093646113592444UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 10072093646113592444UL) + ((uint64_t)op[1] * 16374808041100197840UL) + ((uint64_t)op[2] * 8300450325308716808UL) + ((uint64_t)op[3] * 4920811176427688734UL) + ((((uint64_t)op[4] * 9741493822036522217UL) + ((uint64_t)op[5] * 8411878198808107733UL) + ((uint64_t)op[6] * 5377654014161335158UL) + ((uint64_t)op[7] * 681794268024127868UL) + ((uint64_t)op[8] * 6060957619124398383UL) + ((uint64_t)op[9] * 10538162060124678713UL) + ((uint64_t)op[10] * 9064121168597557845UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 9064121168597557845UL) + ((uint64_t)op[1] * 10072093646113592444UL) + ((uint64_t)op[2] * 16374808041100197840UL) + ((uint64_t)op[3] * 8300450325308716808UL) + ((uint64_t)op[4] * 4920811176427688734UL) + ((((uint64_t)op[5] * 9741493822036522217UL) + ((uint64_t)op[6] * 8411878198808107733UL) + ((uint64_t)op[7] * 5377654014161335158UL) + ((uint64_t)op[8] * 681794268024127868UL) + ((uint64_t)op[9] * 6060957619124398383UL) + ((uint64_t)op[10] * 10538162060124678713UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 10538162060124678713UL) + ((uint64_t)op[1] * 9064121168597557845UL) + ((uint64_t)op[2] * 10072093646113592444UL) + ((uint64_t)op[3] * 16374808041100197840UL) + ((uint64_t)op[4] * 8300450325308716808UL) + ((uint64_t)op[5] * 4920811176427688734UL) + ((((uint64_t)op[6] * 9741493822036522217UL) + ((uint64_t)op[7] * 8411878198808107733UL) + ((uint64_t)op[8] * 5377654014161335158UL) + ((uint64_t)op[9] * 681794268024127868UL) + ((uint64_t)op[10] * 6060957619124398383UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 6060957619124398383UL) + ((uint64_t)op[1] * 10538162060124678713UL) + ((uint64_t)op[2] * 9064121168597557845UL) + ((uint64_t)op[3] * 10072093646113592444UL) + ((uint64_t)op[4] * 16374808041100197840UL) + ((uint64_t)op[5] * 8300450325308716808UL) + ((uint64_t)op[6] * 4920811176427688734UL) + ((((uint64_t)op[7] * 9741493822036522217UL) + ((uint64_t)op[8] * 8411878198808107733UL) + ((uint64_t)op[9] * 5377654014161335158UL) + ((uint64_t)op[10] * 681794268024127868UL)) * 5);
	tmp_q[7] = ((uint64_t)op[0] * 681794268024127868UL) + ((uint64_t)op[1] * 6060957619124398383UL) + ((uint64_t)op[2] * 10538162060124678713UL) + ((uint64_t)op[3] * 9064121168597557845UL) + ((uint64_t)op[4] * 10072093646113592444UL) + ((uint64_t)op[5] * 16374808041100197840UL) + ((uint64_t)op[6] * 8300450325308716808UL) + ((uint64_t)op[7] * 4920811176427688734UL) + ((((uint64_t)op[8] * 9741493822036522217UL) + ((uint64_t)op[9] * 8411878198808107733UL) + ((uint64_t)op[10] * 5377654014161335158UL)) * 5);
	tmp_q[8] = ((uint64_t)op[0] * 5377654014161335158UL) + ((uint64_t)op[1] * 681794268024127868UL) + ((uint64_t)op[2] * 6060957619124398383UL) + ((uint64_t)op[3] * 10538162060124678713UL) + ((uint64_t)op[4] * 9064121168597557845UL) + ((uint64_t)op[5] * 10072093646113592444UL) + ((uint64_t)op[6] * 16374808041100197840UL) + ((uint64_t)op[7] * 8300450325308716808UL) + ((uint64_t)op[8] * 4920811176427688734UL) + ((((uint64_t)op[9] * 9741493822036522217UL) + ((uint64_t)op[10] * 8411878198808107733UL)) * 5);
	tmp_q[9] = ((uint64_t)op[0] * 8411878198808107733UL) + ((uint64_t)op[1] * 5377654014161335158UL) + ((uint64_t)op[2] * 681794268024127868UL) + ((uint64_t)op[3] * 6060957619124398383UL) + ((uint64_t)op[4] * 10538162060124678713UL) + ((uint64_t)op[5] * 9064121168597557845UL) + ((uint64_t)op[6] * 10072093646113592444UL) + ((uint64_t)op[7] * 16374808041100197840UL) + ((uint64_t)op[8] * 8300450325308716808UL) + ((uint64_t)op[9] * 4920811176427688734UL) + ((uint64_t)op[10] * 11813980962763507853UL);
	tmp_q[10] = ((uint64_t)op[0] * 9741493822036522217UL) + ((uint64_t)op[1] * 8411878198808107733UL) + ((uint64_t)op[2] * 5377654014161335158UL) + ((uint64_t)op[3] * 681794268024127868UL) + ((uint64_t)op[4] * 6060957619124398383UL) + ((uint64_t)op[5] * 10538162060124678713UL) + ((uint64_t)op[6] * 9064121168597557845UL) + ((uint64_t)op[7] * 10072093646113592444UL) + ((uint64_t)op[8] * 16374808041100197840UL) + ((uint64_t)op[9] * 8300450325308716808UL) + ((uint64_t)op[10] * 4920811176427688734UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 83820676633283L) + ((((int128)tmp_q[1] * 101107308146455L) - ((int128)tmp_q[2] * 90037082646790L) + ((int128)tmp_q[3] * 10158568668206L) - ((int128)tmp_q[4] * 47147755463972L) - ((int128)tmp_q[5] * 8942258500916L) - ((int128)tmp_q[6] * 10357355016976L) + ((int128)tmp_q[7] * 58927866889179L) - ((int128)tmp_q[8] * 2949931717797L) + ((int128)tmp_q[9] * 81024989924469L) + ((int128)tmp_q[10] * 22439204334350L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 22439204334350L) - ((int128)tmp_q[1] * 83820676633283L) + ((((int128)tmp_q[2] * 101107308146455L) - ((int128)tmp_q[3] * 90037082646790L) + ((int128)tmp_q[4] * 10158568668206L) - ((int128)tmp_q[5] * 47147755463972L) - ((int128)tmp_q[6] * 8942258500916L) - ((int128)tmp_q[7] * 10357355016976L) + ((int128)tmp_q[8] * 58927866889179L) - ((int128)tmp_q[9] * 2949931717797L) + ((int128)tmp_q[10] * 81024989924469L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 81024989924469L) + ((int128)tmp_q[1] * 22439204334350L) - ((int128)tmp_q[2] * 83820676633283L) + ((((int128)tmp_q[3] * 101107308146455L) - ((int128)tmp_q[4] * 90037082646790L) + ((int128)tmp_q[5] * 10158568668206L) - ((int128)tmp_q[6] * 47147755463972L) - ((int128)tmp_q[7] * 8942258500916L) - ((int128)tmp_q[8] * 10357355016976L) + ((int128)tmp_q[9] * 58927866889179L) - ((int128)tmp_q[10] * 2949931717797L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2949931717797L) + ((int128)tmp_q[1] * 81024989924469L) + ((int128)tmp_q[2] * 22439204334350L) - ((int128)tmp_q[3] * 83820676633283L) + ((((int128)tmp_q[4] * 101107308146455L) - ((int128)tmp_q[5] * 90037082646790L) + ((int128)tmp_q[6] * 10158568668206L) - ((int128)tmp_q[7] * 47147755463972L) - ((int128)tmp_q[8] * 8942258500916L) - ((int128)tmp_q[9] * 10357355016976L) + ((int128)tmp_q[10] * 58927866889179L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 58927866889179L) - ((int128)tmp_q[1] * 2949931717797L) + ((int128)tmp_q[2] * 81024989924469L) + ((int128)tmp_q[3] * 22439204334350L) - ((int128)tmp_q[4] * 83820676633283L) + ((((int128)tmp_q[5] * 101107308146455L) - ((int128)tmp_q[6] * 90037082646790L) + ((int128)tmp_q[7] * 10158568668206L) - ((int128)tmp_q[8] * 47147755463972L) - ((int128)tmp_q[9] * 8942258500916L) - ((int128)tmp_q[10] * 10357355016976L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 10357355016976L) + ((int128)tmp_q[1] * 58927866889179L) - ((int128)tmp_q[2] * 2949931717797L) + ((int128)tmp_q[3] * 81024989924469L) + ((int128)tmp_q[4] * 22439204334350L) - ((int128)tmp_q[5] * 83820676633283L) + ((((int128)tmp_q[6] * 101107308146455L) - ((int128)tmp_q[7] * 90037082646790L) + ((int128)tmp_q[8] * 10158568668206L) - ((int128)tmp_q[9] * 47147755463972L) - ((int128)tmp_q[10] * 8942258500916L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 8942258500916L) - ((int128)tmp_q[1] * 10357355016976L) + ((int128)tmp_q[2] * 58927866889179L) - ((int128)tmp_q[3] * 2949931717797L) + ((int128)tmp_q[4] * 81024989924469L) + ((int128)tmp_q[5] * 22439204334350L) - ((int128)tmp_q[6] * 83820676633283L) + ((((int128)tmp_q[7] * 101107308146455L) - ((int128)tmp_q[8] * 90037082646790L) + ((int128)tmp_q[9] * 10158568668206L) - ((int128)tmp_q[10] * 47147755463972L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 47147755463972L) - ((int128)tmp_q[1] * 8942258500916L) - ((int128)tmp_q[2] * 10357355016976L) + ((int128)tmp_q[3] * 58927866889179L) - ((int128)tmp_q[4] * 2949931717797L) + ((int128)tmp_q[5] * 81024989924469L) + ((int128)tmp_q[6] * 22439204334350L) - ((int128)tmp_q[7] * 83820676633283L) + ((((int128)tmp_q[8] * 101107308146455L) - ((int128)tmp_q[9] * 90037082646790L) + ((int128)tmp_q[10] * 10158568668206L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 10158568668206L) - ((int128)tmp_q[1] * 47147755463972L) - ((int128)tmp_q[2] * 8942258500916L) - ((int128)tmp_q[3] * 10357355016976L) + ((int128)tmp_q[4] * 58927866889179L) - ((int128)tmp_q[5] * 2949931717797L) + ((int128)tmp_q[6] * 81024989924469L) + ((int128)tmp_q[7] * 22439204334350L) - ((int128)tmp_q[8] * 83820676633283L) + ((((int128)tmp_q[9] * 101107308146455L) - ((int128)tmp_q[10] * 90037082646790L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 90037082646790L) + ((int128)tmp_q[1] * 10158568668206L) - ((int128)tmp_q[2] * 47147755463972L) - ((int128)tmp_q[3] * 8942258500916L) - ((int128)tmp_q[4] * 10357355016976L) + ((int128)tmp_q[5] * 58927866889179L) - ((int128)tmp_q[6] * 2949931717797L) + ((int128)tmp_q[7] * 81024989924469L) + ((int128)tmp_q[8] * 22439204334350L) - ((int128)tmp_q[9] * 83820676633283L) + ((int128)tmp_q[10] * 505536540732275L);
	tmp_zero[10] = ((int128)tmp_q[0] * 101107308146455L) - ((int128)tmp_q[1] * 90037082646790L) + ((int128)tmp_q[2] * 10158568668206L) - ((int128)tmp_q[3] * 47147755463972L) - ((int128)tmp_q[4] * 8942258500916L) - ((int128)tmp_q[5] * 10357355016976L) + ((int128)tmp_q[6] * 58927866889179L) - ((int128)tmp_q[7] * 2949931717797L) + ((int128)tmp_q[8] * 81024989924469L) + ((int128)tmp_q[9] * 22439204334350L) - ((int128)tmp_q[10] * 83820676633283L);

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
}

