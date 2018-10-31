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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9867783470019540750UL) + ((((uint64_t)op[1] * 3331260327530482918UL) + ((uint64_t)op[2] * 12164009346327202539UL) + ((uint64_t)op[3] * 11444991014004555127UL) + ((uint64_t)op[4] * 14433317654556676061UL) + ((uint64_t)op[5] * 11636605066072523905UL) + ((uint64_t)op[6] * 11275895507608273896UL) + ((uint64_t)op[7] * 11580013113890348735UL) + ((uint64_t)op[8] * 13685204148065345747UL) + ((uint64_t)op[9] * 18014632120131401216UL) + ((uint64_t)op[10] * 3394747536635185377UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 3394747536635185377UL) + ((uint64_t)op[1] * 9867783470019540750UL) + ((((uint64_t)op[2] * 3331260327530482918UL) + ((uint64_t)op[3] * 12164009346327202539UL) + ((uint64_t)op[4] * 11444991014004555127UL) + ((uint64_t)op[5] * 14433317654556676061UL) + ((uint64_t)op[6] * 11636605066072523905UL) + ((uint64_t)op[7] * 11275895507608273896UL) + ((uint64_t)op[8] * 11580013113890348735UL) + ((uint64_t)op[9] * 13685204148065345747UL) + ((uint64_t)op[10] * 18014632120131401216UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 18014632120131401216UL) + ((uint64_t)op[1] * 3394747536635185377UL) + ((uint64_t)op[2] * 9867783470019540750UL) + ((((uint64_t)op[3] * 3331260327530482918UL) + ((uint64_t)op[4] * 12164009346327202539UL) + ((uint64_t)op[5] * 11444991014004555127UL) + ((uint64_t)op[6] * 14433317654556676061UL) + ((uint64_t)op[7] * 11636605066072523905UL) + ((uint64_t)op[8] * 11275895507608273896UL) + ((uint64_t)op[9] * 11580013113890348735UL) + ((uint64_t)op[10] * 13685204148065345747UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 13685204148065345747UL) + ((uint64_t)op[1] * 18014632120131401216UL) + ((uint64_t)op[2] * 3394747536635185377UL) + ((uint64_t)op[3] * 9867783470019540750UL) + ((((uint64_t)op[4] * 3331260327530482918UL) + ((uint64_t)op[5] * 12164009346327202539UL) + ((uint64_t)op[6] * 11444991014004555127UL) + ((uint64_t)op[7] * 14433317654556676061UL) + ((uint64_t)op[8] * 11636605066072523905UL) + ((uint64_t)op[9] * 11275895507608273896UL) + ((uint64_t)op[10] * 11580013113890348735UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 11580013113890348735UL) + ((uint64_t)op[1] * 13685204148065345747UL) + ((uint64_t)op[2] * 18014632120131401216UL) + ((uint64_t)op[3] * 3394747536635185377UL) + ((uint64_t)op[4] * 9867783470019540750UL) + ((((uint64_t)op[5] * 3331260327530482918UL) + ((uint64_t)op[6] * 12164009346327202539UL) + ((uint64_t)op[7] * 11444991014004555127UL) + ((uint64_t)op[8] * 14433317654556676061UL) + ((uint64_t)op[9] * 11636605066072523905UL) + ((uint64_t)op[10] * 11275895507608273896UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 11275895507608273896UL) + ((uint64_t)op[1] * 11580013113890348735UL) + ((uint64_t)op[2] * 13685204148065345747UL) + ((uint64_t)op[3] * 18014632120131401216UL) + ((uint64_t)op[4] * 3394747536635185377UL) + ((uint64_t)op[5] * 9867783470019540750UL) + ((((uint64_t)op[6] * 3331260327530482918UL) + ((uint64_t)op[7] * 12164009346327202539UL) + ((uint64_t)op[8] * 11444991014004555127UL) + ((uint64_t)op[9] * 14433317654556676061UL) + ((uint64_t)op[10] * 11636605066072523905UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 11636605066072523905UL) + ((uint64_t)op[1] * 11275895507608273896UL) + ((uint64_t)op[2] * 11580013113890348735UL) + ((uint64_t)op[3] * 13685204148065345747UL) + ((uint64_t)op[4] * 18014632120131401216UL) + ((uint64_t)op[5] * 3394747536635185377UL) + ((uint64_t)op[6] * 9867783470019540750UL) + ((((uint64_t)op[7] * 3331260327530482918UL) + ((uint64_t)op[8] * 12164009346327202539UL) + ((uint64_t)op[9] * 11444991014004555127UL) + ((uint64_t)op[10] * 14433317654556676061UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 14433317654556676061UL) + ((uint64_t)op[1] * 11636605066072523905UL) + ((uint64_t)op[2] * 11275895507608273896UL) + ((uint64_t)op[3] * 11580013113890348735UL) + ((uint64_t)op[4] * 13685204148065345747UL) + ((uint64_t)op[5] * 18014632120131401216UL) + ((uint64_t)op[6] * 3394747536635185377UL) + ((uint64_t)op[7] * 9867783470019540750UL) + ((((uint64_t)op[8] * 3331260327530482918UL) + ((uint64_t)op[9] * 12164009346327202539UL) + ((uint64_t)op[10] * 11444991014004555127UL)) * 18446744073709551611);
	tmp_q[8] = ((uint64_t)op[0] * 11444991014004555127UL) + ((uint64_t)op[1] * 14433317654556676061UL) + ((uint64_t)op[2] * 11636605066072523905UL) + ((uint64_t)op[3] * 11275895507608273896UL) + ((uint64_t)op[4] * 11580013113890348735UL) + ((uint64_t)op[5] * 13685204148065345747UL) + ((uint64_t)op[6] * 18014632120131401216UL) + ((uint64_t)op[7] * 3394747536635185377UL) + ((uint64_t)op[8] * 9867783470019540750UL) + ((((uint64_t)op[9] * 3331260327530482918UL) + ((uint64_t)op[10] * 12164009346327202539UL)) * 18446744073709551611);
	tmp_q[9] = ((uint64_t)op[0] * 12164009346327202539UL) + ((uint64_t)op[1] * 11444991014004555127UL) + ((uint64_t)op[2] * 14433317654556676061UL) + ((uint64_t)op[3] * 11636605066072523905UL) + ((uint64_t)op[4] * 11275895507608273896UL) + ((uint64_t)op[5] * 11580013113890348735UL) + ((uint64_t)op[6] * 13685204148065345747UL) + ((uint64_t)op[7] * 18014632120131401216UL) + ((uint64_t)op[8] * 3394747536635185377UL) + ((uint64_t)op[9] * 9867783470019540750UL) + ((uint64_t)op[10] * 1790442436057137026UL);
	tmp_q[10] = ((uint64_t)op[0] * 3331260327530482918UL) + ((uint64_t)op[1] * 12164009346327202539UL) + ((uint64_t)op[2] * 11444991014004555127UL) + ((uint64_t)op[3] * 14433317654556676061UL) + ((uint64_t)op[4] * 11636605066072523905UL) + ((uint64_t)op[5] * 11275895507608273896UL) + ((uint64_t)op[6] * 11580013113890348735UL) + ((uint64_t)op[7] * 13685204148065345747UL) + ((uint64_t)op[8] * 18014632120131401216UL) + ((uint64_t)op[9] * 3394747536635185377UL) + ((uint64_t)op[10] * 9867783470019540750UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 22497945635073L) - ((-((int128)tmp_q[1] * 83322240279113L) - ((int128)tmp_q[2] * 19953798831968L) + ((int128)tmp_q[3] * 121810842421665L) - ((int128)tmp_q[4] * 3648999536992L) - ((int128)tmp_q[5] * 47754366921725L) + ((int128)tmp_q[6] * 66871746631876L) + ((int128)tmp_q[7] * 109955879918138L) - ((int128)tmp_q[8] * 62215625099441L) - ((int128)tmp_q[9] * 32508318949714L) + ((int128)tmp_q[10] * 28826195643716L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 28826195643716L) + ((int128)tmp_q[1] * 22497945635073L) - ((-((int128)tmp_q[2] * 83322240279113L) - ((int128)tmp_q[3] * 19953798831968L) + ((int128)tmp_q[4] * 121810842421665L) - ((int128)tmp_q[5] * 3648999536992L) - ((int128)tmp_q[6] * 47754366921725L) + ((int128)tmp_q[7] * 66871746631876L) + ((int128)tmp_q[8] * 109955879918138L) - ((int128)tmp_q[9] * 62215625099441L) - ((int128)tmp_q[10] * 32508318949714L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 32508318949714L) + ((int128)tmp_q[1] * 28826195643716L) + ((int128)tmp_q[2] * 22497945635073L) - ((-((int128)tmp_q[3] * 83322240279113L) - ((int128)tmp_q[4] * 19953798831968L) + ((int128)tmp_q[5] * 121810842421665L) - ((int128)tmp_q[6] * 3648999536992L) - ((int128)tmp_q[7] * 47754366921725L) + ((int128)tmp_q[8] * 66871746631876L) + ((int128)tmp_q[9] * 109955879918138L) - ((int128)tmp_q[10] * 62215625099441L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 62215625099441L) - ((int128)tmp_q[1] * 32508318949714L) + ((int128)tmp_q[2] * 28826195643716L) + ((int128)tmp_q[3] * 22497945635073L) - ((-((int128)tmp_q[4] * 83322240279113L) - ((int128)tmp_q[5] * 19953798831968L) + ((int128)tmp_q[6] * 121810842421665L) - ((int128)tmp_q[7] * 3648999536992L) - ((int128)tmp_q[8] * 47754366921725L) + ((int128)tmp_q[9] * 66871746631876L) + ((int128)tmp_q[10] * 109955879918138L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 109955879918138L) - ((int128)tmp_q[1] * 62215625099441L) - ((int128)tmp_q[2] * 32508318949714L) + ((int128)tmp_q[3] * 28826195643716L) + ((int128)tmp_q[4] * 22497945635073L) - ((-((int128)tmp_q[5] * 83322240279113L) - ((int128)tmp_q[6] * 19953798831968L) + ((int128)tmp_q[7] * 121810842421665L) - ((int128)tmp_q[8] * 3648999536992L) - ((int128)tmp_q[9] * 47754366921725L) + ((int128)tmp_q[10] * 66871746631876L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 66871746631876L) + ((int128)tmp_q[1] * 109955879918138L) - ((int128)tmp_q[2] * 62215625099441L) - ((int128)tmp_q[3] * 32508318949714L) + ((int128)tmp_q[4] * 28826195643716L) + ((int128)tmp_q[5] * 22497945635073L) - ((-((int128)tmp_q[6] * 83322240279113L) - ((int128)tmp_q[7] * 19953798831968L) + ((int128)tmp_q[8] * 121810842421665L) - ((int128)tmp_q[9] * 3648999536992L) - ((int128)tmp_q[10] * 47754366921725L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 47754366921725L) + ((int128)tmp_q[1] * 66871746631876L) + ((int128)tmp_q[2] * 109955879918138L) - ((int128)tmp_q[3] * 62215625099441L) - ((int128)tmp_q[4] * 32508318949714L) + ((int128)tmp_q[5] * 28826195643716L) + ((int128)tmp_q[6] * 22497945635073L) - ((-((int128)tmp_q[7] * 83322240279113L) - ((int128)tmp_q[8] * 19953798831968L) + ((int128)tmp_q[9] * 121810842421665L) - ((int128)tmp_q[10] * 3648999536992L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 3648999536992L) - ((int128)tmp_q[1] * 47754366921725L) + ((int128)tmp_q[2] * 66871746631876L) + ((int128)tmp_q[3] * 109955879918138L) - ((int128)tmp_q[4] * 62215625099441L) - ((int128)tmp_q[5] * 32508318949714L) + ((int128)tmp_q[6] * 28826195643716L) + ((int128)tmp_q[7] * 22497945635073L) - ((-((int128)tmp_q[8] * 83322240279113L) - ((int128)tmp_q[9] * 19953798831968L) + ((int128)tmp_q[10] * 121810842421665L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 121810842421665L) - ((int128)tmp_q[1] * 3648999536992L) - ((int128)tmp_q[2] * 47754366921725L) + ((int128)tmp_q[3] * 66871746631876L) + ((int128)tmp_q[4] * 109955879918138L) - ((int128)tmp_q[5] * 62215625099441L) - ((int128)tmp_q[6] * 32508318949714L) + ((int128)tmp_q[7] * 28826195643716L) + ((int128)tmp_q[8] * 22497945635073L) - ((-((int128)tmp_q[9] * 83322240279113L) - ((int128)tmp_q[10] * 19953798831968L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 19953798831968L) + ((int128)tmp_q[1] * 121810842421665L) - ((int128)tmp_q[2] * 3648999536992L) - ((int128)tmp_q[3] * 47754366921725L) + ((int128)tmp_q[4] * 66871746631876L) + ((int128)tmp_q[5] * 109955879918138L) - ((int128)tmp_q[6] * 62215625099441L) - ((int128)tmp_q[7] * 32508318949714L) + ((int128)tmp_q[8] * 28826195643716L) + ((int128)tmp_q[9] * 22497945635073L) + ((int128)tmp_q[10] * 416611201395565L);
	tmp_zero[10] = -((int128)tmp_q[0] * 83322240279113L) - ((int128)tmp_q[1] * 19953798831968L) + ((int128)tmp_q[2] * 121810842421665L) - ((int128)tmp_q[3] * 3648999536992L) - ((int128)tmp_q[4] * 47754366921725L) + ((int128)tmp_q[5] * 66871746631876L) + ((int128)tmp_q[6] * 109955879918138L) - ((int128)tmp_q[7] * 62215625099441L) - ((int128)tmp_q[8] * 32508318949714L) + ((int128)tmp_q[9] * 28826195643716L) + ((int128)tmp_q[10] * 22497945635073L);

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

