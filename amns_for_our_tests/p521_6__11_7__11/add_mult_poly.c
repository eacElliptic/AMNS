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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9460364029530569355UL) + ((((uint64_t)op[1] * 5999975981379495800UL) + ((uint64_t)op[2] * 545785716931920108UL) + ((uint64_t)op[3] * 2625840924675445375UL) + ((uint64_t)op[4] * 442480745143131249UL) + ((uint64_t)op[5] * 5844803558561112671UL) + ((uint64_t)op[6] * 8950042792753012103UL) + ((uint64_t)op[7] * 8329227487383619288UL) + ((uint64_t)op[8] * 1684258504140901019UL) + ((uint64_t)op[9] * 12184778758109675827UL) + ((uint64_t)op[10] * 14531690733357639896UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 14531690733357639896UL) + ((uint64_t)op[1] * 9460364029530569355UL) + ((((uint64_t)op[2] * 5999975981379495800UL) + ((uint64_t)op[3] * 545785716931920108UL) + ((uint64_t)op[4] * 2625840924675445375UL) + ((uint64_t)op[5] * 442480745143131249UL) + ((uint64_t)op[6] * 5844803558561112671UL) + ((uint64_t)op[7] * 8950042792753012103UL) + ((uint64_t)op[8] * 8329227487383619288UL) + ((uint64_t)op[9] * 1684258504140901019UL) + ((uint64_t)op[10] * 12184778758109675827UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 12184778758109675827UL) + ((uint64_t)op[1] * 14531690733357639896UL) + ((uint64_t)op[2] * 9460364029530569355UL) + ((((uint64_t)op[3] * 5999975981379495800UL) + ((uint64_t)op[4] * 545785716931920108UL) + ((uint64_t)op[5] * 2625840924675445375UL) + ((uint64_t)op[6] * 442480745143131249UL) + ((uint64_t)op[7] * 5844803558561112671UL) + ((uint64_t)op[8] * 8950042792753012103UL) + ((uint64_t)op[9] * 8329227487383619288UL) + ((uint64_t)op[10] * 1684258504140901019UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 1684258504140901019UL) + ((uint64_t)op[1] * 12184778758109675827UL) + ((uint64_t)op[2] * 14531690733357639896UL) + ((uint64_t)op[3] * 9460364029530569355UL) + ((((uint64_t)op[4] * 5999975981379495800UL) + ((uint64_t)op[5] * 545785716931920108UL) + ((uint64_t)op[6] * 2625840924675445375UL) + ((uint64_t)op[7] * 442480745143131249UL) + ((uint64_t)op[8] * 5844803558561112671UL) + ((uint64_t)op[9] * 8950042792753012103UL) + ((uint64_t)op[10] * 8329227487383619288UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 8329227487383619288UL) + ((uint64_t)op[1] * 1684258504140901019UL) + ((uint64_t)op[2] * 12184778758109675827UL) + ((uint64_t)op[3] * 14531690733357639896UL) + ((uint64_t)op[4] * 9460364029530569355UL) + ((((uint64_t)op[5] * 5999975981379495800UL) + ((uint64_t)op[6] * 545785716931920108UL) + ((uint64_t)op[7] * 2625840924675445375UL) + ((uint64_t)op[8] * 442480745143131249UL) + ((uint64_t)op[9] * 5844803558561112671UL) + ((uint64_t)op[10] * 8950042792753012103UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 8950042792753012103UL) + ((uint64_t)op[1] * 8329227487383619288UL) + ((uint64_t)op[2] * 1684258504140901019UL) + ((uint64_t)op[3] * 12184778758109675827UL) + ((uint64_t)op[4] * 14531690733357639896UL) + ((uint64_t)op[5] * 9460364029530569355UL) + ((((uint64_t)op[6] * 5999975981379495800UL) + ((uint64_t)op[7] * 545785716931920108UL) + ((uint64_t)op[8] * 2625840924675445375UL) + ((uint64_t)op[9] * 442480745143131249UL) + ((uint64_t)op[10] * 5844803558561112671UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 5844803558561112671UL) + ((uint64_t)op[1] * 8950042792753012103UL) + ((uint64_t)op[2] * 8329227487383619288UL) + ((uint64_t)op[3] * 1684258504140901019UL) + ((uint64_t)op[4] * 12184778758109675827UL) + ((uint64_t)op[5] * 14531690733357639896UL) + ((uint64_t)op[6] * 9460364029530569355UL) + ((((uint64_t)op[7] * 5999975981379495800UL) + ((uint64_t)op[8] * 545785716931920108UL) + ((uint64_t)op[9] * 2625840924675445375UL) + ((uint64_t)op[10] * 442480745143131249UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 442480745143131249UL) + ((uint64_t)op[1] * 5844803558561112671UL) + ((uint64_t)op[2] * 8950042792753012103UL) + ((uint64_t)op[3] * 8329227487383619288UL) + ((uint64_t)op[4] * 1684258504140901019UL) + ((uint64_t)op[5] * 12184778758109675827UL) + ((uint64_t)op[6] * 14531690733357639896UL) + ((uint64_t)op[7] * 9460364029530569355UL) + ((((uint64_t)op[8] * 5999975981379495800UL) + ((uint64_t)op[9] * 545785716931920108UL) + ((uint64_t)op[10] * 2625840924675445375UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 2625840924675445375UL) + ((uint64_t)op[1] * 442480745143131249UL) + ((uint64_t)op[2] * 5844803558561112671UL) + ((uint64_t)op[3] * 8950042792753012103UL) + ((uint64_t)op[4] * 8329227487383619288UL) + ((uint64_t)op[5] * 1684258504140901019UL) + ((uint64_t)op[6] * 12184778758109675827UL) + ((uint64_t)op[7] * 14531690733357639896UL) + ((uint64_t)op[8] * 9460364029530569355UL) + ((((uint64_t)op[9] * 5999975981379495800UL) + ((uint64_t)op[10] * 545785716931920108UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 545785716931920108UL) + ((uint64_t)op[1] * 2625840924675445375UL) + ((uint64_t)op[2] * 442480745143131249UL) + ((uint64_t)op[3] * 5844803558561112671UL) + ((uint64_t)op[4] * 8950042792753012103UL) + ((uint64_t)op[5] * 8329227487383619288UL) + ((uint64_t)op[6] * 1684258504140901019UL) + ((uint64_t)op[7] * 12184778758109675827UL) + ((uint64_t)op[8] * 14531690733357639896UL) + ((uint64_t)op[9] * 9460364029530569355UL) + ((uint64_t)op[10] * 5106343722237367368UL);
	tmp_q[10] = ((uint64_t)op[0] * 5999975981379495800UL) + ((uint64_t)op[1] * 545785716931920108UL) + ((uint64_t)op[2] * 2625840924675445375UL) + ((uint64_t)op[3] * 442480745143131249UL) + ((uint64_t)op[4] * 5844803558561112671UL) + ((uint64_t)op[5] * 8950042792753012103UL) + ((uint64_t)op[6] * 8329227487383619288UL) + ((uint64_t)op[7] * 1684258504140901019UL) + ((uint64_t)op[8] * 12184778758109675827UL) + ((uint64_t)op[9] * 14531690733357639896UL) + ((uint64_t)op[10] * 9460364029530569355UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 56128548650945L) + ((((int128)tmp_q[1] * 51992880976396L) + ((int128)tmp_q[2] * 32196121540155L) - ((int128)tmp_q[3] * 61926276514266L) + ((int128)tmp_q[4] * 79287654246535L) - ((int128)tmp_q[5] * 43902669172566L) + ((int128)tmp_q[6] * 42720791630002L) - ((int128)tmp_q[7] * 5980688851971L) + ((int128)tmp_q[8] * 60793253649000L) - ((int128)tmp_q[9] * 13053168184176L) - ((int128)tmp_q[10] * 84082629914757L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 84082629914757L) + ((int128)tmp_q[1] * 56128548650945L) + ((((int128)tmp_q[2] * 51992880976396L) + ((int128)tmp_q[3] * 32196121540155L) - ((int128)tmp_q[4] * 61926276514266L) + ((int128)tmp_q[5] * 79287654246535L) - ((int128)tmp_q[6] * 43902669172566L) + ((int128)tmp_q[7] * 42720791630002L) - ((int128)tmp_q[8] * 5980688851971L) + ((int128)tmp_q[9] * 60793253649000L) - ((int128)tmp_q[10] * 13053168184176L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 13053168184176L) - ((int128)tmp_q[1] * 84082629914757L) + ((int128)tmp_q[2] * 56128548650945L) + ((((int128)tmp_q[3] * 51992880976396L) + ((int128)tmp_q[4] * 32196121540155L) - ((int128)tmp_q[5] * 61926276514266L) + ((int128)tmp_q[6] * 79287654246535L) - ((int128)tmp_q[7] * 43902669172566L) + ((int128)tmp_q[8] * 42720791630002L) - ((int128)tmp_q[9] * 5980688851971L) + ((int128)tmp_q[10] * 60793253649000L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 60793253649000L) - ((int128)tmp_q[1] * 13053168184176L) - ((int128)tmp_q[2] * 84082629914757L) + ((int128)tmp_q[3] * 56128548650945L) + ((((int128)tmp_q[4] * 51992880976396L) + ((int128)tmp_q[5] * 32196121540155L) - ((int128)tmp_q[6] * 61926276514266L) + ((int128)tmp_q[7] * 79287654246535L) - ((int128)tmp_q[8] * 43902669172566L) + ((int128)tmp_q[9] * 42720791630002L) - ((int128)tmp_q[10] * 5980688851971L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 5980688851971L) + ((int128)tmp_q[1] * 60793253649000L) - ((int128)tmp_q[2] * 13053168184176L) - ((int128)tmp_q[3] * 84082629914757L) + ((int128)tmp_q[4] * 56128548650945L) + ((((int128)tmp_q[5] * 51992880976396L) + ((int128)tmp_q[6] * 32196121540155L) - ((int128)tmp_q[7] * 61926276514266L) + ((int128)tmp_q[8] * 79287654246535L) - ((int128)tmp_q[9] * 43902669172566L) + ((int128)tmp_q[10] * 42720791630002L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 42720791630002L) - ((int128)tmp_q[1] * 5980688851971L) + ((int128)tmp_q[2] * 60793253649000L) - ((int128)tmp_q[3] * 13053168184176L) - ((int128)tmp_q[4] * 84082629914757L) + ((int128)tmp_q[5] * 56128548650945L) + ((((int128)tmp_q[6] * 51992880976396L) + ((int128)tmp_q[7] * 32196121540155L) - ((int128)tmp_q[8] * 61926276514266L) + ((int128)tmp_q[9] * 79287654246535L) - ((int128)tmp_q[10] * 43902669172566L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 43902669172566L) + ((int128)tmp_q[1] * 42720791630002L) - ((int128)tmp_q[2] * 5980688851971L) + ((int128)tmp_q[3] * 60793253649000L) - ((int128)tmp_q[4] * 13053168184176L) - ((int128)tmp_q[5] * 84082629914757L) + ((int128)tmp_q[6] * 56128548650945L) + ((((int128)tmp_q[7] * 51992880976396L) + ((int128)tmp_q[8] * 32196121540155L) - ((int128)tmp_q[9] * 61926276514266L) + ((int128)tmp_q[10] * 79287654246535L)) * 7);
	tmp_zero[7] = ((int128)tmp_q[0] * 79287654246535L) - ((int128)tmp_q[1] * 43902669172566L) + ((int128)tmp_q[2] * 42720791630002L) - ((int128)tmp_q[3] * 5980688851971L) + ((int128)tmp_q[4] * 60793253649000L) - ((int128)tmp_q[5] * 13053168184176L) - ((int128)tmp_q[6] * 84082629914757L) + ((int128)tmp_q[7] * 56128548650945L) + ((((int128)tmp_q[8] * 51992880976396L) + ((int128)tmp_q[9] * 32196121540155L) - ((int128)tmp_q[10] * 61926276514266L)) * 7);
	tmp_zero[8] = -((int128)tmp_q[0] * 61926276514266L) + ((int128)tmp_q[1] * 79287654246535L) - ((int128)tmp_q[2] * 43902669172566L) + ((int128)tmp_q[3] * 42720791630002L) - ((int128)tmp_q[4] * 5980688851971L) + ((int128)tmp_q[5] * 60793253649000L) - ((int128)tmp_q[6] * 13053168184176L) - ((int128)tmp_q[7] * 84082629914757L) + ((int128)tmp_q[8] * 56128548650945L) + ((((int128)tmp_q[9] * 51992880976396L) + ((int128)tmp_q[10] * 32196121540155L)) * 7);
	tmp_zero[9] = ((int128)tmp_q[0] * 32196121540155L) - ((int128)tmp_q[1] * 61926276514266L) + ((int128)tmp_q[2] * 79287654246535L) - ((int128)tmp_q[3] * 43902669172566L) + ((int128)tmp_q[4] * 42720791630002L) - ((int128)tmp_q[5] * 5980688851971L) + ((int128)tmp_q[6] * 60793253649000L) - ((int128)tmp_q[7] * 13053168184176L) - ((int128)tmp_q[8] * 84082629914757L) + ((int128)tmp_q[9] * 56128548650945L) + ((int128)tmp_q[10] * 363950166834772L);
	tmp_zero[10] = ((int128)tmp_q[0] * 51992880976396L) + ((int128)tmp_q[1] * 32196121540155L) - ((int128)tmp_q[2] * 61926276514266L) + ((int128)tmp_q[3] * 79287654246535L) - ((int128)tmp_q[4] * 43902669172566L) + ((int128)tmp_q[5] * 42720791630002L) - ((int128)tmp_q[6] * 5980688851971L) + ((int128)tmp_q[7] * 60793253649000L) - ((int128)tmp_q[8] * 13053168184176L) - ((int128)tmp_q[9] * 84082629914757L) + ((int128)tmp_q[10] * 56128548650945L);

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

