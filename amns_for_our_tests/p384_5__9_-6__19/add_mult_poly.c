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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2447921189243322703UL) + ((((uint64_t)op[1] * 14559523510741963297UL) + ((uint64_t)op[2] * 18070510644082439544UL) + ((uint64_t)op[3] * 11813410489692549711UL) + ((uint64_t)op[4] * 14365027723236954798UL) + ((uint64_t)op[5] * 10826264243531247909UL) + ((uint64_t)op[6] * 4008428896803197529UL) + ((uint64_t)op[7] * 7915979850047004493UL) + ((uint64_t)op[8] * 1009507632387458537UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 1009507632387458537UL) + ((uint64_t)op[1] * 2447921189243322703UL) + ((((uint64_t)op[2] * 14559523510741963297UL) + ((uint64_t)op[3] * 18070510644082439544UL) + ((uint64_t)op[4] * 11813410489692549711UL) + ((uint64_t)op[5] * 14365027723236954798UL) + ((uint64_t)op[6] * 10826264243531247909UL) + ((uint64_t)op[7] * 4008428896803197529UL) + ((uint64_t)op[8] * 7915979850047004493UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 7915979850047004493UL) + ((uint64_t)op[1] * 1009507632387458537UL) + ((uint64_t)op[2] * 2447921189243322703UL) + ((((uint64_t)op[3] * 14559523510741963297UL) + ((uint64_t)op[4] * 18070510644082439544UL) + ((uint64_t)op[5] * 11813410489692549711UL) + ((uint64_t)op[6] * 14365027723236954798UL) + ((uint64_t)op[7] * 10826264243531247909UL) + ((uint64_t)op[8] * 4008428896803197529UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 4008428896803197529UL) + ((uint64_t)op[1] * 7915979850047004493UL) + ((uint64_t)op[2] * 1009507632387458537UL) + ((uint64_t)op[3] * 2447921189243322703UL) + ((((uint64_t)op[4] * 14559523510741963297UL) + ((uint64_t)op[5] * 18070510644082439544UL) + ((uint64_t)op[6] * 11813410489692549711UL) + ((uint64_t)op[7] * 14365027723236954798UL) + ((uint64_t)op[8] * 10826264243531247909UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 10826264243531247909UL) + ((uint64_t)op[1] * 4008428896803197529UL) + ((uint64_t)op[2] * 7915979850047004493UL) + ((uint64_t)op[3] * 1009507632387458537UL) + ((uint64_t)op[4] * 2447921189243322703UL) + ((((uint64_t)op[5] * 14559523510741963297UL) + ((uint64_t)op[6] * 18070510644082439544UL) + ((uint64_t)op[7] * 11813410489692549711UL) + ((uint64_t)op[8] * 14365027723236954798UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 14365027723236954798UL) + ((uint64_t)op[1] * 10826264243531247909UL) + ((uint64_t)op[2] * 4008428896803197529UL) + ((uint64_t)op[3] * 7915979850047004493UL) + ((uint64_t)op[4] * 1009507632387458537UL) + ((uint64_t)op[5] * 2447921189243322703UL) + ((((uint64_t)op[6] * 14559523510741963297UL) + ((uint64_t)op[7] * 18070510644082439544UL) + ((uint64_t)op[8] * 11813410489692549711UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 11813410489692549711UL) + ((uint64_t)op[1] * 14365027723236954798UL) + ((uint64_t)op[2] * 10826264243531247909UL) + ((uint64_t)op[3] * 4008428896803197529UL) + ((uint64_t)op[4] * 7915979850047004493UL) + ((uint64_t)op[5] * 1009507632387458537UL) + ((uint64_t)op[6] * 2447921189243322703UL) + ((((uint64_t)op[7] * 14559523510741963297UL) + ((uint64_t)op[8] * 18070510644082439544UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 18070510644082439544UL) + ((uint64_t)op[1] * 11813410489692549711UL) + ((uint64_t)op[2] * 14365027723236954798UL) + ((uint64_t)op[3] * 10826264243531247909UL) + ((uint64_t)op[4] * 4008428896803197529UL) + ((uint64_t)op[5] * 7915979850047004493UL) + ((uint64_t)op[6] * 1009507632387458537UL) + ((uint64_t)op[7] * 2447921189243322703UL) + ((uint64_t)op[8] * 4876579304095978298UL);
	tmp_q[8] = ((uint64_t)op[0] * 14559523510741963297UL) + ((uint64_t)op[1] * 18070510644082439544UL) + ((uint64_t)op[2] * 11813410489692549711UL) + ((uint64_t)op[3] * 14365027723236954798UL) + ((uint64_t)op[4] * 10826264243531247909UL) + ((uint64_t)op[5] * 4008428896803197529UL) + ((uint64_t)op[6] * 7915979850047004493UL) + ((uint64_t)op[7] * 1009507632387458537UL) + ((uint64_t)op[8] * 2447921189243322703UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1407422805009L) - ((((int128)tmp_q[1] * 1462692759172L) - ((int128)tmp_q[2] * 1023984117006L) - ((int128)tmp_q[3] * 2268548477880L) + ((int128)tmp_q[4] * 153313605399L) - ((int128)tmp_q[5] * 3136469648870L) + ((int128)tmp_q[6] * 2570131304306L) + ((int128)tmp_q[7] * 291614487148L) - ((int128)tmp_q[8] * 2724219735811L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 2724219735811L) + ((int128)tmp_q[1] * 1407422805009L) - ((((int128)tmp_q[2] * 1462692759172L) - ((int128)tmp_q[3] * 1023984117006L) - ((int128)tmp_q[4] * 2268548477880L) + ((int128)tmp_q[5] * 153313605399L) - ((int128)tmp_q[6] * 3136469648870L) + ((int128)tmp_q[7] * 2570131304306L) + ((int128)tmp_q[8] * 291614487148L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 291614487148L) - ((int128)tmp_q[1] * 2724219735811L) + ((int128)tmp_q[2] * 1407422805009L) - ((((int128)tmp_q[3] * 1462692759172L) - ((int128)tmp_q[4] * 1023984117006L) - ((int128)tmp_q[5] * 2268548477880L) + ((int128)tmp_q[6] * 153313605399L) - ((int128)tmp_q[7] * 3136469648870L) + ((int128)tmp_q[8] * 2570131304306L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 2570131304306L) + ((int128)tmp_q[1] * 291614487148L) - ((int128)tmp_q[2] * 2724219735811L) + ((int128)tmp_q[3] * 1407422805009L) - ((((int128)tmp_q[4] * 1462692759172L) - ((int128)tmp_q[5] * 1023984117006L) - ((int128)tmp_q[6] * 2268548477880L) + ((int128)tmp_q[7] * 153313605399L) - ((int128)tmp_q[8] * 3136469648870L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 3136469648870L) + ((int128)tmp_q[1] * 2570131304306L) + ((int128)tmp_q[2] * 291614487148L) - ((int128)tmp_q[3] * 2724219735811L) + ((int128)tmp_q[4] * 1407422805009L) - ((((int128)tmp_q[5] * 1462692759172L) - ((int128)tmp_q[6] * 1023984117006L) - ((int128)tmp_q[7] * 2268548477880L) + ((int128)tmp_q[8] * 153313605399L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 153313605399L) - ((int128)tmp_q[1] * 3136469648870L) + ((int128)tmp_q[2] * 2570131304306L) + ((int128)tmp_q[3] * 291614487148L) - ((int128)tmp_q[4] * 2724219735811L) + ((int128)tmp_q[5] * 1407422805009L) - ((((int128)tmp_q[6] * 1462692759172L) - ((int128)tmp_q[7] * 1023984117006L) - ((int128)tmp_q[8] * 2268548477880L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 2268548477880L) + ((int128)tmp_q[1] * 153313605399L) - ((int128)tmp_q[2] * 3136469648870L) + ((int128)tmp_q[3] * 2570131304306L) + ((int128)tmp_q[4] * 291614487148L) - ((int128)tmp_q[5] * 2724219735811L) + ((int128)tmp_q[6] * 1407422805009L) - ((((int128)tmp_q[7] * 1462692759172L) - ((int128)tmp_q[8] * 1023984117006L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 1023984117006L) - ((int128)tmp_q[1] * 2268548477880L) + ((int128)tmp_q[2] * 153313605399L) - ((int128)tmp_q[3] * 3136469648870L) + ((int128)tmp_q[4] * 2570131304306L) + ((int128)tmp_q[5] * 291614487148L) - ((int128)tmp_q[6] * 2724219735811L) + ((int128)tmp_q[7] * 1407422805009L) - ((int128)tmp_q[8] * 8776156555032L);
	tmp_zero[8] = ((int128)tmp_q[0] * 1462692759172L) - ((int128)tmp_q[1] * 1023984117006L) - ((int128)tmp_q[2] * 2268548477880L) + ((int128)tmp_q[3] * 153313605399L) - ((int128)tmp_q[4] * 3136469648870L) + ((int128)tmp_q[5] * 2570131304306L) + ((int128)tmp_q[6] * 291614487148L) - ((int128)tmp_q[7] * 2724219735811L) + ((int128)tmp_q[8] * 1407422805009L);

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
}

