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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[9] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1352797712770022256UL) + ((((uint64_t)op[1] * 10699012603819610346UL) + ((uint64_t)op[2] * 5849344524544566368UL) + ((uint64_t)op[3] * 9844142778467427563UL) + ((uint64_t)op[4] * 10710780077722662582UL) + ((uint64_t)op[5] * 6135646568188505687UL) + ((uint64_t)op[6] * 8258249587785445469UL) + ((uint64_t)op[7] * 2835658623693396649UL) + ((uint64_t)op[8] * 5154123103763065209UL) + ((uint64_t)op[9] * 13512480557635958030UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 13512480557635958030UL) + ((uint64_t)op[1] * 1352797712770022256UL) + ((((uint64_t)op[2] * 10699012603819610346UL) + ((uint64_t)op[3] * 5849344524544566368UL) + ((uint64_t)op[4] * 9844142778467427563UL) + ((uint64_t)op[5] * 10710780077722662582UL) + ((uint64_t)op[6] * 6135646568188505687UL) + ((uint64_t)op[7] * 8258249587785445469UL) + ((uint64_t)op[8] * 2835658623693396649UL) + ((uint64_t)op[9] * 5154123103763065209UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 5154123103763065209UL) + ((uint64_t)op[1] * 13512480557635958030UL) + ((uint64_t)op[2] * 1352797712770022256UL) + ((((uint64_t)op[3] * 10699012603819610346UL) + ((uint64_t)op[4] * 5849344524544566368UL) + ((uint64_t)op[5] * 9844142778467427563UL) + ((uint64_t)op[6] * 10710780077722662582UL) + ((uint64_t)op[7] * 6135646568188505687UL) + ((uint64_t)op[8] * 8258249587785445469UL) + ((uint64_t)op[9] * 2835658623693396649UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 2835658623693396649UL) + ((uint64_t)op[1] * 5154123103763065209UL) + ((uint64_t)op[2] * 13512480557635958030UL) + ((uint64_t)op[3] * 1352797712770022256UL) + ((((uint64_t)op[4] * 10699012603819610346UL) + ((uint64_t)op[5] * 5849344524544566368UL) + ((uint64_t)op[6] * 9844142778467427563UL) + ((uint64_t)op[7] * 10710780077722662582UL) + ((uint64_t)op[8] * 6135646568188505687UL) + ((uint64_t)op[9] * 8258249587785445469UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 8258249587785445469UL) + ((uint64_t)op[1] * 2835658623693396649UL) + ((uint64_t)op[2] * 5154123103763065209UL) + ((uint64_t)op[3] * 13512480557635958030UL) + ((uint64_t)op[4] * 1352797712770022256UL) + ((((uint64_t)op[5] * 10699012603819610346UL) + ((uint64_t)op[6] * 5849344524544566368UL) + ((uint64_t)op[7] * 9844142778467427563UL) + ((uint64_t)op[8] * 10710780077722662582UL) + ((uint64_t)op[9] * 6135646568188505687UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 6135646568188505687UL) + ((uint64_t)op[1] * 8258249587785445469UL) + ((uint64_t)op[2] * 2835658623693396649UL) + ((uint64_t)op[3] * 5154123103763065209UL) + ((uint64_t)op[4] * 13512480557635958030UL) + ((uint64_t)op[5] * 1352797712770022256UL) + ((((uint64_t)op[6] * 10699012603819610346UL) + ((uint64_t)op[7] * 5849344524544566368UL) + ((uint64_t)op[8] * 9844142778467427563UL) + ((uint64_t)op[9] * 10710780077722662582UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 10710780077722662582UL) + ((uint64_t)op[1] * 6135646568188505687UL) + ((uint64_t)op[2] * 8258249587785445469UL) + ((uint64_t)op[3] * 2835658623693396649UL) + ((uint64_t)op[4] * 5154123103763065209UL) + ((uint64_t)op[5] * 13512480557635958030UL) + ((uint64_t)op[6] * 1352797712770022256UL) + ((((uint64_t)op[7] * 10699012603819610346UL) + ((uint64_t)op[8] * 5849344524544566368UL) + ((uint64_t)op[9] * 9844142778467427563UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 9844142778467427563UL) + ((uint64_t)op[1] * 10710780077722662582UL) + ((uint64_t)op[2] * 6135646568188505687UL) + ((uint64_t)op[3] * 8258249587785445469UL) + ((uint64_t)op[4] * 2835658623693396649UL) + ((uint64_t)op[5] * 5154123103763065209UL) + ((uint64_t)op[6] * 13512480557635958030UL) + ((uint64_t)op[7] * 1352797712770022256UL) + ((((uint64_t)op[8] * 10699012603819610346UL) + ((uint64_t)op[9] * 5849344524544566368UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 5849344524544566368UL) + ((uint64_t)op[1] * 9844142778467427563UL) + ((uint64_t)op[2] * 10710780077722662582UL) + ((uint64_t)op[3] * 6135646568188505687UL) + ((uint64_t)op[4] * 8258249587785445469UL) + ((uint64_t)op[5] * 2835658623693396649UL) + ((uint64_t)op[6] * 5154123103763065209UL) + ((uint64_t)op[7] * 13512480557635958030UL) + ((uint64_t)op[8] * 1352797712770022256UL) + ((uint64_t)op[9] * 13650293737749279422UL);
	tmp_q[9] = ((uint64_t)op[0] * 10699012603819610346UL) + ((uint64_t)op[1] * 5849344524544566368UL) + ((uint64_t)op[2] * 9844142778467427563UL) + ((uint64_t)op[3] * 10710780077722662582UL) + ((uint64_t)op[4] * 6135646568188505687UL) + ((uint64_t)op[5] * 8258249587785445469UL) + ((uint64_t)op[6] * 2835658623693396649UL) + ((uint64_t)op[7] * 5154123103763065209UL) + ((uint64_t)op[8] * 13512480557635958030UL) + ((uint64_t)op[9] * 1352797712770022256UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 205818738405668L) + ((((int128)tmp_q[1] * 531005014050630L) - ((int128)tmp_q[2] * 1269139445004091L) - ((int128)tmp_q[3] * 2230351078777709L) + ((int128)tmp_q[4] * 413756925522799L) - ((int128)tmp_q[5] * 1570316060621338L) - ((int128)tmp_q[6] * 2976118761151787L) + ((int128)tmp_q[7] * 1415810708167450L) + ((int128)tmp_q[8] * 501868486807353L) - ((int128)tmp_q[9] * 867413771345704L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 867413771345704L) - ((int128)tmp_q[1] * 205818738405668L) + ((((int128)tmp_q[2] * 531005014050630L) - ((int128)tmp_q[3] * 1269139445004091L) - ((int128)tmp_q[4] * 2230351078777709L) + ((int128)tmp_q[5] * 413756925522799L) - ((int128)tmp_q[6] * 1570316060621338L) - ((int128)tmp_q[7] * 2976118761151787L) + ((int128)tmp_q[8] * 1415810708167450L) + ((int128)tmp_q[9] * 501868486807353L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 501868486807353L) - ((int128)tmp_q[1] * 867413771345704L) - ((int128)tmp_q[2] * 205818738405668L) + ((((int128)tmp_q[3] * 531005014050630L) - ((int128)tmp_q[4] * 1269139445004091L) - ((int128)tmp_q[5] * 2230351078777709L) + ((int128)tmp_q[6] * 413756925522799L) - ((int128)tmp_q[7] * 1570316060621338L) - ((int128)tmp_q[8] * 2976118761151787L) + ((int128)tmp_q[9] * 1415810708167450L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1415810708167450L) + ((int128)tmp_q[1] * 501868486807353L) - ((int128)tmp_q[2] * 867413771345704L) - ((int128)tmp_q[3] * 205818738405668L) + ((((int128)tmp_q[4] * 531005014050630L) - ((int128)tmp_q[5] * 1269139445004091L) - ((int128)tmp_q[6] * 2230351078777709L) + ((int128)tmp_q[7] * 413756925522799L) - ((int128)tmp_q[8] * 1570316060621338L) - ((int128)tmp_q[9] * 2976118761151787L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 2976118761151787L) + ((int128)tmp_q[1] * 1415810708167450L) + ((int128)tmp_q[2] * 501868486807353L) - ((int128)tmp_q[3] * 867413771345704L) - ((int128)tmp_q[4] * 205818738405668L) + ((((int128)tmp_q[5] * 531005014050630L) - ((int128)tmp_q[6] * 1269139445004091L) - ((int128)tmp_q[7] * 2230351078777709L) + ((int128)tmp_q[8] * 413756925522799L) - ((int128)tmp_q[9] * 1570316060621338L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 1570316060621338L) - ((int128)tmp_q[1] * 2976118761151787L) + ((int128)tmp_q[2] * 1415810708167450L) + ((int128)tmp_q[3] * 501868486807353L) - ((int128)tmp_q[4] * 867413771345704L) - ((int128)tmp_q[5] * 205818738405668L) + ((((int128)tmp_q[6] * 531005014050630L) - ((int128)tmp_q[7] * 1269139445004091L) - ((int128)tmp_q[8] * 2230351078777709L) + ((int128)tmp_q[9] * 413756925522799L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 413756925522799L) - ((int128)tmp_q[1] * 1570316060621338L) - ((int128)tmp_q[2] * 2976118761151787L) + ((int128)tmp_q[3] * 1415810708167450L) + ((int128)tmp_q[4] * 501868486807353L) - ((int128)tmp_q[5] * 867413771345704L) - ((int128)tmp_q[6] * 205818738405668L) + ((((int128)tmp_q[7] * 531005014050630L) - ((int128)tmp_q[8] * 1269139445004091L) - ((int128)tmp_q[9] * 2230351078777709L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 2230351078777709L) + ((int128)tmp_q[1] * 413756925522799L) - ((int128)tmp_q[2] * 1570316060621338L) - ((int128)tmp_q[3] * 2976118761151787L) + ((int128)tmp_q[4] * 1415810708167450L) + ((int128)tmp_q[5] * 501868486807353L) - ((int128)tmp_q[6] * 867413771345704L) - ((int128)tmp_q[7] * 205818738405668L) + ((((int128)tmp_q[8] * 531005014050630L) - ((int128)tmp_q[9] * 1269139445004091L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 1269139445004091L) - ((int128)tmp_q[1] * 2230351078777709L) + ((int128)tmp_q[2] * 413756925522799L) - ((int128)tmp_q[3] * 1570316060621338L) - ((int128)tmp_q[4] * 2976118761151787L) + ((int128)tmp_q[5] * 1415810708167450L) + ((int128)tmp_q[6] * 501868486807353L) - ((int128)tmp_q[7] * 867413771345704L) - ((int128)tmp_q[8] * 205818738405668L) + ((int128)tmp_q[9] * 1593015042151890L);
	tmp_zero[9] = ((int128)tmp_q[0] * 531005014050630L) - ((int128)tmp_q[1] * 1269139445004091L) - ((int128)tmp_q[2] * 2230351078777709L) + ((int128)tmp_q[3] * 413756925522799L) - ((int128)tmp_q[4] * 1570316060621338L) - ((int128)tmp_q[5] * 2976118761151787L) + ((int128)tmp_q[6] * 1415810708167450L) + ((int128)tmp_q[7] * 501868486807353L) - ((int128)tmp_q[8] * 867413771345704L) - ((int128)tmp_q[9] * 205818738405668L);

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
}

