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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2611166586220033787UL) + ((((uint64_t)op[1] * 7967599796097696751UL) + ((uint64_t)op[2] * 10389698513064249011UL) + ((uint64_t)op[3] * 9259980620943013670UL) + ((uint64_t)op[4] * 13383300303802580502UL) + ((uint64_t)op[5] * 6463218713434573019UL) + ((uint64_t)op[6] * 14141490237880390755UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 14141490237880390755UL) + ((uint64_t)op[1] * 2611166586220033787UL) + ((((uint64_t)op[2] * 7967599796097696751UL) + ((uint64_t)op[3] * 10389698513064249011UL) + ((uint64_t)op[4] * 9259980620943013670UL) + ((uint64_t)op[5] * 13383300303802580502UL) + ((uint64_t)op[6] * 6463218713434573019UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 6463218713434573019UL) + ((uint64_t)op[1] * 14141490237880390755UL) + ((uint64_t)op[2] * 2611166586220033787UL) + ((((uint64_t)op[3] * 7967599796097696751UL) + ((uint64_t)op[4] * 10389698513064249011UL) + ((uint64_t)op[5] * 9259980620943013670UL) + ((uint64_t)op[6] * 13383300303802580502UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 13383300303802580502UL) + ((uint64_t)op[1] * 6463218713434573019UL) + ((uint64_t)op[2] * 14141490237880390755UL) + ((uint64_t)op[3] * 2611166586220033787UL) + ((((uint64_t)op[4] * 7967599796097696751UL) + ((uint64_t)op[5] * 10389698513064249011UL) + ((uint64_t)op[6] * 9259980620943013670UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 9259980620943013670UL) + ((uint64_t)op[1] * 13383300303802580502UL) + ((uint64_t)op[2] * 6463218713434573019UL) + ((uint64_t)op[3] * 14141490237880390755UL) + ((uint64_t)op[4] * 2611166586220033787UL) + ((((uint64_t)op[5] * 7967599796097696751UL) + ((uint64_t)op[6] * 10389698513064249011UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 10389698513064249011UL) + ((uint64_t)op[1] * 9259980620943013670UL) + ((uint64_t)op[2] * 13383300303802580502UL) + ((uint64_t)op[3] * 6463218713434573019UL) + ((uint64_t)op[4] * 14141490237880390755UL) + ((uint64_t)op[5] * 2611166586220033787UL) + ((uint64_t)op[6] * 13423655110681235388UL);
	tmp_q[6] = ((uint64_t)op[0] * 7967599796097696751UL) + ((uint64_t)op[1] * 10389698513064249011UL) + ((uint64_t)op[2] * 9259980620943013670UL) + ((uint64_t)op[3] * 13383300303802580502UL) + ((uint64_t)op[4] * 6463218713434573019UL) + ((uint64_t)op[5] * 14141490237880390755UL) + ((uint64_t)op[6] * 2611166586220033787UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 27219511805L) + ((((int128)tmp_q[1] * 8565813928L) - ((int128)tmp_q[2] * 38850065819L) + ((int128)tmp_q[3] * 39294227485L) - ((int128)tmp_q[4] * 13231950869L) - ((int128)tmp_q[5] * 12045646460L) + ((int128)tmp_q[6] * 25117591903L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 25117591903L) + ((int128)tmp_q[1] * 27219511805L) + ((((int128)tmp_q[2] * 8565813928L) - ((int128)tmp_q[3] * 38850065819L) + ((int128)tmp_q[4] * 39294227485L) - ((int128)tmp_q[5] * 13231950869L) - ((int128)tmp_q[6] * 12045646460L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 12045646460L) + ((int128)tmp_q[1] * 25117591903L) + ((int128)tmp_q[2] * 27219511805L) + ((((int128)tmp_q[3] * 8565813928L) - ((int128)tmp_q[4] * 38850065819L) + ((int128)tmp_q[5] * 39294227485L) - ((int128)tmp_q[6] * 13231950869L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 13231950869L) - ((int128)tmp_q[1] * 12045646460L) + ((int128)tmp_q[2] * 25117591903L) + ((int128)tmp_q[3] * 27219511805L) + ((((int128)tmp_q[4] * 8565813928L) - ((int128)tmp_q[5] * 38850065819L) + ((int128)tmp_q[6] * 39294227485L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 39294227485L) - ((int128)tmp_q[1] * 13231950869L) - ((int128)tmp_q[2] * 12045646460L) + ((int128)tmp_q[3] * 25117591903L) + ((int128)tmp_q[4] * 27219511805L) + ((((int128)tmp_q[5] * 8565813928L) - ((int128)tmp_q[6] * 38850065819L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 38850065819L) + ((int128)tmp_q[1] * 39294227485L) - ((int128)tmp_q[2] * 13231950869L) - ((int128)tmp_q[3] * 12045646460L) + ((int128)tmp_q[4] * 25117591903L) + ((int128)tmp_q[5] * 27219511805L) + ((int128)tmp_q[6] * 34263255712L);
	tmp_zero[6] = ((int128)tmp_q[0] * 8565813928L) - ((int128)tmp_q[1] * 38850065819L) + ((int128)tmp_q[2] * 39294227485L) - ((int128)tmp_q[3] * 13231950869L) - ((int128)tmp_q[4] * 12045646460L) + ((int128)tmp_q[5] * 25117591903L) + ((int128)tmp_q[6] * 27219511805L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

