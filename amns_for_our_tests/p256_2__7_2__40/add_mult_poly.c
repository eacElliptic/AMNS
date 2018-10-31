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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14966340066937021921UL) + ((((uint64_t)op[1] * 17432004668644515185UL) + ((uint64_t)op[2] * 3451323614492803150UL) + ((uint64_t)op[3] * 15138276147248897709UL) + ((uint64_t)op[4] * 4915807186348450191UL) + ((uint64_t)op[5] * 9623583487082550129UL) + ((uint64_t)op[6] * 11041100488758218792UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 11041100488758218792UL) + ((uint64_t)op[1] * 14966340066937021921UL) + ((((uint64_t)op[2] * 17432004668644515185UL) + ((uint64_t)op[3] * 3451323614492803150UL) + ((uint64_t)op[4] * 15138276147248897709UL) + ((uint64_t)op[5] * 4915807186348450191UL) + ((uint64_t)op[6] * 9623583487082550129UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 9623583487082550129UL) + ((uint64_t)op[1] * 11041100488758218792UL) + ((uint64_t)op[2] * 14966340066937021921UL) + ((((uint64_t)op[3] * 17432004668644515185UL) + ((uint64_t)op[4] * 3451323614492803150UL) + ((uint64_t)op[5] * 15138276147248897709UL) + ((uint64_t)op[6] * 4915807186348450191UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 4915807186348450191UL) + ((uint64_t)op[1] * 9623583487082550129UL) + ((uint64_t)op[2] * 11041100488758218792UL) + ((uint64_t)op[3] * 14966340066937021921UL) + ((((uint64_t)op[4] * 17432004668644515185UL) + ((uint64_t)op[5] * 3451323614492803150UL) + ((uint64_t)op[6] * 15138276147248897709UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 15138276147248897709UL) + ((uint64_t)op[1] * 4915807186348450191UL) + ((uint64_t)op[2] * 9623583487082550129UL) + ((uint64_t)op[3] * 11041100488758218792UL) + ((uint64_t)op[4] * 14966340066937021921UL) + ((((uint64_t)op[5] * 17432004668644515185UL) + ((uint64_t)op[6] * 3451323614492803150UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 3451323614492803150UL) + ((uint64_t)op[1] * 15138276147248897709UL) + ((uint64_t)op[2] * 4915807186348450191UL) + ((uint64_t)op[3] * 9623583487082550129UL) + ((uint64_t)op[4] * 11041100488758218792UL) + ((uint64_t)op[5] * 14966340066937021921UL) + ((uint64_t)op[6] * 16417265263579478754UL);
	tmp_q[6] = ((uint64_t)op[0] * 17432004668644515185UL) + ((uint64_t)op[1] * 3451323614492803150UL) + ((uint64_t)op[2] * 15138276147248897709UL) + ((uint64_t)op[3] * 4915807186348450191UL) + ((uint64_t)op[4] * 9623583487082550129UL) + ((uint64_t)op[5] * 11041100488758218792UL) + ((uint64_t)op[6] * 14966340066937021921UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 34783750549L) + ((((int128)tmp_q[1] * 4828442803L) - ((int128)tmp_q[2] * 44735075032L) + ((int128)tmp_q[3] * 51642172808L) + ((int128)tmp_q[4] * 2096039303L) + ((int128)tmp_q[5] * 35611535683L) + ((int128)tmp_q[6] * 17099166664L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 17099166664L) + ((int128)tmp_q[1] * 34783750549L) + ((((int128)tmp_q[2] * 4828442803L) - ((int128)tmp_q[3] * 44735075032L) + ((int128)tmp_q[4] * 51642172808L) + ((int128)tmp_q[5] * 2096039303L) + ((int128)tmp_q[6] * 35611535683L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 35611535683L) + ((int128)tmp_q[1] * 17099166664L) + ((int128)tmp_q[2] * 34783750549L) + ((((int128)tmp_q[3] * 4828442803L) - ((int128)tmp_q[4] * 44735075032L) + ((int128)tmp_q[5] * 51642172808L) + ((int128)tmp_q[6] * 2096039303L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 2096039303L) + ((int128)tmp_q[1] * 35611535683L) + ((int128)tmp_q[2] * 17099166664L) + ((int128)tmp_q[3] * 34783750549L) + ((((int128)tmp_q[4] * 4828442803L) - ((int128)tmp_q[5] * 44735075032L) + ((int128)tmp_q[6] * 51642172808L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 51642172808L) + ((int128)tmp_q[1] * 2096039303L) + ((int128)tmp_q[2] * 35611535683L) + ((int128)tmp_q[3] * 17099166664L) + ((int128)tmp_q[4] * 34783750549L) + ((((int128)tmp_q[5] * 4828442803L) - ((int128)tmp_q[6] * 44735075032L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 44735075032L) + ((int128)tmp_q[1] * 51642172808L) + ((int128)tmp_q[2] * 2096039303L) + ((int128)tmp_q[3] * 35611535683L) + ((int128)tmp_q[4] * 17099166664L) + ((int128)tmp_q[5] * 34783750549L) + ((int128)tmp_q[6] * 9656885606L);
	tmp_zero[6] = ((int128)tmp_q[0] * 4828442803L) - ((int128)tmp_q[1] * 44735075032L) + ((int128)tmp_q[2] * 51642172808L) + ((int128)tmp_q[3] * 2096039303L) + ((int128)tmp_q[4] * 35611535683L) + ((int128)tmp_q[5] * 17099166664L) + ((int128)tmp_q[6] * 34783750549L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

