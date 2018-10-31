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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17802198286452356153UL) + ((((uint64_t)op[1] * 17272763891109815442UL) + ((uint64_t)op[2] * 18414468712676813519UL) + ((uint64_t)op[3] * 744986908316070046UL) + ((uint64_t)op[4] * 15911083834084484964UL) + ((uint64_t)op[5] * 7053412200005882876UL) + ((uint64_t)op[6] * 7823804806614365429UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 7823804806614365429UL) + ((uint64_t)op[1] * 17802198286452356153UL) + ((((uint64_t)op[2] * 17272763891109815442UL) + ((uint64_t)op[3] * 18414468712676813519UL) + ((uint64_t)op[4] * 744986908316070046UL) + ((uint64_t)op[5] * 15911083834084484964UL) + ((uint64_t)op[6] * 7053412200005882876UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 7053412200005882876UL) + ((uint64_t)op[1] * 7823804806614365429UL) + ((uint64_t)op[2] * 17802198286452356153UL) + ((((uint64_t)op[3] * 17272763891109815442UL) + ((uint64_t)op[4] * 18414468712676813519UL) + ((uint64_t)op[5] * 744986908316070046UL) + ((uint64_t)op[6] * 15911083834084484964UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 15911083834084484964UL) + ((uint64_t)op[1] * 7053412200005882876UL) + ((uint64_t)op[2] * 7823804806614365429UL) + ((uint64_t)op[3] * 17802198286452356153UL) + ((((uint64_t)op[4] * 17272763891109815442UL) + ((uint64_t)op[5] * 18414468712676813519UL) + ((uint64_t)op[6] * 744986908316070046UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 744986908316070046UL) + ((uint64_t)op[1] * 15911083834084484964UL) + ((uint64_t)op[2] * 7053412200005882876UL) + ((uint64_t)op[3] * 7823804806614365429UL) + ((uint64_t)op[4] * 17802198286452356153UL) + ((((uint64_t)op[5] * 17272763891109815442UL) + ((uint64_t)op[6] * 18414468712676813519UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 18414468712676813519UL) + ((uint64_t)op[1] * 744986908316070046UL) + ((uint64_t)op[2] * 15911083834084484964UL) + ((uint64_t)op[3] * 7053412200005882876UL) + ((uint64_t)op[4] * 7823804806614365429UL) + ((uint64_t)op[5] * 17802198286452356153UL) + ((uint64_t)op[6] * 9391841460797889392UL);
	tmp_q[6] = ((uint64_t)op[0] * 17272763891109815442UL) + ((uint64_t)op[1] * 18414468712676813519UL) + ((uint64_t)op[2] * 744986908316070046UL) + ((uint64_t)op[3] * 15911083834084484964UL) + ((uint64_t)op[4] * 7053412200005882876UL) + ((uint64_t)op[5] * 7823804806614365429UL) + ((uint64_t)op[6] * 17802198286452356153UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 41056485721L) - ((-((int128)tmp_q[1] * 36225101887L) - ((int128)tmp_q[2] * 37505666836L) + ((int128)tmp_q[3] * 26322985657L) - ((int128)tmp_q[4] * 12688030663L) - ((int128)tmp_q[5] * 6972143381L) - ((int128)tmp_q[6] * 50513110259L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 50513110259L) - ((int128)tmp_q[1] * 41056485721L) - ((-((int128)tmp_q[2] * 36225101887L) - ((int128)tmp_q[3] * 37505666836L) + ((int128)tmp_q[4] * 26322985657L) - ((int128)tmp_q[5] * 12688030663L) - ((int128)tmp_q[6] * 6972143381L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 6972143381L) - ((int128)tmp_q[1] * 50513110259L) - ((int128)tmp_q[2] * 41056485721L) - ((-((int128)tmp_q[3] * 36225101887L) - ((int128)tmp_q[4] * 37505666836L) + ((int128)tmp_q[5] * 26322985657L) - ((int128)tmp_q[6] * 12688030663L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 12688030663L) - ((int128)tmp_q[1] * 6972143381L) - ((int128)tmp_q[2] * 50513110259L) - ((int128)tmp_q[3] * 41056485721L) - ((-((int128)tmp_q[4] * 36225101887L) - ((int128)tmp_q[5] * 37505666836L) + ((int128)tmp_q[6] * 26322985657L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 26322985657L) - ((int128)tmp_q[1] * 12688030663L) - ((int128)tmp_q[2] * 6972143381L) - ((int128)tmp_q[3] * 50513110259L) - ((int128)tmp_q[4] * 41056485721L) - ((-((int128)tmp_q[5] * 36225101887L) - ((int128)tmp_q[6] * 37505666836L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 37505666836L) + ((int128)tmp_q[1] * 26322985657L) - ((int128)tmp_q[2] * 12688030663L) - ((int128)tmp_q[3] * 6972143381L) - ((int128)tmp_q[4] * 50513110259L) - ((int128)tmp_q[5] * 41056485721L) + ((int128)tmp_q[6] * 289800815096L);
	tmp_zero[6] = -((int128)tmp_q[0] * 36225101887L) - ((int128)tmp_q[1] * 37505666836L) + ((int128)tmp_q[2] * 26322985657L) - ((int128)tmp_q[3] * 12688030663L) - ((int128)tmp_q[4] * 6972143381L) - ((int128)tmp_q[5] * 50513110259L) - ((int128)tmp_q[6] * 41056485721L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

