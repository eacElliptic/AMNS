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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1342392087573307299UL) + ((((uint64_t)op[1] * 785005592044375762UL) + ((uint64_t)op[2] * 4105989362715071515UL) + ((uint64_t)op[3] * 4474165088651961327UL) + ((uint64_t)op[4] * 6513430505743457400UL) + ((uint64_t)op[5] * 6372140846422677635UL) + ((uint64_t)op[6] * 8575676787710100694UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 8575676787710100694UL) + ((uint64_t)op[1] * 1342392087573307299UL) + ((((uint64_t)op[2] * 785005592044375762UL) + ((uint64_t)op[3] * 4105989362715071515UL) + ((uint64_t)op[4] * 4474165088651961327UL) + ((uint64_t)op[5] * 6513430505743457400UL) + ((uint64_t)op[6] * 6372140846422677635UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 6372140846422677635UL) + ((uint64_t)op[1] * 8575676787710100694UL) + ((uint64_t)op[2] * 1342392087573307299UL) + ((((uint64_t)op[3] * 785005592044375762UL) + ((uint64_t)op[4] * 4105989362715071515UL) + ((uint64_t)op[5] * 4474165088651961327UL) + ((uint64_t)op[6] * 6513430505743457400UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 6513430505743457400UL) + ((uint64_t)op[1] * 6372140846422677635UL) + ((uint64_t)op[2] * 8575676787710100694UL) + ((uint64_t)op[3] * 1342392087573307299UL) + ((((uint64_t)op[4] * 785005592044375762UL) + ((uint64_t)op[5] * 4105989362715071515UL) + ((uint64_t)op[6] * 4474165088651961327UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 4474165088651961327UL) + ((uint64_t)op[1] * 6513430505743457400UL) + ((uint64_t)op[2] * 6372140846422677635UL) + ((uint64_t)op[3] * 8575676787710100694UL) + ((uint64_t)op[4] * 1342392087573307299UL) + ((((uint64_t)op[5] * 785005592044375762UL) + ((uint64_t)op[6] * 4105989362715071515UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 4105989362715071515UL) + ((uint64_t)op[1] * 4474165088651961327UL) + ((uint64_t)op[2] * 6513430505743457400UL) + ((uint64_t)op[3] * 6372140846422677635UL) + ((uint64_t)op[4] * 8575676787710100694UL) + ((uint64_t)op[5] * 1342392087573307299UL) + ((uint64_t)op[6] * 4710033552266254572UL);
	tmp_q[6] = ((uint64_t)op[0] * 785005592044375762UL) + ((uint64_t)op[1] * 4105989362715071515UL) + ((uint64_t)op[2] * 4474165088651961327UL) + ((uint64_t)op[3] * 6513430505743457400UL) + ((uint64_t)op[4] * 6372140846422677635UL) + ((uint64_t)op[5] * 8575676787710100694UL) + ((uint64_t)op[6] * 1342392087573307299UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 14801619713L) + ((((int128)tmp_q[1] * 41456764399L) - ((int128)tmp_q[2] * 5207877299L) - ((int128)tmp_q[3] * 14672030016L) - ((int128)tmp_q[4] * 32210764958L) + ((int128)tmp_q[5] * 25710372249L) - ((int128)tmp_q[6] * 78228429164L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 78228429164L) + ((int128)tmp_q[1] * 14801619713L) + ((((int128)tmp_q[2] * 41456764399L) - ((int128)tmp_q[3] * 5207877299L) - ((int128)tmp_q[4] * 14672030016L) - ((int128)tmp_q[5] * 32210764958L) + ((int128)tmp_q[6] * 25710372249L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 25710372249L) - ((int128)tmp_q[1] * 78228429164L) + ((int128)tmp_q[2] * 14801619713L) + ((((int128)tmp_q[3] * 41456764399L) - ((int128)tmp_q[4] * 5207877299L) - ((int128)tmp_q[5] * 14672030016L) - ((int128)tmp_q[6] * 32210764958L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 32210764958L) + ((int128)tmp_q[1] * 25710372249L) - ((int128)tmp_q[2] * 78228429164L) + ((int128)tmp_q[3] * 14801619713L) + ((((int128)tmp_q[4] * 41456764399L) - ((int128)tmp_q[5] * 5207877299L) - ((int128)tmp_q[6] * 14672030016L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 14672030016L) - ((int128)tmp_q[1] * 32210764958L) + ((int128)tmp_q[2] * 25710372249L) - ((int128)tmp_q[3] * 78228429164L) + ((int128)tmp_q[4] * 14801619713L) + ((((int128)tmp_q[5] * 41456764399L) - ((int128)tmp_q[6] * 5207877299L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 5207877299L) - ((int128)tmp_q[1] * 14672030016L) - ((int128)tmp_q[2] * 32210764958L) + ((int128)tmp_q[3] * 25710372249L) - ((int128)tmp_q[4] * 78228429164L) + ((int128)tmp_q[5] * 14801619713L) + ((int128)tmp_q[6] * 248740586394L);
	tmp_zero[6] = ((int128)tmp_q[0] * 41456764399L) - ((int128)tmp_q[1] * 5207877299L) - ((int128)tmp_q[2] * 14672030016L) - ((int128)tmp_q[3] * 32210764958L) + ((int128)tmp_q[4] * 25710372249L) - ((int128)tmp_q[5] * 78228429164L) + ((int128)tmp_q[6] * 14801619713L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

