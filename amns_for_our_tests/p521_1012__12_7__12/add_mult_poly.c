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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 7);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 14);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 14);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 7);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15004686155748297641UL) + ((((uint64_t)op[1] * 9551910069294305011UL) + ((uint64_t)op[2] * 5873042267339820718UL) + ((uint64_t)op[3] * 15679845951070022575UL) + ((uint64_t)op[4] * 15904532596375521215UL) + ((uint64_t)op[5] * 18327041132185567847UL) + ((uint64_t)op[6] * 7535261698953887950UL) + ((uint64_t)op[7] * 13966099157140589716UL) + ((uint64_t)op[8] * 1415028219404108801UL) + ((uint64_t)op[9] * 10028544907167696858UL) + ((uint64_t)op[10] * 17117751258239820495UL) + ((uint64_t)op[11] * 7869062352608197294UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 7869062352608197294UL) + ((uint64_t)op[1] * 15004686155748297641UL) + ((((uint64_t)op[2] * 9551910069294305011UL) + ((uint64_t)op[3] * 5873042267339820718UL) + ((uint64_t)op[4] * 15679845951070022575UL) + ((uint64_t)op[5] * 15904532596375521215UL) + ((uint64_t)op[6] * 18327041132185567847UL) + ((uint64_t)op[7] * 7535261698953887950UL) + ((uint64_t)op[8] * 13966099157140589716UL) + ((uint64_t)op[9] * 1415028219404108801UL) + ((uint64_t)op[10] * 10028544907167696858UL) + ((uint64_t)op[11] * 17117751258239820495UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 17117751258239820495UL) + ((uint64_t)op[1] * 7869062352608197294UL) + ((uint64_t)op[2] * 15004686155748297641UL) + ((((uint64_t)op[3] * 9551910069294305011UL) + ((uint64_t)op[4] * 5873042267339820718UL) + ((uint64_t)op[5] * 15679845951070022575UL) + ((uint64_t)op[6] * 15904532596375521215UL) + ((uint64_t)op[7] * 18327041132185567847UL) + ((uint64_t)op[8] * 7535261698953887950UL) + ((uint64_t)op[9] * 13966099157140589716UL) + ((uint64_t)op[10] * 1415028219404108801UL) + ((uint64_t)op[11] * 10028544907167696858UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 10028544907167696858UL) + ((uint64_t)op[1] * 17117751258239820495UL) + ((uint64_t)op[2] * 7869062352608197294UL) + ((uint64_t)op[3] * 15004686155748297641UL) + ((((uint64_t)op[4] * 9551910069294305011UL) + ((uint64_t)op[5] * 5873042267339820718UL) + ((uint64_t)op[6] * 15679845951070022575UL) + ((uint64_t)op[7] * 15904532596375521215UL) + ((uint64_t)op[8] * 18327041132185567847UL) + ((uint64_t)op[9] * 7535261698953887950UL) + ((uint64_t)op[10] * 13966099157140589716UL) + ((uint64_t)op[11] * 1415028219404108801UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 1415028219404108801UL) + ((uint64_t)op[1] * 10028544907167696858UL) + ((uint64_t)op[2] * 17117751258239820495UL) + ((uint64_t)op[3] * 7869062352608197294UL) + ((uint64_t)op[4] * 15004686155748297641UL) + ((((uint64_t)op[5] * 9551910069294305011UL) + ((uint64_t)op[6] * 5873042267339820718UL) + ((uint64_t)op[7] * 15679845951070022575UL) + ((uint64_t)op[8] * 15904532596375521215UL) + ((uint64_t)op[9] * 18327041132185567847UL) + ((uint64_t)op[10] * 7535261698953887950UL) + ((uint64_t)op[11] * 13966099157140589716UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 13966099157140589716UL) + ((uint64_t)op[1] * 1415028219404108801UL) + ((uint64_t)op[2] * 10028544907167696858UL) + ((uint64_t)op[3] * 17117751258239820495UL) + ((uint64_t)op[4] * 7869062352608197294UL) + ((uint64_t)op[5] * 15004686155748297641UL) + ((((uint64_t)op[6] * 9551910069294305011UL) + ((uint64_t)op[7] * 5873042267339820718UL) + ((uint64_t)op[8] * 15679845951070022575UL) + ((uint64_t)op[9] * 15904532596375521215UL) + ((uint64_t)op[10] * 18327041132185567847UL) + ((uint64_t)op[11] * 7535261698953887950UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 7535261698953887950UL) + ((uint64_t)op[1] * 13966099157140589716UL) + ((uint64_t)op[2] * 1415028219404108801UL) + ((uint64_t)op[3] * 10028544907167696858UL) + ((uint64_t)op[4] * 17117751258239820495UL) + ((uint64_t)op[5] * 7869062352608197294UL) + ((uint64_t)op[6] * 15004686155748297641UL) + ((((uint64_t)op[7] * 9551910069294305011UL) + ((uint64_t)op[8] * 5873042267339820718UL) + ((uint64_t)op[9] * 15679845951070022575UL) + ((uint64_t)op[10] * 15904532596375521215UL) + ((uint64_t)op[11] * 18327041132185567847UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 18327041132185567847UL) + ((uint64_t)op[1] * 7535261698953887950UL) + ((uint64_t)op[2] * 13966099157140589716UL) + ((uint64_t)op[3] * 1415028219404108801UL) + ((uint64_t)op[4] * 10028544907167696858UL) + ((uint64_t)op[5] * 17117751258239820495UL) + ((uint64_t)op[6] * 7869062352608197294UL) + ((uint64_t)op[7] * 15004686155748297641UL) + ((((uint64_t)op[8] * 9551910069294305011UL) + ((uint64_t)op[9] * 5873042267339820718UL) + ((uint64_t)op[10] * 15679845951070022575UL) + ((uint64_t)op[11] * 15904532596375521215UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 15904532596375521215UL) + ((uint64_t)op[1] * 18327041132185567847UL) + ((uint64_t)op[2] * 7535261698953887950UL) + ((uint64_t)op[3] * 13966099157140589716UL) + ((uint64_t)op[4] * 1415028219404108801UL) + ((uint64_t)op[5] * 10028544907167696858UL) + ((uint64_t)op[6] * 17117751258239820495UL) + ((uint64_t)op[7] * 7869062352608197294UL) + ((uint64_t)op[8] * 15004686155748297641UL) + ((((uint64_t)op[9] * 9551910069294305011UL) + ((uint64_t)op[10] * 5873042267339820718UL) + ((uint64_t)op[11] * 15679845951070022575UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 15679845951070022575UL) + ((uint64_t)op[1] * 15904532596375521215UL) + ((uint64_t)op[2] * 18327041132185567847UL) + ((uint64_t)op[3] * 7535261698953887950UL) + ((uint64_t)op[4] * 13966099157140589716UL) + ((uint64_t)op[5] * 1415028219404108801UL) + ((uint64_t)op[6] * 10028544907167696858UL) + ((uint64_t)op[7] * 17117751258239820495UL) + ((uint64_t)op[8] * 7869062352608197294UL) + ((uint64_t)op[9] * 15004686155748297641UL) + ((((uint64_t)op[10] * 9551910069294305011UL) + ((uint64_t)op[11] * 5873042267339820718UL)) * 7);
	tmp_q[10] = ((uint64_t)op[0] * 5873042267339820718UL) + ((uint64_t)op[1] * 15679845951070022575UL) + ((uint64_t)op[2] * 15904532596375521215UL) + ((uint64_t)op[3] * 18327041132185567847UL) + ((uint64_t)op[4] * 7535261698953887950UL) + ((uint64_t)op[5] * 13966099157140589716UL) + ((uint64_t)op[6] * 1415028219404108801UL) + ((uint64_t)op[7] * 10028544907167696858UL) + ((uint64_t)op[8] * 17117751258239820495UL) + ((uint64_t)op[9] * 7869062352608197294UL) + ((uint64_t)op[10] * 15004686155748297641UL) + ((uint64_t)op[11] * 11523138263931480229UL);
	tmp_q[11] = ((uint64_t)op[0] * 9551910069294305011UL) + ((uint64_t)op[1] * 5873042267339820718UL) + ((uint64_t)op[2] * 15679845951070022575UL) + ((uint64_t)op[3] * 15904532596375521215UL) + ((uint64_t)op[4] * 18327041132185567847UL) + ((uint64_t)op[5] * 7535261698953887950UL) + ((uint64_t)op[6] * 13966099157140589716UL) + ((uint64_t)op[7] * 1415028219404108801UL) + ((uint64_t)op[8] * 10028544907167696858UL) + ((uint64_t)op[9] * 17117751258239820495UL) + ((uint64_t)op[10] * 7869062352608197294UL) + ((uint64_t)op[11] * 15004686155748297641UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4727038440499L) + ((((int128)tmp_q[1] * 2908338110367L) + ((int128)tmp_q[2] * 2594579387201L) - ((int128)tmp_q[3] * 779605575587L) - ((int128)tmp_q[4] * 4295741551757L) + ((int128)tmp_q[5] * 952137505856L) + ((int128)tmp_q[6] * 2386794526610L) - ((int128)tmp_q[7] * 908524468376L) + ((int128)tmp_q[8] * 3059492547223L) - ((int128)tmp_q[9] * 4701607637524L) + ((int128)tmp_q[10] * 46643499920L) - ((int128)tmp_q[11] * 1843550262155L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 1843550262155L) - ((int128)tmp_q[1] * 4727038440499L) + ((((int128)tmp_q[2] * 2908338110367L) + ((int128)tmp_q[3] * 2594579387201L) - ((int128)tmp_q[4] * 779605575587L) - ((int128)tmp_q[5] * 4295741551757L) + ((int128)tmp_q[6] * 952137505856L) + ((int128)tmp_q[7] * 2386794526610L) - ((int128)tmp_q[8] * 908524468376L) + ((int128)tmp_q[9] * 3059492547223L) - ((int128)tmp_q[10] * 4701607637524L) + ((int128)tmp_q[11] * 46643499920L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 46643499920L) - ((int128)tmp_q[1] * 1843550262155L) - ((int128)tmp_q[2] * 4727038440499L) + ((((int128)tmp_q[3] * 2908338110367L) + ((int128)tmp_q[4] * 2594579387201L) - ((int128)tmp_q[5] * 779605575587L) - ((int128)tmp_q[6] * 4295741551757L) + ((int128)tmp_q[7] * 952137505856L) + ((int128)tmp_q[8] * 2386794526610L) - ((int128)tmp_q[9] * 908524468376L) + ((int128)tmp_q[10] * 3059492547223L) - ((int128)tmp_q[11] * 4701607637524L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 4701607637524L) + ((int128)tmp_q[1] * 46643499920L) - ((int128)tmp_q[2] * 1843550262155L) - ((int128)tmp_q[3] * 4727038440499L) + ((((int128)tmp_q[4] * 2908338110367L) + ((int128)tmp_q[5] * 2594579387201L) - ((int128)tmp_q[6] * 779605575587L) - ((int128)tmp_q[7] * 4295741551757L) + ((int128)tmp_q[8] * 952137505856L) + ((int128)tmp_q[9] * 2386794526610L) - ((int128)tmp_q[10] * 908524468376L) + ((int128)tmp_q[11] * 3059492547223L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 3059492547223L) - ((int128)tmp_q[1] * 4701607637524L) + ((int128)tmp_q[2] * 46643499920L) - ((int128)tmp_q[3] * 1843550262155L) - ((int128)tmp_q[4] * 4727038440499L) + ((((int128)tmp_q[5] * 2908338110367L) + ((int128)tmp_q[6] * 2594579387201L) - ((int128)tmp_q[7] * 779605575587L) - ((int128)tmp_q[8] * 4295741551757L) + ((int128)tmp_q[9] * 952137505856L) + ((int128)tmp_q[10] * 2386794526610L) - ((int128)tmp_q[11] * 908524468376L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 908524468376L) + ((int128)tmp_q[1] * 3059492547223L) - ((int128)tmp_q[2] * 4701607637524L) + ((int128)tmp_q[3] * 46643499920L) - ((int128)tmp_q[4] * 1843550262155L) - ((int128)tmp_q[5] * 4727038440499L) + ((((int128)tmp_q[6] * 2908338110367L) + ((int128)tmp_q[7] * 2594579387201L) - ((int128)tmp_q[8] * 779605575587L) - ((int128)tmp_q[9] * 4295741551757L) + ((int128)tmp_q[10] * 952137505856L) + ((int128)tmp_q[11] * 2386794526610L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 2386794526610L) - ((int128)tmp_q[1] * 908524468376L) + ((int128)tmp_q[2] * 3059492547223L) - ((int128)tmp_q[3] * 4701607637524L) + ((int128)tmp_q[4] * 46643499920L) - ((int128)tmp_q[5] * 1843550262155L) - ((int128)tmp_q[6] * 4727038440499L) + ((((int128)tmp_q[7] * 2908338110367L) + ((int128)tmp_q[8] * 2594579387201L) - ((int128)tmp_q[9] * 779605575587L) - ((int128)tmp_q[10] * 4295741551757L) + ((int128)tmp_q[11] * 952137505856L)) * 7);
	tmp_zero[7] = ((int128)tmp_q[0] * 952137505856L) + ((int128)tmp_q[1] * 2386794526610L) - ((int128)tmp_q[2] * 908524468376L) + ((int128)tmp_q[3] * 3059492547223L) - ((int128)tmp_q[4] * 4701607637524L) + ((int128)tmp_q[5] * 46643499920L) - ((int128)tmp_q[6] * 1843550262155L) - ((int128)tmp_q[7] * 4727038440499L) + ((((int128)tmp_q[8] * 2908338110367L) + ((int128)tmp_q[9] * 2594579387201L) - ((int128)tmp_q[10] * 779605575587L) - ((int128)tmp_q[11] * 4295741551757L)) * 7);
	tmp_zero[8] = -((int128)tmp_q[0] * 4295741551757L) + ((int128)tmp_q[1] * 952137505856L) + ((int128)tmp_q[2] * 2386794526610L) - ((int128)tmp_q[3] * 908524468376L) + ((int128)tmp_q[4] * 3059492547223L) - ((int128)tmp_q[5] * 4701607637524L) + ((int128)tmp_q[6] * 46643499920L) - ((int128)tmp_q[7] * 1843550262155L) - ((int128)tmp_q[8] * 4727038440499L) + ((((int128)tmp_q[9] * 2908338110367L) + ((int128)tmp_q[10] * 2594579387201L) - ((int128)tmp_q[11] * 779605575587L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 779605575587L) - ((int128)tmp_q[1] * 4295741551757L) + ((int128)tmp_q[2] * 952137505856L) + ((int128)tmp_q[3] * 2386794526610L) - ((int128)tmp_q[4] * 908524468376L) + ((int128)tmp_q[5] * 3059492547223L) - ((int128)tmp_q[6] * 4701607637524L) + ((int128)tmp_q[7] * 46643499920L) - ((int128)tmp_q[8] * 1843550262155L) - ((int128)tmp_q[9] * 4727038440499L) + ((((int128)tmp_q[10] * 2908338110367L) + ((int128)tmp_q[11] * 2594579387201L)) * 7);
	tmp_zero[10] = ((int128)tmp_q[0] * 2594579387201L) - ((int128)tmp_q[1] * 779605575587L) - ((int128)tmp_q[2] * 4295741551757L) + ((int128)tmp_q[3] * 952137505856L) + ((int128)tmp_q[4] * 2386794526610L) - ((int128)tmp_q[5] * 908524468376L) + ((int128)tmp_q[6] * 3059492547223L) - ((int128)tmp_q[7] * 4701607637524L) + ((int128)tmp_q[8] * 46643499920L) - ((int128)tmp_q[9] * 1843550262155L) - ((int128)tmp_q[10] * 4727038440499L) + ((int128)tmp_q[11] * 20358366772569L);
	tmp_zero[11] = ((int128)tmp_q[0] * 2908338110367L) + ((int128)tmp_q[1] * 2594579387201L) - ((int128)tmp_q[2] * 779605575587L) - ((int128)tmp_q[3] * 4295741551757L) + ((int128)tmp_q[4] * 952137505856L) + ((int128)tmp_q[5] * 2386794526610L) - ((int128)tmp_q[6] * 908524468376L) + ((int128)tmp_q[7] * 3059492547223L) - ((int128)tmp_q[8] * 4701607637524L) + ((int128)tmp_q[9] * 46643499920L) - ((int128)tmp_q[10] * 1843550262155L) - ((int128)tmp_q[11] * 4727038440499L);

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

