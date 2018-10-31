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
	tmp_q[0] = ((uint64_t)op[0] * 10593182930817444001UL) + ((((uint64_t)op[1] * 6991430668053757292UL) + ((uint64_t)op[2] * 8585680752715444991UL) + ((uint64_t)op[3] * 10436420730228955740UL) + ((uint64_t)op[4] * 8845751890425418217UL) + ((uint64_t)op[5] * 13379103670735250557UL) + ((uint64_t)op[6] * 8687714705853980167UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 8687714705853980167UL) + ((uint64_t)op[1] * 10593182930817444001UL) + ((((uint64_t)op[2] * 6991430668053757292UL) + ((uint64_t)op[3] * 8585680752715444991UL) + ((uint64_t)op[4] * 10436420730228955740UL) + ((uint64_t)op[5] * 8845751890425418217UL) + ((uint64_t)op[6] * 13379103670735250557UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 13379103670735250557UL) + ((uint64_t)op[1] * 8687714705853980167UL) + ((uint64_t)op[2] * 10593182930817444001UL) + ((((uint64_t)op[3] * 6991430668053757292UL) + ((uint64_t)op[4] * 8585680752715444991UL) + ((uint64_t)op[5] * 10436420730228955740UL) + ((uint64_t)op[6] * 8845751890425418217UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 8845751890425418217UL) + ((uint64_t)op[1] * 13379103670735250557UL) + ((uint64_t)op[2] * 8687714705853980167UL) + ((uint64_t)op[3] * 10593182930817444001UL) + ((((uint64_t)op[4] * 6991430668053757292UL) + ((uint64_t)op[5] * 8585680752715444991UL) + ((uint64_t)op[6] * 10436420730228955740UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 10436420730228955740UL) + ((uint64_t)op[1] * 8845751890425418217UL) + ((uint64_t)op[2] * 13379103670735250557UL) + ((uint64_t)op[3] * 8687714705853980167UL) + ((uint64_t)op[4] * 10593182930817444001UL) + ((((uint64_t)op[5] * 6991430668053757292UL) + ((uint64_t)op[6] * 8585680752715444991UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 8585680752715444991UL) + ((uint64_t)op[1] * 10436420730228955740UL) + ((uint64_t)op[2] * 8845751890425418217UL) + ((uint64_t)op[3] * 13379103670735250557UL) + ((uint64_t)op[4] * 8687714705853980167UL) + ((uint64_t)op[5] * 10593182930817444001UL) + ((uint64_t)op[6] * 13982861336107514584UL);
	tmp_q[6] = ((uint64_t)op[0] * 6991430668053757292UL) + ((uint64_t)op[1] * 8585680752715444991UL) + ((uint64_t)op[2] * 10436420730228955740UL) + ((uint64_t)op[3] * 8845751890425418217UL) + ((uint64_t)op[4] * 13379103670735250557UL) + ((uint64_t)op[5] * 8687714705853980167UL) + ((uint64_t)op[6] * 10593182930817444001UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1760624279642671L) + ((-((int128)tmp_q[1] * 3564052639746776L) + ((int128)tmp_q[2] * 14759740875657580L) - ((int128)tmp_q[3] * 12523005840602003L) - ((int128)tmp_q[4] * 13496187834758096L) + ((int128)tmp_q[5] * 9776379938745998L) - ((int128)tmp_q[6] * 2298589562429383L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 2298589562429383L) - ((int128)tmp_q[1] * 1760624279642671L) + ((-((int128)tmp_q[2] * 3564052639746776L) + ((int128)tmp_q[3] * 14759740875657580L) - ((int128)tmp_q[4] * 12523005840602003L) - ((int128)tmp_q[5] * 13496187834758096L) + ((int128)tmp_q[6] * 9776379938745998L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 9776379938745998L) - ((int128)tmp_q[1] * 2298589562429383L) - ((int128)tmp_q[2] * 1760624279642671L) + ((-((int128)tmp_q[3] * 3564052639746776L) + ((int128)tmp_q[4] * 14759740875657580L) - ((int128)tmp_q[5] * 12523005840602003L) - ((int128)tmp_q[6] * 13496187834758096L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 13496187834758096L) + ((int128)tmp_q[1] * 9776379938745998L) - ((int128)tmp_q[2] * 2298589562429383L) - ((int128)tmp_q[3] * 1760624279642671L) + ((-((int128)tmp_q[4] * 3564052639746776L) + ((int128)tmp_q[5] * 14759740875657580L) - ((int128)tmp_q[6] * 12523005840602003L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 12523005840602003L) - ((int128)tmp_q[1] * 13496187834758096L) + ((int128)tmp_q[2] * 9776379938745998L) - ((int128)tmp_q[3] * 2298589562429383L) - ((int128)tmp_q[4] * 1760624279642671L) + ((-((int128)tmp_q[5] * 3564052639746776L) + ((int128)tmp_q[6] * 14759740875657580L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 14759740875657580L) - ((int128)tmp_q[1] * 12523005840602003L) - ((int128)tmp_q[2] * 13496187834758096L) + ((int128)tmp_q[3] * 9776379938745998L) - ((int128)tmp_q[4] * 2298589562429383L) - ((int128)tmp_q[5] * 1760624279642671L) - ((int128)tmp_q[6] * 7128105279493552L);
	tmp_zero[6] = -((int128)tmp_q[0] * 3564052639746776L) + ((int128)tmp_q[1] * 14759740875657580L) - ((int128)tmp_q[2] * 12523005840602003L) - ((int128)tmp_q[3] * 13496187834758096L) + ((int128)tmp_q[4] * 9776379938745998L) - ((int128)tmp_q[5] * 2298589562429383L) - ((int128)tmp_q[6] * 1760624279642671L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

