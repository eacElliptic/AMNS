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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14207373603602980307UL) + ((((uint64_t)op[1] * 3984550111724488985UL) + ((uint64_t)op[2] * 12109492509314845416UL) + ((uint64_t)op[3] * 14897158686247977606UL) + ((uint64_t)op[4] * 15346283002996813592UL) + ((uint64_t)op[5] * 1159552766442702835UL) + ((uint64_t)op[6] * 15106703103783156746UL) + ((uint64_t)op[7] * 5909205327278115223UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 5909205327278115223UL) + ((uint64_t)op[1] * 14207373603602980307UL) + ((((uint64_t)op[2] * 3984550111724488985UL) + ((uint64_t)op[3] * 12109492509314845416UL) + ((uint64_t)op[4] * 14897158686247977606UL) + ((uint64_t)op[5] * 15346283002996813592UL) + ((uint64_t)op[6] * 1159552766442702835UL) + ((uint64_t)op[7] * 15106703103783156746UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 15106703103783156746UL) + ((uint64_t)op[1] * 5909205327278115223UL) + ((uint64_t)op[2] * 14207373603602980307UL) + ((((uint64_t)op[3] * 3984550111724488985UL) + ((uint64_t)op[4] * 12109492509314845416UL) + ((uint64_t)op[5] * 14897158686247977606UL) + ((uint64_t)op[6] * 15346283002996813592UL) + ((uint64_t)op[7] * 1159552766442702835UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 1159552766442702835UL) + ((uint64_t)op[1] * 15106703103783156746UL) + ((uint64_t)op[2] * 5909205327278115223UL) + ((uint64_t)op[3] * 14207373603602980307UL) + ((((uint64_t)op[4] * 3984550111724488985UL) + ((uint64_t)op[5] * 12109492509314845416UL) + ((uint64_t)op[6] * 14897158686247977606UL) + ((uint64_t)op[7] * 15346283002996813592UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 15346283002996813592UL) + ((uint64_t)op[1] * 1159552766442702835UL) + ((uint64_t)op[2] * 15106703103783156746UL) + ((uint64_t)op[3] * 5909205327278115223UL) + ((uint64_t)op[4] * 14207373603602980307UL) + ((((uint64_t)op[5] * 3984550111724488985UL) + ((uint64_t)op[6] * 12109492509314845416UL) + ((uint64_t)op[7] * 14897158686247977606UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 14897158686247977606UL) + ((uint64_t)op[1] * 15346283002996813592UL) + ((uint64_t)op[2] * 1159552766442702835UL) + ((uint64_t)op[3] * 15106703103783156746UL) + ((uint64_t)op[4] * 5909205327278115223UL) + ((uint64_t)op[5] * 14207373603602980307UL) + ((((uint64_t)op[6] * 3984550111724488985UL) + ((uint64_t)op[7] * 12109492509314845416UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 12109492509314845416UL) + ((uint64_t)op[1] * 14897158686247977606UL) + ((uint64_t)op[2] * 15346283002996813592UL) + ((uint64_t)op[3] * 1159552766442702835UL) + ((uint64_t)op[4] * 15106703103783156746UL) + ((uint64_t)op[5] * 5909205327278115223UL) + ((uint64_t)op[6] * 14207373603602980307UL) + ((uint64_t)op[7] * 5460556596637382294UL);
	tmp_q[7] = ((uint64_t)op[0] * 3984550111724488985UL) + ((uint64_t)op[1] * 12109492509314845416UL) + ((uint64_t)op[2] * 14897158686247977606UL) + ((uint64_t)op[3] * 15346283002996813592UL) + ((uint64_t)op[4] * 1159552766442702835UL) + ((uint64_t)op[5] * 15106703103783156746UL) + ((uint64_t)op[6] * 5909205327278115223UL) + ((uint64_t)op[7] * 14207373603602980307UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 112977682936971L) + ((((int128)tmp_q[1] * 45078279607728L) - ((int128)tmp_q[2] * 106514114339288L) - ((int128)tmp_q[3] * 108215926500914L) + ((int128)tmp_q[4] * 124024656867741L) - ((int128)tmp_q[5] * 81357654908066L) + ((int128)tmp_q[6] * 67553458943191L) - ((int128)tmp_q[7] * 84621785157017L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 84621785157017L) + ((int128)tmp_q[1] * 112977682936971L) + ((((int128)tmp_q[2] * 45078279607728L) - ((int128)tmp_q[3] * 106514114339288L) - ((int128)tmp_q[4] * 108215926500914L) + ((int128)tmp_q[5] * 124024656867741L) - ((int128)tmp_q[6] * 81357654908066L) + ((int128)tmp_q[7] * 67553458943191L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 67553458943191L) - ((int128)tmp_q[1] * 84621785157017L) + ((int128)tmp_q[2] * 112977682936971L) + ((((int128)tmp_q[3] * 45078279607728L) - ((int128)tmp_q[4] * 106514114339288L) - ((int128)tmp_q[5] * 108215926500914L) + ((int128)tmp_q[6] * 124024656867741L) - ((int128)tmp_q[7] * 81357654908066L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 81357654908066L) + ((int128)tmp_q[1] * 67553458943191L) - ((int128)tmp_q[2] * 84621785157017L) + ((int128)tmp_q[3] * 112977682936971L) + ((((int128)tmp_q[4] * 45078279607728L) - ((int128)tmp_q[5] * 106514114339288L) - ((int128)tmp_q[6] * 108215926500914L) + ((int128)tmp_q[7] * 124024656867741L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 124024656867741L) - ((int128)tmp_q[1] * 81357654908066L) + ((int128)tmp_q[2] * 67553458943191L) - ((int128)tmp_q[3] * 84621785157017L) + ((int128)tmp_q[4] * 112977682936971L) + ((((int128)tmp_q[5] * 45078279607728L) - ((int128)tmp_q[6] * 106514114339288L) - ((int128)tmp_q[7] * 108215926500914L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 108215926500914L) + ((int128)tmp_q[1] * 124024656867741L) - ((int128)tmp_q[2] * 81357654908066L) + ((int128)tmp_q[3] * 67553458943191L) - ((int128)tmp_q[4] * 84621785157017L) + ((int128)tmp_q[5] * 112977682936971L) + ((((int128)tmp_q[6] * 45078279607728L) - ((int128)tmp_q[7] * 106514114339288L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 106514114339288L) - ((int128)tmp_q[1] * 108215926500914L) + ((int128)tmp_q[2] * 124024656867741L) - ((int128)tmp_q[3] * 81357654908066L) + ((int128)tmp_q[4] * 67553458943191L) - ((int128)tmp_q[5] * 84621785157017L) + ((int128)tmp_q[6] * 112977682936971L) + ((int128)tmp_q[7] * 270469677646368L);
	tmp_zero[7] = ((int128)tmp_q[0] * 45078279607728L) - ((int128)tmp_q[1] * 106514114339288L) - ((int128)tmp_q[2] * 108215926500914L) + ((int128)tmp_q[3] * 124024656867741L) - ((int128)tmp_q[4] * 81357654908066L) + ((int128)tmp_q[5] * 67553458943191L) - ((int128)tmp_q[6] * 84621785157017L) + ((int128)tmp_q[7] * 112977682936971L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

