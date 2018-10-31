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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2065473124114430057UL) + ((((uint64_t)op[1] * 10488027524736910128UL) + ((uint64_t)op[2] * 12457959752871040234UL) + ((uint64_t)op[3] * 11408867784380029229UL) + ((uint64_t)op[4] * 5013678683758370127UL) + ((uint64_t)op[5] * 16125149467419237595UL) + ((uint64_t)op[6] * 222469143476351708UL) + ((uint64_t)op[7] * 12667832759597063407UL) + ((uint64_t)op[8] * 10718469084878803840UL) + ((uint64_t)op[9] * 4929736269080567192UL) + ((uint64_t)op[10] * 2457931586228274314UL) + ((uint64_t)op[11] * 6110991910097450496UL)) * 14);
	tmp_q[1] = ((uint64_t)op[0] * 6110991910097450496UL) + ((uint64_t)op[1] * 2065473124114430057UL) + ((((uint64_t)op[2] * 10488027524736910128UL) + ((uint64_t)op[3] * 12457959752871040234UL) + ((uint64_t)op[4] * 11408867784380029229UL) + ((uint64_t)op[5] * 5013678683758370127UL) + ((uint64_t)op[6] * 16125149467419237595UL) + ((uint64_t)op[7] * 222469143476351708UL) + ((uint64_t)op[8] * 12667832759597063407UL) + ((uint64_t)op[9] * 10718469084878803840UL) + ((uint64_t)op[10] * 4929736269080567192UL) + ((uint64_t)op[11] * 2457931586228274314UL)) * 14);
	tmp_q[2] = ((uint64_t)op[0] * 2457931586228274314UL) + ((uint64_t)op[1] * 6110991910097450496UL) + ((uint64_t)op[2] * 2065473124114430057UL) + ((((uint64_t)op[3] * 10488027524736910128UL) + ((uint64_t)op[4] * 12457959752871040234UL) + ((uint64_t)op[5] * 11408867784380029229UL) + ((uint64_t)op[6] * 5013678683758370127UL) + ((uint64_t)op[7] * 16125149467419237595UL) + ((uint64_t)op[8] * 222469143476351708UL) + ((uint64_t)op[9] * 12667832759597063407UL) + ((uint64_t)op[10] * 10718469084878803840UL) + ((uint64_t)op[11] * 4929736269080567192UL)) * 14);
	tmp_q[3] = ((uint64_t)op[0] * 4929736269080567192UL) + ((uint64_t)op[1] * 2457931586228274314UL) + ((uint64_t)op[2] * 6110991910097450496UL) + ((uint64_t)op[3] * 2065473124114430057UL) + ((((uint64_t)op[4] * 10488027524736910128UL) + ((uint64_t)op[5] * 12457959752871040234UL) + ((uint64_t)op[6] * 11408867784380029229UL) + ((uint64_t)op[7] * 5013678683758370127UL) + ((uint64_t)op[8] * 16125149467419237595UL) + ((uint64_t)op[9] * 222469143476351708UL) + ((uint64_t)op[10] * 12667832759597063407UL) + ((uint64_t)op[11] * 10718469084878803840UL)) * 14);
	tmp_q[4] = ((uint64_t)op[0] * 10718469084878803840UL) + ((uint64_t)op[1] * 4929736269080567192UL) + ((uint64_t)op[2] * 2457931586228274314UL) + ((uint64_t)op[3] * 6110991910097450496UL) + ((uint64_t)op[4] * 2065473124114430057UL) + ((((uint64_t)op[5] * 10488027524736910128UL) + ((uint64_t)op[6] * 12457959752871040234UL) + ((uint64_t)op[7] * 11408867784380029229UL) + ((uint64_t)op[8] * 5013678683758370127UL) + ((uint64_t)op[9] * 16125149467419237595UL) + ((uint64_t)op[10] * 222469143476351708UL) + ((uint64_t)op[11] * 12667832759597063407UL)) * 14);
	tmp_q[5] = ((uint64_t)op[0] * 12667832759597063407UL) + ((uint64_t)op[1] * 10718469084878803840UL) + ((uint64_t)op[2] * 4929736269080567192UL) + ((uint64_t)op[3] * 2457931586228274314UL) + ((uint64_t)op[4] * 6110991910097450496UL) + ((uint64_t)op[5] * 2065473124114430057UL) + ((((uint64_t)op[6] * 10488027524736910128UL) + ((uint64_t)op[7] * 12457959752871040234UL) + ((uint64_t)op[8] * 11408867784380029229UL) + ((uint64_t)op[9] * 5013678683758370127UL) + ((uint64_t)op[10] * 16125149467419237595UL) + ((uint64_t)op[11] * 222469143476351708UL)) * 14);
	tmp_q[6] = ((uint64_t)op[0] * 222469143476351708UL) + ((uint64_t)op[1] * 12667832759597063407UL) + ((uint64_t)op[2] * 10718469084878803840UL) + ((uint64_t)op[3] * 4929736269080567192UL) + ((uint64_t)op[4] * 2457931586228274314UL) + ((uint64_t)op[5] * 6110991910097450496UL) + ((uint64_t)op[6] * 2065473124114430057UL) + ((((uint64_t)op[7] * 10488027524736910128UL) + ((uint64_t)op[8] * 12457959752871040234UL) + ((uint64_t)op[9] * 11408867784380029229UL) + ((uint64_t)op[10] * 5013678683758370127UL) + ((uint64_t)op[11] * 16125149467419237595UL)) * 14);
	tmp_q[7] = ((uint64_t)op[0] * 16125149467419237595UL) + ((uint64_t)op[1] * 222469143476351708UL) + ((uint64_t)op[2] * 12667832759597063407UL) + ((uint64_t)op[3] * 10718469084878803840UL) + ((uint64_t)op[4] * 4929736269080567192UL) + ((uint64_t)op[5] * 2457931586228274314UL) + ((uint64_t)op[6] * 6110991910097450496UL) + ((uint64_t)op[7] * 2065473124114430057UL) + ((((uint64_t)op[8] * 10488027524736910128UL) + ((uint64_t)op[9] * 12457959752871040234UL) + ((uint64_t)op[10] * 11408867784380029229UL) + ((uint64_t)op[11] * 5013678683758370127UL)) * 14);
	tmp_q[8] = ((uint64_t)op[0] * 5013678683758370127UL) + ((uint64_t)op[1] * 16125149467419237595UL) + ((uint64_t)op[2] * 222469143476351708UL) + ((uint64_t)op[3] * 12667832759597063407UL) + ((uint64_t)op[4] * 10718469084878803840UL) + ((uint64_t)op[5] * 4929736269080567192UL) + ((uint64_t)op[6] * 2457931586228274314UL) + ((uint64_t)op[7] * 6110991910097450496UL) + ((uint64_t)op[8] * 2065473124114430057UL) + ((((uint64_t)op[9] * 10488027524736910128UL) + ((uint64_t)op[10] * 12457959752871040234UL) + ((uint64_t)op[11] * 11408867784380029229UL)) * 14);
	tmp_q[9] = ((uint64_t)op[0] * 11408867784380029229UL) + ((uint64_t)op[1] * 5013678683758370127UL) + ((uint64_t)op[2] * 16125149467419237595UL) + ((uint64_t)op[3] * 222469143476351708UL) + ((uint64_t)op[4] * 12667832759597063407UL) + ((uint64_t)op[5] * 10718469084878803840UL) + ((uint64_t)op[6] * 4929736269080567192UL) + ((uint64_t)op[7] * 2457931586228274314UL) + ((uint64_t)op[8] * 6110991910097450496UL) + ((uint64_t)op[9] * 2065473124114430057UL) + ((((uint64_t)op[10] * 10488027524736910128UL) + ((uint64_t)op[11] * 12457959752871040234UL)) * 14);
	tmp_q[10] = ((uint64_t)op[0] * 12457959752871040234UL) + ((uint64_t)op[1] * 11408867784380029229UL) + ((uint64_t)op[2] * 5013678683758370127UL) + ((uint64_t)op[3] * 16125149467419237595UL) + ((uint64_t)op[4] * 222469143476351708UL) + ((uint64_t)op[5] * 12667832759597063407UL) + ((uint64_t)op[6] * 10718469084878803840UL) + ((uint64_t)op[7] * 4929736269080567192UL) + ((uint64_t)op[8] * 2457931586228274314UL) + ((uint64_t)op[9] * 6110991910097450496UL) + ((uint64_t)op[10] * 2065473124114430057UL) + ((uint64_t)op[11] * 17705176830349880480UL);
	tmp_q[11] = ((uint64_t)op[0] * 10488027524736910128UL) + ((uint64_t)op[1] * 12457959752871040234UL) + ((uint64_t)op[2] * 11408867784380029229UL) + ((uint64_t)op[3] * 5013678683758370127UL) + ((uint64_t)op[4] * 16125149467419237595UL) + ((uint64_t)op[5] * 222469143476351708UL) + ((uint64_t)op[6] * 12667832759597063407UL) + ((uint64_t)op[7] * 10718469084878803840UL) + ((uint64_t)op[8] * 4929736269080567192UL) + ((uint64_t)op[9] * 2457931586228274314UL) + ((uint64_t)op[10] * 6110991910097450496UL) + ((uint64_t)op[11] * 2065473124114430057UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4366695850131L) + ((-((int128)tmp_q[1] * 1053978248386L) - ((int128)tmp_q[2] * 1411074960785L) - ((int128)tmp_q[3] * 2613831889171L) - ((int128)tmp_q[4] * 5063450455047L) + ((int128)tmp_q[5] * 3361270032191L) + ((int128)tmp_q[6] * 4092500033260L) - ((int128)tmp_q[7] * 6777175800315L) + ((int128)tmp_q[8] * 872349556534L) - ((int128)tmp_q[9] * 3463534281386L) + ((int128)tmp_q[10] * 7797799032924L) - ((int128)tmp_q[11] * 6346060870288L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 6346060870288L) + ((int128)tmp_q[1] * 4366695850131L) + ((-((int128)tmp_q[2] * 1053978248386L) - ((int128)tmp_q[3] * 1411074960785L) - ((int128)tmp_q[4] * 2613831889171L) - ((int128)tmp_q[5] * 5063450455047L) + ((int128)tmp_q[6] * 3361270032191L) + ((int128)tmp_q[7] * 4092500033260L) - ((int128)tmp_q[8] * 6777175800315L) + ((int128)tmp_q[9] * 872349556534L) - ((int128)tmp_q[10] * 3463534281386L) + ((int128)tmp_q[11] * 7797799032924L)) * 14);
	tmp_zero[2] = ((int128)tmp_q[0] * 7797799032924L) - ((int128)tmp_q[1] * 6346060870288L) + ((int128)tmp_q[2] * 4366695850131L) + ((-((int128)tmp_q[3] * 1053978248386L) - ((int128)tmp_q[4] * 1411074960785L) - ((int128)tmp_q[5] * 2613831889171L) - ((int128)tmp_q[6] * 5063450455047L) + ((int128)tmp_q[7] * 3361270032191L) + ((int128)tmp_q[8] * 4092500033260L) - ((int128)tmp_q[9] * 6777175800315L) + ((int128)tmp_q[10] * 872349556534L) - ((int128)tmp_q[11] * 3463534281386L)) * 14);
	tmp_zero[3] = -((int128)tmp_q[0] * 3463534281386L) + ((int128)tmp_q[1] * 7797799032924L) - ((int128)tmp_q[2] * 6346060870288L) + ((int128)tmp_q[3] * 4366695850131L) + ((-((int128)tmp_q[4] * 1053978248386L) - ((int128)tmp_q[5] * 1411074960785L) - ((int128)tmp_q[6] * 2613831889171L) - ((int128)tmp_q[7] * 5063450455047L) + ((int128)tmp_q[8] * 3361270032191L) + ((int128)tmp_q[9] * 4092500033260L) - ((int128)tmp_q[10] * 6777175800315L) + ((int128)tmp_q[11] * 872349556534L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 872349556534L) - ((int128)tmp_q[1] * 3463534281386L) + ((int128)tmp_q[2] * 7797799032924L) - ((int128)tmp_q[3] * 6346060870288L) + ((int128)tmp_q[4] * 4366695850131L) + ((-((int128)tmp_q[5] * 1053978248386L) - ((int128)tmp_q[6] * 1411074960785L) - ((int128)tmp_q[7] * 2613831889171L) - ((int128)tmp_q[8] * 5063450455047L) + ((int128)tmp_q[9] * 3361270032191L) + ((int128)tmp_q[10] * 4092500033260L) - ((int128)tmp_q[11] * 6777175800315L)) * 14);
	tmp_zero[5] = -((int128)tmp_q[0] * 6777175800315L) + ((int128)tmp_q[1] * 872349556534L) - ((int128)tmp_q[2] * 3463534281386L) + ((int128)tmp_q[3] * 7797799032924L) - ((int128)tmp_q[4] * 6346060870288L) + ((int128)tmp_q[5] * 4366695850131L) + ((-((int128)tmp_q[6] * 1053978248386L) - ((int128)tmp_q[7] * 1411074960785L) - ((int128)tmp_q[8] * 2613831889171L) - ((int128)tmp_q[9] * 5063450455047L) + ((int128)tmp_q[10] * 3361270032191L) + ((int128)tmp_q[11] * 4092500033260L)) * 14);
	tmp_zero[6] = ((int128)tmp_q[0] * 4092500033260L) - ((int128)tmp_q[1] * 6777175800315L) + ((int128)tmp_q[2] * 872349556534L) - ((int128)tmp_q[3] * 3463534281386L) + ((int128)tmp_q[4] * 7797799032924L) - ((int128)tmp_q[5] * 6346060870288L) + ((int128)tmp_q[6] * 4366695850131L) + ((-((int128)tmp_q[7] * 1053978248386L) - ((int128)tmp_q[8] * 1411074960785L) - ((int128)tmp_q[9] * 2613831889171L) - ((int128)tmp_q[10] * 5063450455047L) + ((int128)tmp_q[11] * 3361270032191L)) * 14);
	tmp_zero[7] = ((int128)tmp_q[0] * 3361270032191L) + ((int128)tmp_q[1] * 4092500033260L) - ((int128)tmp_q[2] * 6777175800315L) + ((int128)tmp_q[3] * 872349556534L) - ((int128)tmp_q[4] * 3463534281386L) + ((int128)tmp_q[5] * 7797799032924L) - ((int128)tmp_q[6] * 6346060870288L) + ((int128)tmp_q[7] * 4366695850131L) + ((-((int128)tmp_q[8] * 1053978248386L) - ((int128)tmp_q[9] * 1411074960785L) - ((int128)tmp_q[10] * 2613831889171L) - ((int128)tmp_q[11] * 5063450455047L)) * 14);
	tmp_zero[8] = -((int128)tmp_q[0] * 5063450455047L) + ((int128)tmp_q[1] * 3361270032191L) + ((int128)tmp_q[2] * 4092500033260L) - ((int128)tmp_q[3] * 6777175800315L) + ((int128)tmp_q[4] * 872349556534L) - ((int128)tmp_q[5] * 3463534281386L) + ((int128)tmp_q[6] * 7797799032924L) - ((int128)tmp_q[7] * 6346060870288L) + ((int128)tmp_q[8] * 4366695850131L) + ((-((int128)tmp_q[9] * 1053978248386L) - ((int128)tmp_q[10] * 1411074960785L) - ((int128)tmp_q[11] * 2613831889171L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 2613831889171L) - ((int128)tmp_q[1] * 5063450455047L) + ((int128)tmp_q[2] * 3361270032191L) + ((int128)tmp_q[3] * 4092500033260L) - ((int128)tmp_q[4] * 6777175800315L) + ((int128)tmp_q[5] * 872349556534L) - ((int128)tmp_q[6] * 3463534281386L) + ((int128)tmp_q[7] * 7797799032924L) - ((int128)tmp_q[8] * 6346060870288L) + ((int128)tmp_q[9] * 4366695850131L) + ((-((int128)tmp_q[10] * 1053978248386L) - ((int128)tmp_q[11] * 1411074960785L)) * 14);
	tmp_zero[10] = -((int128)tmp_q[0] * 1411074960785L) - ((int128)tmp_q[1] * 2613831889171L) - ((int128)tmp_q[2] * 5063450455047L) + ((int128)tmp_q[3] * 3361270032191L) + ((int128)tmp_q[4] * 4092500033260L) - ((int128)tmp_q[5] * 6777175800315L) + ((int128)tmp_q[6] * 872349556534L) - ((int128)tmp_q[7] * 3463534281386L) + ((int128)tmp_q[8] * 7797799032924L) - ((int128)tmp_q[9] * 6346060870288L) + ((int128)tmp_q[10] * 4366695850131L) - ((int128)tmp_q[11] * 14755695477404L);
	tmp_zero[11] = -((int128)tmp_q[0] * 1053978248386L) - ((int128)tmp_q[1] * 1411074960785L) - ((int128)tmp_q[2] * 2613831889171L) - ((int128)tmp_q[3] * 5063450455047L) + ((int128)tmp_q[4] * 3361270032191L) + ((int128)tmp_q[5] * 4092500033260L) - ((int128)tmp_q[6] * 6777175800315L) + ((int128)tmp_q[7] * 872349556534L) - ((int128)tmp_q[8] * 3463534281386L) + ((int128)tmp_q[9] * 7797799032924L) - ((int128)tmp_q[10] * 6346060870288L) + ((int128)tmp_q[11] * 4366695850131L);

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

