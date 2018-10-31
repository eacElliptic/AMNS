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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6149340217972058864UL) + ((((uint64_t)op[1] * 14077083719797523868UL) + ((uint64_t)op[2] * 7595537026148688443UL) + ((uint64_t)op[3] * 11817531521354123920UL) + ((uint64_t)op[4] * 17463456366273860722UL) + ((uint64_t)op[5] * 12029707225318025862UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 12029707225318025862UL) + ((uint64_t)op[1] * 6149340217972058864UL) + ((((uint64_t)op[2] * 14077083719797523868UL) + ((uint64_t)op[3] * 7595537026148688443UL) + ((uint64_t)op[4] * 11817531521354123920UL) + ((uint64_t)op[5] * 17463456366273860722UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 17463456366273860722UL) + ((uint64_t)op[1] * 12029707225318025862UL) + ((uint64_t)op[2] * 6149340217972058864UL) + ((((uint64_t)op[3] * 14077083719797523868UL) + ((uint64_t)op[4] * 7595537026148688443UL) + ((uint64_t)op[5] * 11817531521354123920UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 11817531521354123920UL) + ((uint64_t)op[1] * 17463456366273860722UL) + ((uint64_t)op[2] * 12029707225318025862UL) + ((uint64_t)op[3] * 6149340217972058864UL) + ((((uint64_t)op[4] * 14077083719797523868UL) + ((uint64_t)op[5] * 7595537026148688443UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 7595537026148688443UL) + ((uint64_t)op[1] * 11817531521354123920UL) + ((uint64_t)op[2] * 17463456366273860722UL) + ((uint64_t)op[3] * 12029707225318025862UL) + ((uint64_t)op[4] * 6149340217972058864UL) + ((uint64_t)op[5] * 3401557695850587124UL);
	tmp_q[5] = ((uint64_t)op[0] * 14077083719797523868UL) + ((uint64_t)op[1] * 7595537026148688443UL) + ((uint64_t)op[2] * 11817531521354123920UL) + ((uint64_t)op[3] * 17463456366273860722UL) + ((uint64_t)op[4] * 12029707225318025862UL) + ((uint64_t)op[5] * 6149340217972058864UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1501996599586L) - ((-((int128)tmp_q[1] * 3977954070498L) - ((int128)tmp_q[2] * 391367868108L) + ((int128)tmp_q[3] * 2993005574364L) - ((int128)tmp_q[4] * 1383579711933L) + ((int128)tmp_q[5] * 1409337735080L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1409337735080L) - ((int128)tmp_q[1] * 1501996599586L) - ((-((int128)tmp_q[2] * 3977954070498L) - ((int128)tmp_q[3] * 391367868108L) + ((int128)tmp_q[4] * 2993005574364L) - ((int128)tmp_q[5] * 1383579711933L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 1383579711933L) + ((int128)tmp_q[1] * 1409337735080L) - ((int128)tmp_q[2] * 1501996599586L) - ((-((int128)tmp_q[3] * 3977954070498L) - ((int128)tmp_q[4] * 391367868108L) + ((int128)tmp_q[5] * 2993005574364L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 2993005574364L) - ((int128)tmp_q[1] * 1383579711933L) + ((int128)tmp_q[2] * 1409337735080L) - ((int128)tmp_q[3] * 1501996599586L) - ((-((int128)tmp_q[4] * 3977954070498L) - ((int128)tmp_q[5] * 391367868108L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 391367868108L) + ((int128)tmp_q[1] * 2993005574364L) - ((int128)tmp_q[2] * 1383579711933L) + ((int128)tmp_q[3] * 1409337735080L) - ((int128)tmp_q[4] * 1501996599586L) + ((int128)tmp_q[5] * 19889770352490L);
	tmp_zero[5] = -((int128)tmp_q[0] * 3977954070498L) - ((int128)tmp_q[1] * 391367868108L) + ((int128)tmp_q[2] * 2993005574364L) - ((int128)tmp_q[3] * 1383579711933L) + ((int128)tmp_q[4] * 1409337735080L) - ((int128)tmp_q[5] * 1501996599586L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

