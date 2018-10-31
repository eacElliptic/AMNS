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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15365219625122286681UL) + ((((uint64_t)op[1] * 16180710607903497142UL) + ((uint64_t)op[2] * 17539810394739483238UL) + ((uint64_t)op[3] * 15323456766100265487UL) + ((uint64_t)op[4] * 2033630235560191047UL) + ((uint64_t)op[5] * 8581910008632968224UL) + ((uint64_t)op[6] * 1291443711998124290UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 1291443711998124290UL) + ((uint64_t)op[1] * 15365219625122286681UL) + ((((uint64_t)op[2] * 16180710607903497142UL) + ((uint64_t)op[3] * 17539810394739483238UL) + ((uint64_t)op[4] * 15323456766100265487UL) + ((uint64_t)op[5] * 2033630235560191047UL) + ((uint64_t)op[6] * 8581910008632968224UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 8581910008632968224UL) + ((uint64_t)op[1] * 1291443711998124290UL) + ((uint64_t)op[2] * 15365219625122286681UL) + ((((uint64_t)op[3] * 16180710607903497142UL) + ((uint64_t)op[4] * 17539810394739483238UL) + ((uint64_t)op[5] * 15323456766100265487UL) + ((uint64_t)op[6] * 2033630235560191047UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 2033630235560191047UL) + ((uint64_t)op[1] * 8581910008632968224UL) + ((uint64_t)op[2] * 1291443711998124290UL) + ((uint64_t)op[3] * 15365219625122286681UL) + ((((uint64_t)op[4] * 16180710607903497142UL) + ((uint64_t)op[5] * 17539810394739483238UL) + ((uint64_t)op[6] * 15323456766100265487UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 15323456766100265487UL) + ((uint64_t)op[1] * 2033630235560191047UL) + ((uint64_t)op[2] * 8581910008632968224UL) + ((uint64_t)op[3] * 1291443711998124290UL) + ((uint64_t)op[4] * 15365219625122286681UL) + ((((uint64_t)op[5] * 16180710607903497142UL) + ((uint64_t)op[6] * 17539810394739483238UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 17539810394739483238UL) + ((uint64_t)op[1] * 15323456766100265487UL) + ((uint64_t)op[2] * 2033630235560191047UL) + ((uint64_t)op[3] * 8581910008632968224UL) + ((uint64_t)op[4] * 1291443711998124290UL) + ((uint64_t)op[5] * 15365219625122286681UL) + ((uint64_t)op[6] * 11330167329030272370UL);
	tmp_q[6] = ((uint64_t)op[0] * 16180710607903497142UL) + ((uint64_t)op[1] * 17539810394739483238UL) + ((uint64_t)op[2] * 15323456766100265487UL) + ((uint64_t)op[3] * 2033630235560191047UL) + ((uint64_t)op[4] * 8581910008632968224UL) + ((uint64_t)op[5] * 1291443711998124290UL) + ((uint64_t)op[6] * 15365219625122286681UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 15318438337L) - ((-((int128)tmp_q[1] * 18679975054L) + ((int128)tmp_q[2] * 61049111L) - ((int128)tmp_q[3] * 20973037101L) - ((int128)tmp_q[4] * 23354802509L) + ((int128)tmp_q[5] * 47251249159L) - ((int128)tmp_q[6] * 41004092652L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 41004092652L) - ((int128)tmp_q[1] * 15318438337L) - ((-((int128)tmp_q[2] * 18679975054L) + ((int128)tmp_q[3] * 61049111L) - ((int128)tmp_q[4] * 20973037101L) - ((int128)tmp_q[5] * 23354802509L) + ((int128)tmp_q[6] * 47251249159L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 47251249159L) - ((int128)tmp_q[1] * 41004092652L) - ((int128)tmp_q[2] * 15318438337L) - ((-((int128)tmp_q[3] * 18679975054L) + ((int128)tmp_q[4] * 61049111L) - ((int128)tmp_q[5] * 20973037101L) - ((int128)tmp_q[6] * 23354802509L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 23354802509L) + ((int128)tmp_q[1] * 47251249159L) - ((int128)tmp_q[2] * 41004092652L) - ((int128)tmp_q[3] * 15318438337L) - ((-((int128)tmp_q[4] * 18679975054L) + ((int128)tmp_q[5] * 61049111L) - ((int128)tmp_q[6] * 20973037101L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 20973037101L) - ((int128)tmp_q[1] * 23354802509L) + ((int128)tmp_q[2] * 47251249159L) - ((int128)tmp_q[3] * 41004092652L) - ((int128)tmp_q[4] * 15318438337L) - ((-((int128)tmp_q[5] * 18679975054L) + ((int128)tmp_q[6] * 61049111L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 61049111L) - ((int128)tmp_q[1] * 20973037101L) - ((int128)tmp_q[2] * 23354802509L) + ((int128)tmp_q[3] * 47251249159L) - ((int128)tmp_q[4] * 41004092652L) - ((int128)tmp_q[5] * 15318438337L) + ((int128)tmp_q[6] * 93399875270L);
	tmp_zero[6] = -((int128)tmp_q[0] * 18679975054L) + ((int128)tmp_q[1] * 61049111L) - ((int128)tmp_q[2] * 20973037101L) - ((int128)tmp_q[3] * 23354802509L) + ((int128)tmp_q[4] * 47251249159L) - ((int128)tmp_q[5] * 41004092652L) - ((int128)tmp_q[6] * 15318438337L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

