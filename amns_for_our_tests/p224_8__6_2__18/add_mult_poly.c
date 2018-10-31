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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14727889283387733275UL) + ((((uint64_t)op[1] * 17823868372462604731UL) + ((uint64_t)op[2] * 8338435313615136790UL) + ((uint64_t)op[3] * 3974250347261750107UL) + ((uint64_t)op[4] * 7678759986880546833UL) + ((uint64_t)op[5] * 18053231171055650487UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 18053231171055650487UL) + ((uint64_t)op[1] * 14727889283387733275UL) + ((((uint64_t)op[2] * 17823868372462604731UL) + ((uint64_t)op[3] * 8338435313615136790UL) + ((uint64_t)op[4] * 3974250347261750107UL) + ((uint64_t)op[5] * 7678759986880546833UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 7678759986880546833UL) + ((uint64_t)op[1] * 18053231171055650487UL) + ((uint64_t)op[2] * 14727889283387733275UL) + ((((uint64_t)op[3] * 17823868372462604731UL) + ((uint64_t)op[4] * 8338435313615136790UL) + ((uint64_t)op[5] * 3974250347261750107UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 3974250347261750107UL) + ((uint64_t)op[1] * 7678759986880546833UL) + ((uint64_t)op[2] * 18053231171055650487UL) + ((uint64_t)op[3] * 14727889283387733275UL) + ((((uint64_t)op[4] * 17823868372462604731UL) + ((uint64_t)op[5] * 8338435313615136790UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 8338435313615136790UL) + ((uint64_t)op[1] * 3974250347261750107UL) + ((uint64_t)op[2] * 7678759986880546833UL) + ((uint64_t)op[3] * 18053231171055650487UL) + ((uint64_t)op[4] * 14727889283387733275UL) + ((uint64_t)op[5] * 17200992671215657846UL);
	tmp_q[5] = ((uint64_t)op[0] * 17823868372462604731UL) + ((uint64_t)op[1] * 8338435313615136790UL) + ((uint64_t)op[2] * 3974250347261750107UL) + ((uint64_t)op[3] * 7678759986880546833UL) + ((uint64_t)op[4] * 18053231171055650487UL) + ((uint64_t)op[5] * 14727889283387733275UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 139921679163L) + ((-((int128)tmp_q[1] * 10805067012L) - ((int128)tmp_q[2] * 35561910389L) - ((int128)tmp_q[3] * 18186473114L) + ((int128)tmp_q[4] * 9201477620L) + ((int128)tmp_q[5] * 28320374745L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 28320374745L) - ((int128)tmp_q[1] * 139921679163L) + ((-((int128)tmp_q[2] * 10805067012L) - ((int128)tmp_q[3] * 35561910389L) - ((int128)tmp_q[4] * 18186473114L) + ((int128)tmp_q[5] * 9201477620L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 9201477620L) + ((int128)tmp_q[1] * 28320374745L) - ((int128)tmp_q[2] * 139921679163L) + ((-((int128)tmp_q[3] * 10805067012L) - ((int128)tmp_q[4] * 35561910389L) - ((int128)tmp_q[5] * 18186473114L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 18186473114L) + ((int128)tmp_q[1] * 9201477620L) + ((int128)tmp_q[2] * 28320374745L) - ((int128)tmp_q[3] * 139921679163L) + ((-((int128)tmp_q[4] * 10805067012L) - ((int128)tmp_q[5] * 35561910389L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 35561910389L) - ((int128)tmp_q[1] * 18186473114L) + ((int128)tmp_q[2] * 9201477620L) + ((int128)tmp_q[3] * 28320374745L) - ((int128)tmp_q[4] * 139921679163L) - ((int128)tmp_q[5] * 21610134024L);
	tmp_zero[5] = -((int128)tmp_q[0] * 10805067012L) - ((int128)tmp_q[1] * 35561910389L) - ((int128)tmp_q[2] * 18186473114L) + ((int128)tmp_q[3] * 9201477620L) + ((int128)tmp_q[4] * 28320374745L) - ((int128)tmp_q[5] * 139921679163L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

