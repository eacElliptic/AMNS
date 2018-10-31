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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12957636528698460687UL) + ((((uint64_t)op[1] * 10291059704563999853UL) + ((uint64_t)op[2] * 16769542762213218476UL) + ((uint64_t)op[3] * 8176236704410221776UL) + ((uint64_t)op[4] * 7799912804870791090UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 7799912804870791090UL) + ((uint64_t)op[1] * 12957636528698460687UL) + ((((uint64_t)op[2] * 10291059704563999853UL) + ((uint64_t)op[3] * 16769542762213218476UL) + ((uint64_t)op[4] * 8176236704410221776UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 8176236704410221776UL) + ((uint64_t)op[1] * 7799912804870791090UL) + ((uint64_t)op[2] * 12957636528698460687UL) + ((((uint64_t)op[3] * 10291059704563999853UL) + ((uint64_t)op[4] * 16769542762213218476UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 16769542762213218476UL) + ((uint64_t)op[1] * 8176236704410221776UL) + ((uint64_t)op[2] * 7799912804870791090UL) + ((uint64_t)op[3] * 12957636528698460687UL) + ((uint64_t)op[4] * 8541501341673792360UL);
	tmp_q[4] = ((uint64_t)op[0] * 10291059704563999853UL) + ((uint64_t)op[1] * 16769542762213218476UL) + ((uint64_t)op[2] * 8176236704410221776UL) + ((uint64_t)op[3] * 7799912804870791090UL) + ((uint64_t)op[4] * 12957636528698460687UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 26644764819151L) + ((-((int128)tmp_q[1] * 1730403686595L) + ((int128)tmp_q[2] * 15897241496956L) - ((int128)tmp_q[3] * 9564091297708L) + ((int128)tmp_q[4] * 15983457341778L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 15983457341778L) - ((int128)tmp_q[1] * 26644764819151L) + ((-((int128)tmp_q[2] * 1730403686595L) + ((int128)tmp_q[3] * 15897241496956L) - ((int128)tmp_q[4] * 9564091297708L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 9564091297708L) + ((int128)tmp_q[1] * 15983457341778L) - ((int128)tmp_q[2] * 26644764819151L) + ((-((int128)tmp_q[3] * 1730403686595L) + ((int128)tmp_q[4] * 15897241496956L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 15897241496956L) - ((int128)tmp_q[1] * 9564091297708L) + ((int128)tmp_q[2] * 15983457341778L) - ((int128)tmp_q[3] * 26644764819151L) - ((int128)tmp_q[4] * 13843229492760L);
	tmp_zero[4] = -((int128)tmp_q[0] * 1730403686595L) + ((int128)tmp_q[1] * 15897241496956L) - ((int128)tmp_q[2] * 9564091297708L) + ((int128)tmp_q[3] * 15983457341778L) - ((int128)tmp_q[4] * 26644764819151L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

