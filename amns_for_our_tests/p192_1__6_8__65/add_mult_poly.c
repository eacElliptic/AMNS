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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17654853304438899917UL) + ((((uint64_t)op[1] * 16860802650017784739UL) + ((uint64_t)op[2] * 4481289927874956161UL) + ((uint64_t)op[3] * 18137323254876948952UL) + ((uint64_t)op[4] * 14221329960637640328UL) + ((uint64_t)op[5] * 13486995788705491102UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 13486995788705491102UL) + ((uint64_t)op[1] * 17654853304438899917UL) + ((((uint64_t)op[2] * 16860802650017784739UL) + ((uint64_t)op[3] * 4481289927874956161UL) + ((uint64_t)op[4] * 18137323254876948952UL) + ((uint64_t)op[5] * 14221329960637640328UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 14221329960637640328UL) + ((uint64_t)op[1] * 13486995788705491102UL) + ((uint64_t)op[2] * 17654853304438899917UL) + ((((uint64_t)op[3] * 16860802650017784739UL) + ((uint64_t)op[4] * 4481289927874956161UL) + ((uint64_t)op[5] * 18137323254876948952UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 18137323254876948952UL) + ((uint64_t)op[1] * 14221329960637640328UL) + ((uint64_t)op[2] * 13486995788705491102UL) + ((uint64_t)op[3] * 17654853304438899917UL) + ((((uint64_t)op[4] * 16860802650017784739UL) + ((uint64_t)op[5] * 4481289927874956161UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 4481289927874956161UL) + ((uint64_t)op[1] * 18137323254876948952UL) + ((uint64_t)op[2] * 14221329960637640328UL) + ((uint64_t)op[3] * 13486995788705491102UL) + ((uint64_t)op[4] * 17654853304438899917UL) + ((uint64_t)op[5] * 5759212684175416600UL);
	tmp_q[5] = ((uint64_t)op[0] * 16860802650017784739UL) + ((uint64_t)op[1] * 4481289927874956161UL) + ((uint64_t)op[2] * 18137323254876948952UL) + ((uint64_t)op[3] * 14221329960637640328UL) + ((uint64_t)op[4] * 13486995788705491102UL) + ((uint64_t)op[5] * 17654853304438899917UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 21465198499067L) + ((-((int128)tmp_q[1] * 4657271203281L) + ((int128)tmp_q[2] * 30522210317505L) + ((int128)tmp_q[3] * 52133076382544L) - ((int128)tmp_q[4] * 58191502156884L) - ((int128)tmp_q[5] * 24902220390450L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 24902220390450L) + ((int128)tmp_q[1] * 21465198499067L) + ((-((int128)tmp_q[2] * 4657271203281L) + ((int128)tmp_q[3] * 30522210317505L) + ((int128)tmp_q[4] * 52133076382544L) - ((int128)tmp_q[5] * 58191502156884L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 58191502156884L) - ((int128)tmp_q[1] * 24902220390450L) + ((int128)tmp_q[2] * 21465198499067L) + ((-((int128)tmp_q[3] * 4657271203281L) + ((int128)tmp_q[4] * 30522210317505L) + ((int128)tmp_q[5] * 52133076382544L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 52133076382544L) - ((int128)tmp_q[1] * 58191502156884L) - ((int128)tmp_q[2] * 24902220390450L) + ((int128)tmp_q[3] * 21465198499067L) + ((-((int128)tmp_q[4] * 4657271203281L) + ((int128)tmp_q[5] * 30522210317505L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 30522210317505L) + ((int128)tmp_q[1] * 52133076382544L) - ((int128)tmp_q[2] * 58191502156884L) - ((int128)tmp_q[3] * 24902220390450L) + ((int128)tmp_q[4] * 21465198499067L) - ((int128)tmp_q[5] * 37258169626248L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4657271203281L) + ((int128)tmp_q[1] * 30522210317505L) + ((int128)tmp_q[2] * 52133076382544L) - ((int128)tmp_q[3] * 58191502156884L) - ((int128)tmp_q[4] * 24902220390450L) + ((int128)tmp_q[5] * 21465198499067L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

