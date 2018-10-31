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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14950765889432979947UL) + ((((uint64_t)op[1] * 4127767336669811286UL) + ((uint64_t)op[2] * 13029272736891401933UL) + ((uint64_t)op[3] * 5712093235025719889UL) + ((uint64_t)op[4] * 11371157011405589591UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 11371157011405589591UL) + ((uint64_t)op[1] * 14950765889432979947UL) + ((((uint64_t)op[2] * 4127767336669811286UL) + ((uint64_t)op[3] * 13029272736891401933UL) + ((uint64_t)op[4] * 5712093235025719889UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 5712093235025719889UL) + ((uint64_t)op[1] * 11371157011405589591UL) + ((uint64_t)op[2] * 14950765889432979947UL) + ((((uint64_t)op[3] * 4127767336669811286UL) + ((uint64_t)op[4] * 13029272736891401933UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 13029272736891401933UL) + ((uint64_t)op[1] * 5712093235025719889UL) + ((uint64_t)op[2] * 11371157011405589591UL) + ((uint64_t)op[3] * 14950765889432979947UL) + ((uint64_t)op[4] * 6319859946309316100UL);
	tmp_q[4] = ((uint64_t)op[0] * 4127767336669811286UL) + ((uint64_t)op[1] * 13029272736891401933UL) + ((uint64_t)op[2] * 5712093235025719889UL) + ((uint64_t)op[3] * 11371157011405589591UL) + ((uint64_t)op[4] * 14950765889432979947UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 15722207331L) + ((-((int128)tmp_q[1] * 124208186165L) - ((int128)tmp_q[2] * 115566848228L) - ((int128)tmp_q[3] * 191658085566L) + ((int128)tmp_q[4] * 42151263995L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 42151263995L) + ((int128)tmp_q[1] * 15722207331L) + ((-((int128)tmp_q[2] * 124208186165L) - ((int128)tmp_q[3] * 115566848228L) - ((int128)tmp_q[4] * 191658085566L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 191658085566L) + ((int128)tmp_q[1] * 42151263995L) + ((int128)tmp_q[2] * 15722207331L) + ((-((int128)tmp_q[3] * 124208186165L) - ((int128)tmp_q[4] * 115566848228L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 115566848228L) - ((int128)tmp_q[1] * 191658085566L) + ((int128)tmp_q[2] * 42151263995L) + ((int128)tmp_q[3] * 15722207331L) - ((int128)tmp_q[4] * 745249116990L);
	tmp_zero[4] = -((int128)tmp_q[0] * 124208186165L) - ((int128)tmp_q[1] * 115566848228L) - ((int128)tmp_q[2] * 191658085566L) + ((int128)tmp_q[3] * 42151263995L) + ((int128)tmp_q[4] * 15722207331L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

