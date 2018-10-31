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
	tmp_q[0] = ((uint64_t)op[0] * 678407028740099093UL) + ((((uint64_t)op[1] * 1876566986423599654UL) + ((uint64_t)op[2] * 18247261468792710590UL) + ((uint64_t)op[3] * 6581418314320118754UL) + ((uint64_t)op[4] * 16474497495296353808UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 16474497495296353808UL) + ((uint64_t)op[1] * 678407028740099093UL) + ((((uint64_t)op[2] * 1876566986423599654UL) + ((uint64_t)op[3] * 18247261468792710590UL) + ((uint64_t)op[4] * 6581418314320118754UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 6581418314320118754UL) + ((uint64_t)op[1] * 16474497495296353808UL) + ((uint64_t)op[2] * 678407028740099093UL) + ((((uint64_t)op[3] * 1876566986423599654UL) + ((uint64_t)op[4] * 18247261468792710590UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 18247261468792710590UL) + ((uint64_t)op[1] * 6581418314320118754UL) + ((uint64_t)op[2] * 16474497495296353808UL) + ((uint64_t)op[3] * 678407028740099093UL) + ((uint64_t)op[4] * 11259401918541597924UL);
	tmp_q[4] = ((uint64_t)op[0] * 1876566986423599654UL) + ((uint64_t)op[1] * 18247261468792710590UL) + ((uint64_t)op[2] * 6581418314320118754UL) + ((uint64_t)op[3] * 16474497495296353808UL) + ((uint64_t)op[4] * 678407028740099093UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5247725613619L) + ((-((int128)tmp_q[1] * 6146201649934L) - ((int128)tmp_q[2] * 7510747559530L) - ((int128)tmp_q[3] * 3022594840366L) + ((int128)tmp_q[4] * 4485943553048L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 4485943553048L) + ((int128)tmp_q[1] * 5247725613619L) + ((-((int128)tmp_q[2] * 6146201649934L) - ((int128)tmp_q[3] * 7510747559530L) - ((int128)tmp_q[4] * 3022594840366L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 3022594840366L) + ((int128)tmp_q[1] * 4485943553048L) + ((int128)tmp_q[2] * 5247725613619L) + ((-((int128)tmp_q[3] * 6146201649934L) - ((int128)tmp_q[4] * 7510747559530L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 7510747559530L) - ((int128)tmp_q[1] * 3022594840366L) + ((int128)tmp_q[2] * 4485943553048L) + ((int128)tmp_q[3] * 5247725613619L) - ((int128)tmp_q[4] * 36877209899604L);
	tmp_zero[4] = -((int128)tmp_q[0] * 6146201649934L) - ((int128)tmp_q[1] * 7510747559530L) - ((int128)tmp_q[2] * 3022594840366L) + ((int128)tmp_q[3] * 4485943553048L) + ((int128)tmp_q[4] * 5247725613619L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

