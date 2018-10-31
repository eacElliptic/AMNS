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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15173511578223190832UL) + ((((uint64_t)op[1] * 10539566962992917254UL) + ((uint64_t)op[2] * 24996961155742319UL) + ((uint64_t)op[3] * 7804288293212756619UL) + ((uint64_t)op[4] * 3002451959037135620UL) + ((uint64_t)op[5] * 16202058555204142013UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16202058555204142013UL) + ((uint64_t)op[1] * 15173511578223190832UL) + ((((uint64_t)op[2] * 10539566962992917254UL) + ((uint64_t)op[3] * 24996961155742319UL) + ((uint64_t)op[4] * 7804288293212756619UL) + ((uint64_t)op[5] * 3002451959037135620UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 3002451959037135620UL) + ((uint64_t)op[1] * 16202058555204142013UL) + ((uint64_t)op[2] * 15173511578223190832UL) + ((((uint64_t)op[3] * 10539566962992917254UL) + ((uint64_t)op[4] * 24996961155742319UL) + ((uint64_t)op[5] * 7804288293212756619UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7804288293212756619UL) + ((uint64_t)op[1] * 3002451959037135620UL) + ((uint64_t)op[2] * 16202058555204142013UL) + ((uint64_t)op[3] * 15173511578223190832UL) + ((((uint64_t)op[4] * 10539566962992917254UL) + ((uint64_t)op[5] * 24996961155742319UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 24996961155742319UL) + ((uint64_t)op[1] * 7804288293212756619UL) + ((uint64_t)op[2] * 3002451959037135620UL) + ((uint64_t)op[3] * 16202058555204142013UL) + ((uint64_t)op[4] * 15173511578223190832UL) + ((uint64_t)op[5] * 5274787258440351470UL);
	tmp_q[5] = ((uint64_t)op[0] * 10539566962992917254UL) + ((uint64_t)op[1] * 24996961155742319UL) + ((uint64_t)op[2] * 7804288293212756619UL) + ((uint64_t)op[3] * 3002451959037135620UL) + ((uint64_t)op[4] * 16202058555204142013UL) + ((uint64_t)op[5] * 15173511578223190832UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 525594884L) - ((-((int128)tmp_q[1] * 571268248L) - ((int128)tmp_q[2] * 1970204273L) + ((int128)tmp_q[3] * 167651901L) - ((int128)tmp_q[4] * 1259814450L) + ((int128)tmp_q[5] * 1565114345L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 1565114345L) + ((int128)tmp_q[1] * 525594884L) - ((-((int128)tmp_q[2] * 571268248L) - ((int128)tmp_q[3] * 1970204273L) + ((int128)tmp_q[4] * 167651901L) - ((int128)tmp_q[5] * 1259814450L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1259814450L) + ((int128)tmp_q[1] * 1565114345L) + ((int128)tmp_q[2] * 525594884L) - ((-((int128)tmp_q[3] * 571268248L) - ((int128)tmp_q[4] * 1970204273L) + ((int128)tmp_q[5] * 167651901L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 167651901L) - ((int128)tmp_q[1] * 1259814450L) + ((int128)tmp_q[2] * 1565114345L) + ((int128)tmp_q[3] * 525594884L) - ((-((int128)tmp_q[4] * 571268248L) - ((int128)tmp_q[5] * 1970204273L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 1970204273L) + ((int128)tmp_q[1] * 167651901L) - ((int128)tmp_q[2] * 1259814450L) + ((int128)tmp_q[3] * 1565114345L) + ((int128)tmp_q[4] * 525594884L) + ((int128)tmp_q[5] * 1713804744L);
	tmp_zero[5] = -((int128)tmp_q[0] * 571268248L) - ((int128)tmp_q[1] * 1970204273L) + ((int128)tmp_q[2] * 167651901L) - ((int128)tmp_q[3] * 1259814450L) + ((int128)tmp_q[4] * 1565114345L) + ((int128)tmp_q[5] * 525594884L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

