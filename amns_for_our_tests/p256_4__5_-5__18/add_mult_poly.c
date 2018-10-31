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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3547827637207658505UL) + ((((uint64_t)op[1] * 10320796016895631575UL) + ((uint64_t)op[2] * 4303665599163437800UL) + ((uint64_t)op[3] * 10251702371212028380UL) + ((uint64_t)op[4] * 1788688970983450585UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 1788688970983450585UL) + ((uint64_t)op[1] * 3547827637207658505UL) + ((((uint64_t)op[2] * 10320796016895631575UL) + ((uint64_t)op[3] * 4303665599163437800UL) + ((uint64_t)op[4] * 10251702371212028380UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 10251702371212028380UL) + ((uint64_t)op[1] * 1788688970983450585UL) + ((uint64_t)op[2] * 3547827637207658505UL) + ((((uint64_t)op[3] * 10320796016895631575UL) + ((uint64_t)op[4] * 4303665599163437800UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 4303665599163437800UL) + ((uint64_t)op[1] * 10251702371212028380UL) + ((uint64_t)op[2] * 1788688970983450585UL) + ((uint64_t)op[3] * 3547827637207658505UL) + ((uint64_t)op[4] * 3736252136650496973UL);
	tmp_q[4] = ((uint64_t)op[0] * 10320796016895631575UL) + ((uint64_t)op[1] * 4303665599163437800UL) + ((uint64_t)op[2] * 10251702371212028380UL) + ((uint64_t)op[3] * 1788688970983450585UL) + ((uint64_t)op[4] * 3547827637207658505UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 819155610235339L) - ((((int128)tmp_q[1] * 97450118623816L) - ((int128)tmp_q[2] * 696345306444905L) - ((int128)tmp_q[3] * 840864531696775L) + ((int128)tmp_q[4] * 1465302440019360L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1465302440019360L) + ((int128)tmp_q[1] * 819155610235339L) - ((((int128)tmp_q[2] * 97450118623816L) - ((int128)tmp_q[3] * 696345306444905L) - ((int128)tmp_q[4] * 840864531696775L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 840864531696775L) + ((int128)tmp_q[1] * 1465302440019360L) + ((int128)tmp_q[2] * 819155610235339L) - ((((int128)tmp_q[3] * 97450118623816L) - ((int128)tmp_q[4] * 696345306444905L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 696345306444905L) - ((int128)tmp_q[1] * 840864531696775L) + ((int128)tmp_q[2] * 1465302440019360L) + ((int128)tmp_q[3] * 819155610235339L) - ((int128)tmp_q[4] * 487250593119080L);
	tmp_zero[4] = ((int128)tmp_q[0] * 97450118623816L) - ((int128)tmp_q[1] * 696345306444905L) - ((int128)tmp_q[2] * 840864531696775L) + ((int128)tmp_q[3] * 1465302440019360L) + ((int128)tmp_q[4] * 819155610235339L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

