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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7433702981718464447UL) + ((((uint64_t)op[1] * 427418805664619295UL) + ((uint64_t)op[2] * 11847357205148539210UL) + ((uint64_t)op[3] * 3495100076938877122UL) + ((uint64_t)op[4] * 4487746700900384781UL) + ((uint64_t)op[5] * 12001803729351731310UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 12001803729351731310UL) + ((uint64_t)op[1] * 7433702981718464447UL) + ((((uint64_t)op[2] * 427418805664619295UL) + ((uint64_t)op[3] * 11847357205148539210UL) + ((uint64_t)op[4] * 3495100076938877122UL) + ((uint64_t)op[5] * 4487746700900384781UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 4487746700900384781UL) + ((uint64_t)op[1] * 12001803729351731310UL) + ((uint64_t)op[2] * 7433702981718464447UL) + ((((uint64_t)op[3] * 427418805664619295UL) + ((uint64_t)op[4] * 11847357205148539210UL) + ((uint64_t)op[5] * 3495100076938877122UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 3495100076938877122UL) + ((uint64_t)op[1] * 4487746700900384781UL) + ((uint64_t)op[2] * 12001803729351731310UL) + ((uint64_t)op[3] * 7433702981718464447UL) + ((((uint64_t)op[4] * 427418805664619295UL) + ((uint64_t)op[5] * 11847357205148539210UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 11847357205148539210UL) + ((uint64_t)op[1] * 3495100076938877122UL) + ((uint64_t)op[2] * 4487746700900384781UL) + ((uint64_t)op[3] * 12001803729351731310UL) + ((uint64_t)op[4] * 7433702981718464447UL) + ((uint64_t)op[5] * 15454812434057216551UL);
	tmp_q[5] = ((uint64_t)op[0] * 427418805664619295UL) + ((uint64_t)op[1] * 11847357205148539210UL) + ((uint64_t)op[2] * 3495100076938877122UL) + ((uint64_t)op[3] * 4487746700900384781UL) + ((uint64_t)op[4] * 12001803729351731310UL) + ((uint64_t)op[5] * 7433702981718464447UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1217372239L) - ((((int128)tmp_q[1] * 2038371567L) + ((int128)tmp_q[2] * 769266160L) + ((int128)tmp_q[3] * 448153310L) + ((int128)tmp_q[4] * 1546541181L) - ((int128)tmp_q[5] * 3349198896L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 3349198896L) - ((int128)tmp_q[1] * 1217372239L) - ((((int128)tmp_q[2] * 2038371567L) + ((int128)tmp_q[3] * 769266160L) + ((int128)tmp_q[4] * 448153310L) + ((int128)tmp_q[5] * 1546541181L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 1546541181L) - ((int128)tmp_q[1] * 3349198896L) - ((int128)tmp_q[2] * 1217372239L) - ((((int128)tmp_q[3] * 2038371567L) + ((int128)tmp_q[4] * 769266160L) + ((int128)tmp_q[5] * 448153310L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 448153310L) + ((int128)tmp_q[1] * 1546541181L) - ((int128)tmp_q[2] * 3349198896L) - ((int128)tmp_q[3] * 1217372239L) - ((((int128)tmp_q[4] * 2038371567L) + ((int128)tmp_q[5] * 769266160L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 769266160L) + ((int128)tmp_q[1] * 448153310L) + ((int128)tmp_q[2] * 1546541181L) - ((int128)tmp_q[3] * 3349198896L) - ((int128)tmp_q[4] * 1217372239L) - ((int128)tmp_q[5] * 14268600969L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2038371567L) + ((int128)tmp_q[1] * 769266160L) + ((int128)tmp_q[2] * 448153310L) + ((int128)tmp_q[3] * 1546541181L) - ((int128)tmp_q[4] * 3349198896L) - ((int128)tmp_q[5] * 1217372239L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

