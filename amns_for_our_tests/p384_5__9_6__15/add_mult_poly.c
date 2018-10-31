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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[8] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9311829605497953335UL) + ((((uint64_t)op[1] * 14716928100050876944UL) + ((uint64_t)op[2] * 8154556347845020470UL) + ((uint64_t)op[3] * 11940551012166874464UL) + ((uint64_t)op[4] * 297704423584373111UL) + ((uint64_t)op[5] * 2818149301471947856UL) + ((uint64_t)op[6] * 13094740998564192732UL) + ((uint64_t)op[7] * 11903080324992329428UL) + ((uint64_t)op[8] * 1415954847831113611UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 1415954847831113611UL) + ((uint64_t)op[1] * 9311829605497953335UL) + ((((uint64_t)op[2] * 14716928100050876944UL) + ((uint64_t)op[3] * 8154556347845020470UL) + ((uint64_t)op[4] * 11940551012166874464UL) + ((uint64_t)op[5] * 297704423584373111UL) + ((uint64_t)op[6] * 2818149301471947856UL) + ((uint64_t)op[7] * 13094740998564192732UL) + ((uint64_t)op[8] * 11903080324992329428UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 11903080324992329428UL) + ((uint64_t)op[1] * 1415954847831113611UL) + ((uint64_t)op[2] * 9311829605497953335UL) + ((((uint64_t)op[3] * 14716928100050876944UL) + ((uint64_t)op[4] * 8154556347845020470UL) + ((uint64_t)op[5] * 11940551012166874464UL) + ((uint64_t)op[6] * 297704423584373111UL) + ((uint64_t)op[7] * 2818149301471947856UL) + ((uint64_t)op[8] * 13094740998564192732UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 13094740998564192732UL) + ((uint64_t)op[1] * 11903080324992329428UL) + ((uint64_t)op[2] * 1415954847831113611UL) + ((uint64_t)op[3] * 9311829605497953335UL) + ((((uint64_t)op[4] * 14716928100050876944UL) + ((uint64_t)op[5] * 8154556347845020470UL) + ((uint64_t)op[6] * 11940551012166874464UL) + ((uint64_t)op[7] * 297704423584373111UL) + ((uint64_t)op[8] * 2818149301471947856UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 2818149301471947856UL) + ((uint64_t)op[1] * 13094740998564192732UL) + ((uint64_t)op[2] * 11903080324992329428UL) + ((uint64_t)op[3] * 1415954847831113611UL) + ((uint64_t)op[4] * 9311829605497953335UL) + ((((uint64_t)op[5] * 14716928100050876944UL) + ((uint64_t)op[6] * 8154556347845020470UL) + ((uint64_t)op[7] * 11940551012166874464UL) + ((uint64_t)op[8] * 297704423584373111UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 297704423584373111UL) + ((uint64_t)op[1] * 2818149301471947856UL) + ((uint64_t)op[2] * 13094740998564192732UL) + ((uint64_t)op[3] * 11903080324992329428UL) + ((uint64_t)op[4] * 1415954847831113611UL) + ((uint64_t)op[5] * 9311829605497953335UL) + ((((uint64_t)op[6] * 14716928100050876944UL) + ((uint64_t)op[7] * 8154556347845020470UL) + ((uint64_t)op[8] * 11940551012166874464UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 11940551012166874464UL) + ((uint64_t)op[1] * 297704423584373111UL) + ((uint64_t)op[2] * 2818149301471947856UL) + ((uint64_t)op[3] * 13094740998564192732UL) + ((uint64_t)op[4] * 11903080324992329428UL) + ((uint64_t)op[5] * 1415954847831113611UL) + ((uint64_t)op[6] * 9311829605497953335UL) + ((((uint64_t)op[7] * 14716928100050876944UL) + ((uint64_t)op[8] * 8154556347845020470UL)) * 6);
	tmp_q[7] = ((uint64_t)op[0] * 8154556347845020470UL) + ((uint64_t)op[1] * 11940551012166874464UL) + ((uint64_t)op[2] * 297704423584373111UL) + ((uint64_t)op[3] * 2818149301471947856UL) + ((uint64_t)op[4] * 13094740998564192732UL) + ((uint64_t)op[5] * 11903080324992329428UL) + ((uint64_t)op[6] * 1415954847831113611UL) + ((uint64_t)op[7] * 9311829605497953335UL) + ((uint64_t)op[8] * 14514592305467055200UL);
	tmp_q[8] = ((uint64_t)op[0] * 14716928100050876944UL) + ((uint64_t)op[1] * 8154556347845020470UL) + ((uint64_t)op[2] * 11940551012166874464UL) + ((uint64_t)op[3] * 297704423584373111UL) + ((uint64_t)op[4] * 2818149301471947856UL) + ((uint64_t)op[5] * 13094740998564192732UL) + ((uint64_t)op[6] * 11903080324992329428UL) + ((uint64_t)op[7] * 1415954847831113611UL) + ((uint64_t)op[8] * 9311829605497953335UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1585682667735L) + ((((int128)tmp_q[1] * 382371179829L) - ((int128)tmp_q[2] * 3364146668268L) - ((int128)tmp_q[3] * 106167936149L) + ((int128)tmp_q[4] * 2720449969278L) - ((int128)tmp_q[5] * 90808782203L) - ((int128)tmp_q[6] * 3705897509603L) - ((int128)tmp_q[7] * 1336405184505L) - ((int128)tmp_q[8] * 531557470205L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 531557470205L) - ((int128)tmp_q[1] * 1585682667735L) + ((((int128)tmp_q[2] * 382371179829L) - ((int128)tmp_q[3] * 3364146668268L) - ((int128)tmp_q[4] * 106167936149L) + ((int128)tmp_q[5] * 2720449969278L) - ((int128)tmp_q[6] * 90808782203L) - ((int128)tmp_q[7] * 3705897509603L) - ((int128)tmp_q[8] * 1336405184505L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 1336405184505L) - ((int128)tmp_q[1] * 531557470205L) - ((int128)tmp_q[2] * 1585682667735L) + ((((int128)tmp_q[3] * 382371179829L) - ((int128)tmp_q[4] * 3364146668268L) - ((int128)tmp_q[5] * 106167936149L) + ((int128)tmp_q[6] * 2720449969278L) - ((int128)tmp_q[7] * 90808782203L) - ((int128)tmp_q[8] * 3705897509603L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 3705897509603L) - ((int128)tmp_q[1] * 1336405184505L) - ((int128)tmp_q[2] * 531557470205L) - ((int128)tmp_q[3] * 1585682667735L) + ((((int128)tmp_q[4] * 382371179829L) - ((int128)tmp_q[5] * 3364146668268L) - ((int128)tmp_q[6] * 106167936149L) + ((int128)tmp_q[7] * 2720449969278L) - ((int128)tmp_q[8] * 90808782203L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 90808782203L) - ((int128)tmp_q[1] * 3705897509603L) - ((int128)tmp_q[2] * 1336405184505L) - ((int128)tmp_q[3] * 531557470205L) - ((int128)tmp_q[4] * 1585682667735L) + ((((int128)tmp_q[5] * 382371179829L) - ((int128)tmp_q[6] * 3364146668268L) - ((int128)tmp_q[7] * 106167936149L) + ((int128)tmp_q[8] * 2720449969278L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 2720449969278L) - ((int128)tmp_q[1] * 90808782203L) - ((int128)tmp_q[2] * 3705897509603L) - ((int128)tmp_q[3] * 1336405184505L) - ((int128)tmp_q[4] * 531557470205L) - ((int128)tmp_q[5] * 1585682667735L) + ((((int128)tmp_q[6] * 382371179829L) - ((int128)tmp_q[7] * 3364146668268L) - ((int128)tmp_q[8] * 106167936149L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 106167936149L) + ((int128)tmp_q[1] * 2720449969278L) - ((int128)tmp_q[2] * 90808782203L) - ((int128)tmp_q[3] * 3705897509603L) - ((int128)tmp_q[4] * 1336405184505L) - ((int128)tmp_q[5] * 531557470205L) - ((int128)tmp_q[6] * 1585682667735L) + ((((int128)tmp_q[7] * 382371179829L) - ((int128)tmp_q[8] * 3364146668268L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 3364146668268L) - ((int128)tmp_q[1] * 106167936149L) + ((int128)tmp_q[2] * 2720449969278L) - ((int128)tmp_q[3] * 90808782203L) - ((int128)tmp_q[4] * 3705897509603L) - ((int128)tmp_q[5] * 1336405184505L) - ((int128)tmp_q[6] * 531557470205L) - ((int128)tmp_q[7] * 1585682667735L) + ((int128)tmp_q[8] * 2294227078974L);
	tmp_zero[8] = ((int128)tmp_q[0] * 382371179829L) - ((int128)tmp_q[1] * 3364146668268L) - ((int128)tmp_q[2] * 106167936149L) + ((int128)tmp_q[3] * 2720449969278L) - ((int128)tmp_q[4] * 90808782203L) - ((int128)tmp_q[5] * 3705897509603L) - ((int128)tmp_q[6] * 1336405184505L) - ((int128)tmp_q[7] * 531557470205L) - ((int128)tmp_q[8] * 1585682667735L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
}

