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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[8] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4264517792475626133UL) + ((((uint64_t)op[1] * 12196859527393994929UL) + ((uint64_t)op[2] * 8516492286996961033UL) + ((uint64_t)op[3] * 10667286426085857484UL) + ((uint64_t)op[4] * 12260574400394934749UL) + ((uint64_t)op[5] * 6099854023424131419UL) + ((uint64_t)op[6] * 8166338908514903703UL) + ((uint64_t)op[7] * 4990480991241087354UL) + ((uint64_t)op[8] * 1635460715821022011UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 1635460715821022011UL) + ((uint64_t)op[1] * 4264517792475626133UL) + ((((uint64_t)op[2] * 12196859527393994929UL) + ((uint64_t)op[3] * 8516492286996961033UL) + ((uint64_t)op[4] * 10667286426085857484UL) + ((uint64_t)op[5] * 12260574400394934749UL) + ((uint64_t)op[6] * 6099854023424131419UL) + ((uint64_t)op[7] * 8166338908514903703UL) + ((uint64_t)op[8] * 4990480991241087354UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 4990480991241087354UL) + ((uint64_t)op[1] * 1635460715821022011UL) + ((uint64_t)op[2] * 4264517792475626133UL) + ((((uint64_t)op[3] * 12196859527393994929UL) + ((uint64_t)op[4] * 8516492286996961033UL) + ((uint64_t)op[5] * 10667286426085857484UL) + ((uint64_t)op[6] * 12260574400394934749UL) + ((uint64_t)op[7] * 6099854023424131419UL) + ((uint64_t)op[8] * 8166338908514903703UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 8166338908514903703UL) + ((uint64_t)op[1] * 4990480991241087354UL) + ((uint64_t)op[2] * 1635460715821022011UL) + ((uint64_t)op[3] * 4264517792475626133UL) + ((((uint64_t)op[4] * 12196859527393994929UL) + ((uint64_t)op[5] * 8516492286996961033UL) + ((uint64_t)op[6] * 10667286426085857484UL) + ((uint64_t)op[7] * 12260574400394934749UL) + ((uint64_t)op[8] * 6099854023424131419UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 6099854023424131419UL) + ((uint64_t)op[1] * 8166338908514903703UL) + ((uint64_t)op[2] * 4990480991241087354UL) + ((uint64_t)op[3] * 1635460715821022011UL) + ((uint64_t)op[4] * 4264517792475626133UL) + ((((uint64_t)op[5] * 12196859527393994929UL) + ((uint64_t)op[6] * 8516492286996961033UL) + ((uint64_t)op[7] * 10667286426085857484UL) + ((uint64_t)op[8] * 12260574400394934749UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 12260574400394934749UL) + ((uint64_t)op[1] * 6099854023424131419UL) + ((uint64_t)op[2] * 8166338908514903703UL) + ((uint64_t)op[3] * 4990480991241087354UL) + ((uint64_t)op[4] * 1635460715821022011UL) + ((uint64_t)op[5] * 4264517792475626133UL) + ((((uint64_t)op[6] * 12196859527393994929UL) + ((uint64_t)op[7] * 8516492286996961033UL) + ((uint64_t)op[8] * 10667286426085857484UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 10667286426085857484UL) + ((uint64_t)op[1] * 12260574400394934749UL) + ((uint64_t)op[2] * 6099854023424131419UL) + ((uint64_t)op[3] * 8166338908514903703UL) + ((uint64_t)op[4] * 4990480991241087354UL) + ((uint64_t)op[5] * 1635460715821022011UL) + ((uint64_t)op[6] * 4264517792475626133UL) + ((((uint64_t)op[7] * 12196859527393994929UL) + ((uint64_t)op[8] * 8516492286996961033UL)) * 5);
	tmp_q[7] = ((uint64_t)op[0] * 8516492286996961033UL) + ((uint64_t)op[1] * 10667286426085857484UL) + ((uint64_t)op[2] * 12260574400394934749UL) + ((uint64_t)op[3] * 6099854023424131419UL) + ((uint64_t)op[4] * 8166338908514903703UL) + ((uint64_t)op[5] * 4990480991241087354UL) + ((uint64_t)op[6] * 1635460715821022011UL) + ((uint64_t)op[7] * 4264517792475626133UL) + ((uint64_t)op[8] * 5644065415841319797UL);
	tmp_q[8] = ((uint64_t)op[0] * 12196859527393994929UL) + ((uint64_t)op[1] * 8516492286996961033UL) + ((uint64_t)op[2] * 10667286426085857484UL) + ((uint64_t)op[3] * 12260574400394934749UL) + ((uint64_t)op[4] * 6099854023424131419UL) + ((uint64_t)op[5] * 8166338908514903703UL) + ((uint64_t)op[6] * 4990480991241087354UL) + ((uint64_t)op[7] * 1635460715821022011UL) + ((uint64_t)op[8] * 4264517792475626133UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 138182144818L) + ((-((int128)tmp_q[1] * 1227401239056L) - ((int128)tmp_q[2] * 3532182571245L) + ((int128)tmp_q[3] * 4517879666009L) - ((int128)tmp_q[4] * 930031572683L) - ((int128)tmp_q[5] * 628987692451L) - ((int128)tmp_q[6] * 1769453038419L) + ((int128)tmp_q[7] * 1949649700864L) + ((int128)tmp_q[8] * 1591654883080L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1591654883080L) + ((int128)tmp_q[1] * 138182144818L) + ((-((int128)tmp_q[2] * 1227401239056L) - ((int128)tmp_q[3] * 3532182571245L) + ((int128)tmp_q[4] * 4517879666009L) - ((int128)tmp_q[5] * 930031572683L) - ((int128)tmp_q[6] * 628987692451L) - ((int128)tmp_q[7] * 1769453038419L) + ((int128)tmp_q[8] * 1949649700864L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 1949649700864L) + ((int128)tmp_q[1] * 1591654883080L) + ((int128)tmp_q[2] * 138182144818L) + ((-((int128)tmp_q[3] * 1227401239056L) - ((int128)tmp_q[4] * 3532182571245L) + ((int128)tmp_q[5] * 4517879666009L) - ((int128)tmp_q[6] * 930031572683L) - ((int128)tmp_q[7] * 628987692451L) - ((int128)tmp_q[8] * 1769453038419L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 1769453038419L) + ((int128)tmp_q[1] * 1949649700864L) + ((int128)tmp_q[2] * 1591654883080L) + ((int128)tmp_q[3] * 138182144818L) + ((-((int128)tmp_q[4] * 1227401239056L) - ((int128)tmp_q[5] * 3532182571245L) + ((int128)tmp_q[6] * 4517879666009L) - ((int128)tmp_q[7] * 930031572683L) - ((int128)tmp_q[8] * 628987692451L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 628987692451L) - ((int128)tmp_q[1] * 1769453038419L) + ((int128)tmp_q[2] * 1949649700864L) + ((int128)tmp_q[3] * 1591654883080L) + ((int128)tmp_q[4] * 138182144818L) + ((-((int128)tmp_q[5] * 1227401239056L) - ((int128)tmp_q[6] * 3532182571245L) + ((int128)tmp_q[7] * 4517879666009L) - ((int128)tmp_q[8] * 930031572683L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 930031572683L) - ((int128)tmp_q[1] * 628987692451L) - ((int128)tmp_q[2] * 1769453038419L) + ((int128)tmp_q[3] * 1949649700864L) + ((int128)tmp_q[4] * 1591654883080L) + ((int128)tmp_q[5] * 138182144818L) + ((-((int128)tmp_q[6] * 1227401239056L) - ((int128)tmp_q[7] * 3532182571245L) + ((int128)tmp_q[8] * 4517879666009L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 4517879666009L) - ((int128)tmp_q[1] * 930031572683L) - ((int128)tmp_q[2] * 628987692451L) - ((int128)tmp_q[3] * 1769453038419L) + ((int128)tmp_q[4] * 1949649700864L) + ((int128)tmp_q[5] * 1591654883080L) + ((int128)tmp_q[6] * 138182144818L) + ((-((int128)tmp_q[7] * 1227401239056L) - ((int128)tmp_q[8] * 3532182571245L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 3532182571245L) + ((int128)tmp_q[1] * 4517879666009L) - ((int128)tmp_q[2] * 930031572683L) - ((int128)tmp_q[3] * 628987692451L) - ((int128)tmp_q[4] * 1769453038419L) + ((int128)tmp_q[5] * 1949649700864L) + ((int128)tmp_q[6] * 1591654883080L) + ((int128)tmp_q[7] * 138182144818L) - ((int128)tmp_q[8] * 6137006195280L);
	tmp_zero[8] = -((int128)tmp_q[0] * 1227401239056L) - ((int128)tmp_q[1] * 3532182571245L) + ((int128)tmp_q[2] * 4517879666009L) - ((int128)tmp_q[3] * 930031572683L) - ((int128)tmp_q[4] * 628987692451L) - ((int128)tmp_q[5] * 1769453038419L) + ((int128)tmp_q[6] * 1949649700864L) + ((int128)tmp_q[7] * 1591654883080L) + ((int128)tmp_q[8] * 138182144818L);

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

