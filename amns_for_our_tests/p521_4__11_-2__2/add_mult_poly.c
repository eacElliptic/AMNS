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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2206005580371624619UL) + ((((uint64_t)op[1] * 10623114139826810515UL) + ((uint64_t)op[2] * 11390912356276053911UL) + ((uint64_t)op[3] * 7999781594882481118UL) + ((uint64_t)op[4] * 16533595601628990356UL) + ((uint64_t)op[5] * 1892271431477559865UL) + ((uint64_t)op[6] * 11665282908257879664UL) + ((uint64_t)op[7] * 9581095481545774617UL) + ((uint64_t)op[8] * 4529343932609515736UL) + ((uint64_t)op[9] * 12450654550130001431UL) + ((uint64_t)op[10] * 8849542654538680931UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 8849542654538680931UL) + ((uint64_t)op[1] * 2206005580371624619UL) + ((((uint64_t)op[2] * 10623114139826810515UL) + ((uint64_t)op[3] * 11390912356276053911UL) + ((uint64_t)op[4] * 7999781594882481118UL) + ((uint64_t)op[5] * 16533595601628990356UL) + ((uint64_t)op[6] * 1892271431477559865UL) + ((uint64_t)op[7] * 11665282908257879664UL) + ((uint64_t)op[8] * 9581095481545774617UL) + ((uint64_t)op[9] * 4529343932609515736UL) + ((uint64_t)op[10] * 12450654550130001431UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 12450654550130001431UL) + ((uint64_t)op[1] * 8849542654538680931UL) + ((uint64_t)op[2] * 2206005580371624619UL) + ((((uint64_t)op[3] * 10623114139826810515UL) + ((uint64_t)op[4] * 11390912356276053911UL) + ((uint64_t)op[5] * 7999781594882481118UL) + ((uint64_t)op[6] * 16533595601628990356UL) + ((uint64_t)op[7] * 1892271431477559865UL) + ((uint64_t)op[8] * 11665282908257879664UL) + ((uint64_t)op[9] * 9581095481545774617UL) + ((uint64_t)op[10] * 4529343932609515736UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 4529343932609515736UL) + ((uint64_t)op[1] * 12450654550130001431UL) + ((uint64_t)op[2] * 8849542654538680931UL) + ((uint64_t)op[3] * 2206005580371624619UL) + ((((uint64_t)op[4] * 10623114139826810515UL) + ((uint64_t)op[5] * 11390912356276053911UL) + ((uint64_t)op[6] * 7999781594882481118UL) + ((uint64_t)op[7] * 16533595601628990356UL) + ((uint64_t)op[8] * 1892271431477559865UL) + ((uint64_t)op[9] * 11665282908257879664UL) + ((uint64_t)op[10] * 9581095481545774617UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 9581095481545774617UL) + ((uint64_t)op[1] * 4529343932609515736UL) + ((uint64_t)op[2] * 12450654550130001431UL) + ((uint64_t)op[3] * 8849542654538680931UL) + ((uint64_t)op[4] * 2206005580371624619UL) + ((((uint64_t)op[5] * 10623114139826810515UL) + ((uint64_t)op[6] * 11390912356276053911UL) + ((uint64_t)op[7] * 7999781594882481118UL) + ((uint64_t)op[8] * 16533595601628990356UL) + ((uint64_t)op[9] * 1892271431477559865UL) + ((uint64_t)op[10] * 11665282908257879664UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 11665282908257879664UL) + ((uint64_t)op[1] * 9581095481545774617UL) + ((uint64_t)op[2] * 4529343932609515736UL) + ((uint64_t)op[3] * 12450654550130001431UL) + ((uint64_t)op[4] * 8849542654538680931UL) + ((uint64_t)op[5] * 2206005580371624619UL) + ((((uint64_t)op[6] * 10623114139826810515UL) + ((uint64_t)op[7] * 11390912356276053911UL) + ((uint64_t)op[8] * 7999781594882481118UL) + ((uint64_t)op[9] * 16533595601628990356UL) + ((uint64_t)op[10] * 1892271431477559865UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 1892271431477559865UL) + ((uint64_t)op[1] * 11665282908257879664UL) + ((uint64_t)op[2] * 9581095481545774617UL) + ((uint64_t)op[3] * 4529343932609515736UL) + ((uint64_t)op[4] * 12450654550130001431UL) + ((uint64_t)op[5] * 8849542654538680931UL) + ((uint64_t)op[6] * 2206005580371624619UL) + ((((uint64_t)op[7] * 10623114139826810515UL) + ((uint64_t)op[8] * 11390912356276053911UL) + ((uint64_t)op[9] * 7999781594882481118UL) + ((uint64_t)op[10] * 16533595601628990356UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 16533595601628990356UL) + ((uint64_t)op[1] * 1892271431477559865UL) + ((uint64_t)op[2] * 11665282908257879664UL) + ((uint64_t)op[3] * 9581095481545774617UL) + ((uint64_t)op[4] * 4529343932609515736UL) + ((uint64_t)op[5] * 12450654550130001431UL) + ((uint64_t)op[6] * 8849542654538680931UL) + ((uint64_t)op[7] * 2206005580371624619UL) + ((((uint64_t)op[8] * 10623114139826810515UL) + ((uint64_t)op[9] * 11390912356276053911UL) + ((uint64_t)op[10] * 7999781594882481118UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 7999781594882481118UL) + ((uint64_t)op[1] * 16533595601628990356UL) + ((uint64_t)op[2] * 1892271431477559865UL) + ((uint64_t)op[3] * 11665282908257879664UL) + ((uint64_t)op[4] * 9581095481545774617UL) + ((uint64_t)op[5] * 4529343932609515736UL) + ((uint64_t)op[6] * 12450654550130001431UL) + ((uint64_t)op[7] * 8849542654538680931UL) + ((uint64_t)op[8] * 2206005580371624619UL) + ((((uint64_t)op[9] * 10623114139826810515UL) + ((uint64_t)op[10] * 11390912356276053911UL)) * 18446744073709551614);
	tmp_q[9] = ((uint64_t)op[0] * 11390912356276053911UL) + ((uint64_t)op[1] * 7999781594882481118UL) + ((uint64_t)op[2] * 16533595601628990356UL) + ((uint64_t)op[3] * 1892271431477559865UL) + ((uint64_t)op[4] * 11665282908257879664UL) + ((uint64_t)op[5] * 9581095481545774617UL) + ((uint64_t)op[6] * 4529343932609515736UL) + ((uint64_t)op[7] * 12450654550130001431UL) + ((uint64_t)op[8] * 8849542654538680931UL) + ((uint64_t)op[9] * 2206005580371624619UL) + ((uint64_t)op[10] * 15647259867765482202UL);
	tmp_q[10] = ((uint64_t)op[0] * 10623114139826810515UL) + ((uint64_t)op[1] * 11390912356276053911UL) + ((uint64_t)op[2] * 7999781594882481118UL) + ((uint64_t)op[3] * 16533595601628990356UL) + ((uint64_t)op[4] * 1892271431477559865UL) + ((uint64_t)op[5] * 11665282908257879664UL) + ((uint64_t)op[6] * 9581095481545774617UL) + ((uint64_t)op[7] * 4529343932609515736UL) + ((uint64_t)op[8] * 12450654550130001431UL) + ((uint64_t)op[9] * 8849542654538680931UL) + ((uint64_t)op[10] * 2206005580371624619UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 94473840284775L) - ((-((int128)tmp_q[1] * 50328613663276L) - ((int128)tmp_q[2] * 46124838709299L) - ((int128)tmp_q[3] * 98671693676426L) - ((int128)tmp_q[4] * 62911759865559L) - ((int128)tmp_q[5] * 7227152500441L) + ((int128)tmp_q[6] * 38743706706200L) - ((int128)tmp_q[7] * 99619338434700L) - ((int128)tmp_q[8] * 29311114672031L) + ((int128)tmp_q[9] * 68558203687790L) - ((int128)tmp_q[10] * 60349790004083L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 60349790004083L) + ((int128)tmp_q[1] * 94473840284775L) - ((-((int128)tmp_q[2] * 50328613663276L) - ((int128)tmp_q[3] * 46124838709299L) - ((int128)tmp_q[4] * 98671693676426L) - ((int128)tmp_q[5] * 62911759865559L) - ((int128)tmp_q[6] * 7227152500441L) + ((int128)tmp_q[7] * 38743706706200L) - ((int128)tmp_q[8] * 99619338434700L) - ((int128)tmp_q[9] * 29311114672031L) + ((int128)tmp_q[10] * 68558203687790L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 68558203687790L) - ((int128)tmp_q[1] * 60349790004083L) + ((int128)tmp_q[2] * 94473840284775L) - ((-((int128)tmp_q[3] * 50328613663276L) - ((int128)tmp_q[4] * 46124838709299L) - ((int128)tmp_q[5] * 98671693676426L) - ((int128)tmp_q[6] * 62911759865559L) - ((int128)tmp_q[7] * 7227152500441L) + ((int128)tmp_q[8] * 38743706706200L) - ((int128)tmp_q[9] * 99619338434700L) - ((int128)tmp_q[10] * 29311114672031L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 29311114672031L) + ((int128)tmp_q[1] * 68558203687790L) - ((int128)tmp_q[2] * 60349790004083L) + ((int128)tmp_q[3] * 94473840284775L) - ((-((int128)tmp_q[4] * 50328613663276L) - ((int128)tmp_q[5] * 46124838709299L) - ((int128)tmp_q[6] * 98671693676426L) - ((int128)tmp_q[7] * 62911759865559L) - ((int128)tmp_q[8] * 7227152500441L) + ((int128)tmp_q[9] * 38743706706200L) - ((int128)tmp_q[10] * 99619338434700L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 99619338434700L) - ((int128)tmp_q[1] * 29311114672031L) + ((int128)tmp_q[2] * 68558203687790L) - ((int128)tmp_q[3] * 60349790004083L) + ((int128)tmp_q[4] * 94473840284775L) - ((-((int128)tmp_q[5] * 50328613663276L) - ((int128)tmp_q[6] * 46124838709299L) - ((int128)tmp_q[7] * 98671693676426L) - ((int128)tmp_q[8] * 62911759865559L) - ((int128)tmp_q[9] * 7227152500441L) + ((int128)tmp_q[10] * 38743706706200L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 38743706706200L) - ((int128)tmp_q[1] * 99619338434700L) - ((int128)tmp_q[2] * 29311114672031L) + ((int128)tmp_q[3] * 68558203687790L) - ((int128)tmp_q[4] * 60349790004083L) + ((int128)tmp_q[5] * 94473840284775L) - ((-((int128)tmp_q[6] * 50328613663276L) - ((int128)tmp_q[7] * 46124838709299L) - ((int128)tmp_q[8] * 98671693676426L) - ((int128)tmp_q[9] * 62911759865559L) - ((int128)tmp_q[10] * 7227152500441L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 7227152500441L) + ((int128)tmp_q[1] * 38743706706200L) - ((int128)tmp_q[2] * 99619338434700L) - ((int128)tmp_q[3] * 29311114672031L) + ((int128)tmp_q[4] * 68558203687790L) - ((int128)tmp_q[5] * 60349790004083L) + ((int128)tmp_q[6] * 94473840284775L) - ((-((int128)tmp_q[7] * 50328613663276L) - ((int128)tmp_q[8] * 46124838709299L) - ((int128)tmp_q[9] * 98671693676426L) - ((int128)tmp_q[10] * 62911759865559L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 62911759865559L) - ((int128)tmp_q[1] * 7227152500441L) + ((int128)tmp_q[2] * 38743706706200L) - ((int128)tmp_q[3] * 99619338434700L) - ((int128)tmp_q[4] * 29311114672031L) + ((int128)tmp_q[5] * 68558203687790L) - ((int128)tmp_q[6] * 60349790004083L) + ((int128)tmp_q[7] * 94473840284775L) - ((-((int128)tmp_q[8] * 50328613663276L) - ((int128)tmp_q[9] * 46124838709299L) - ((int128)tmp_q[10] * 98671693676426L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 98671693676426L) - ((int128)tmp_q[1] * 62911759865559L) - ((int128)tmp_q[2] * 7227152500441L) + ((int128)tmp_q[3] * 38743706706200L) - ((int128)tmp_q[4] * 99619338434700L) - ((int128)tmp_q[5] * 29311114672031L) + ((int128)tmp_q[6] * 68558203687790L) - ((int128)tmp_q[7] * 60349790004083L) + ((int128)tmp_q[8] * 94473840284775L) - ((-((int128)tmp_q[9] * 50328613663276L) - ((int128)tmp_q[10] * 46124838709299L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 46124838709299L) - ((int128)tmp_q[1] * 98671693676426L) - ((int128)tmp_q[2] * 62911759865559L) - ((int128)tmp_q[3] * 7227152500441L) + ((int128)tmp_q[4] * 38743706706200L) - ((int128)tmp_q[5] * 99619338434700L) - ((int128)tmp_q[6] * 29311114672031L) + ((int128)tmp_q[7] * 68558203687790L) - ((int128)tmp_q[8] * 60349790004083L) + ((int128)tmp_q[9] * 94473840284775L) + ((int128)tmp_q[10] * 100657227326552L);
	tmp_zero[10] = -((int128)tmp_q[0] * 50328613663276L) - ((int128)tmp_q[1] * 46124838709299L) - ((int128)tmp_q[2] * 98671693676426L) - ((int128)tmp_q[3] * 62911759865559L) - ((int128)tmp_q[4] * 7227152500441L) + ((int128)tmp_q[5] * 38743706706200L) - ((int128)tmp_q[6] * 99619338434700L) - ((int128)tmp_q[7] * 29311114672031L) + ((int128)tmp_q[8] * 68558203687790L) - ((int128)tmp_q[9] * 60349790004083L) + ((int128)tmp_q[10] * 94473840284775L);

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
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
}

