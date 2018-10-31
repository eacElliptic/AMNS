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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2325747854744013109UL) + ((((uint64_t)op[1] * 886056309200235167UL) + ((uint64_t)op[2] * 13286684366955195329UL) + ((uint64_t)op[3] * 732015700635118907UL) + ((uint64_t)op[4] * 16372329120823111132UL) + ((uint64_t)op[5] * 13600466792135336697UL) + ((uint64_t)op[6] * 14335668095179292004UL) + ((uint64_t)op[7] * 17621548976691877982UL) + ((uint64_t)op[8] * 4911398843066453487UL) + ((uint64_t)op[9] * 5864199075027693142UL) + ((uint64_t)op[10] * 8814664365608660142UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 8814664365608660142UL) + ((uint64_t)op[1] * 2325747854744013109UL) + ((((uint64_t)op[2] * 886056309200235167UL) + ((uint64_t)op[3] * 13286684366955195329UL) + ((uint64_t)op[4] * 732015700635118907UL) + ((uint64_t)op[5] * 16372329120823111132UL) + ((uint64_t)op[6] * 13600466792135336697UL) + ((uint64_t)op[7] * 14335668095179292004UL) + ((uint64_t)op[8] * 17621548976691877982UL) + ((uint64_t)op[9] * 4911398843066453487UL) + ((uint64_t)op[10] * 5864199075027693142UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 5864199075027693142UL) + ((uint64_t)op[1] * 8814664365608660142UL) + ((uint64_t)op[2] * 2325747854744013109UL) + ((((uint64_t)op[3] * 886056309200235167UL) + ((uint64_t)op[4] * 13286684366955195329UL) + ((uint64_t)op[5] * 732015700635118907UL) + ((uint64_t)op[6] * 16372329120823111132UL) + ((uint64_t)op[7] * 13600466792135336697UL) + ((uint64_t)op[8] * 14335668095179292004UL) + ((uint64_t)op[9] * 17621548976691877982UL) + ((uint64_t)op[10] * 4911398843066453487UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 4911398843066453487UL) + ((uint64_t)op[1] * 5864199075027693142UL) + ((uint64_t)op[2] * 8814664365608660142UL) + ((uint64_t)op[3] * 2325747854744013109UL) + ((((uint64_t)op[4] * 886056309200235167UL) + ((uint64_t)op[5] * 13286684366955195329UL) + ((uint64_t)op[6] * 732015700635118907UL) + ((uint64_t)op[7] * 16372329120823111132UL) + ((uint64_t)op[8] * 13600466792135336697UL) + ((uint64_t)op[9] * 14335668095179292004UL) + ((uint64_t)op[10] * 17621548976691877982UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 17621548976691877982UL) + ((uint64_t)op[1] * 4911398843066453487UL) + ((uint64_t)op[2] * 5864199075027693142UL) + ((uint64_t)op[3] * 8814664365608660142UL) + ((uint64_t)op[4] * 2325747854744013109UL) + ((((uint64_t)op[5] * 886056309200235167UL) + ((uint64_t)op[6] * 13286684366955195329UL) + ((uint64_t)op[7] * 732015700635118907UL) + ((uint64_t)op[8] * 16372329120823111132UL) + ((uint64_t)op[9] * 13600466792135336697UL) + ((uint64_t)op[10] * 14335668095179292004UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 14335668095179292004UL) + ((uint64_t)op[1] * 17621548976691877982UL) + ((uint64_t)op[2] * 4911398843066453487UL) + ((uint64_t)op[3] * 5864199075027693142UL) + ((uint64_t)op[4] * 8814664365608660142UL) + ((uint64_t)op[5] * 2325747854744013109UL) + ((((uint64_t)op[6] * 886056309200235167UL) + ((uint64_t)op[7] * 13286684366955195329UL) + ((uint64_t)op[8] * 732015700635118907UL) + ((uint64_t)op[9] * 16372329120823111132UL) + ((uint64_t)op[10] * 13600466792135336697UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 13600466792135336697UL) + ((uint64_t)op[1] * 14335668095179292004UL) + ((uint64_t)op[2] * 17621548976691877982UL) + ((uint64_t)op[3] * 4911398843066453487UL) + ((uint64_t)op[4] * 5864199075027693142UL) + ((uint64_t)op[5] * 8814664365608660142UL) + ((uint64_t)op[6] * 2325747854744013109UL) + ((((uint64_t)op[7] * 886056309200235167UL) + ((uint64_t)op[8] * 13286684366955195329UL) + ((uint64_t)op[9] * 732015700635118907UL) + ((uint64_t)op[10] * 16372329120823111132UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 16372329120823111132UL) + ((uint64_t)op[1] * 13600466792135336697UL) + ((uint64_t)op[2] * 14335668095179292004UL) + ((uint64_t)op[3] * 17621548976691877982UL) + ((uint64_t)op[4] * 4911398843066453487UL) + ((uint64_t)op[5] * 5864199075027693142UL) + ((uint64_t)op[6] * 8814664365608660142UL) + ((uint64_t)op[7] * 2325747854744013109UL) + ((((uint64_t)op[8] * 886056309200235167UL) + ((uint64_t)op[9] * 13286684366955195329UL) + ((uint64_t)op[10] * 732015700635118907UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 732015700635118907UL) + ((uint64_t)op[1] * 16372329120823111132UL) + ((uint64_t)op[2] * 13600466792135336697UL) + ((uint64_t)op[3] * 14335668095179292004UL) + ((uint64_t)op[4] * 17621548976691877982UL) + ((uint64_t)op[5] * 4911398843066453487UL) + ((uint64_t)op[6] * 5864199075027693142UL) + ((uint64_t)op[7] * 8814664365608660142UL) + ((uint64_t)op[8] * 2325747854744013109UL) + ((((uint64_t)op[9] * 886056309200235167UL) + ((uint64_t)op[10] * 13286684366955195329UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 13286684366955195329UL) + ((uint64_t)op[1] * 732015700635118907UL) + ((uint64_t)op[2] * 16372329120823111132UL) + ((uint64_t)op[3] * 13600466792135336697UL) + ((uint64_t)op[4] * 14335668095179292004UL) + ((uint64_t)op[5] * 17621548976691877982UL) + ((uint64_t)op[6] * 4911398843066453487UL) + ((uint64_t)op[7] * 5864199075027693142UL) + ((uint64_t)op[8] * 8814664365608660142UL) + ((uint64_t)op[9] * 2325747854744013109UL) + ((uint64_t)op[10] * 1772112618400470334UL);
	tmp_q[10] = ((uint64_t)op[0] * 886056309200235167UL) + ((uint64_t)op[1] * 13286684366955195329UL) + ((uint64_t)op[2] * 732015700635118907UL) + ((uint64_t)op[3] * 16372329120823111132UL) + ((uint64_t)op[4] * 13600466792135336697UL) + ((uint64_t)op[5] * 14335668095179292004UL) + ((uint64_t)op[6] * 17621548976691877982UL) + ((uint64_t)op[7] * 4911398843066453487UL) + ((uint64_t)op[8] * 5864199075027693142UL) + ((uint64_t)op[9] * 8814664365608660142UL) + ((uint64_t)op[10] * 2325747854744013109UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 45711843815181L) + ((-((int128)tmp_q[1] * 22904268989603L) - ((int128)tmp_q[2] * 9553601434968L) + ((int128)tmp_q[3] * 36420883907103L) + ((int128)tmp_q[4] * 66654373226494L) - ((int128)tmp_q[5] * 67758429304040L) + ((int128)tmp_q[6] * 29294333443620L) + ((int128)tmp_q[7] * 59747763969240L) + ((int128)tmp_q[8] * 8945956931177L) - ((int128)tmp_q[9] * 62314421255818L) + ((int128)tmp_q[10] * 77605883461296L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 77605883461296L) - ((int128)tmp_q[1] * 45711843815181L) + ((-((int128)tmp_q[2] * 22904268989603L) - ((int128)tmp_q[3] * 9553601434968L) + ((int128)tmp_q[4] * 36420883907103L) + ((int128)tmp_q[5] * 66654373226494L) - ((int128)tmp_q[6] * 67758429304040L) + ((int128)tmp_q[7] * 29294333443620L) + ((int128)tmp_q[8] * 59747763969240L) + ((int128)tmp_q[9] * 8945956931177L) - ((int128)tmp_q[10] * 62314421255818L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 62314421255818L) + ((int128)tmp_q[1] * 77605883461296L) - ((int128)tmp_q[2] * 45711843815181L) + ((-((int128)tmp_q[3] * 22904268989603L) - ((int128)tmp_q[4] * 9553601434968L) + ((int128)tmp_q[5] * 36420883907103L) + ((int128)tmp_q[6] * 66654373226494L) - ((int128)tmp_q[7] * 67758429304040L) + ((int128)tmp_q[8] * 29294333443620L) + ((int128)tmp_q[9] * 59747763969240L) + ((int128)tmp_q[10] * 8945956931177L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 8945956931177L) - ((int128)tmp_q[1] * 62314421255818L) + ((int128)tmp_q[2] * 77605883461296L) - ((int128)tmp_q[3] * 45711843815181L) + ((-((int128)tmp_q[4] * 22904268989603L) - ((int128)tmp_q[5] * 9553601434968L) + ((int128)tmp_q[6] * 36420883907103L) + ((int128)tmp_q[7] * 66654373226494L) - ((int128)tmp_q[8] * 67758429304040L) + ((int128)tmp_q[9] * 29294333443620L) + ((int128)tmp_q[10] * 59747763969240L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 59747763969240L) + ((int128)tmp_q[1] * 8945956931177L) - ((int128)tmp_q[2] * 62314421255818L) + ((int128)tmp_q[3] * 77605883461296L) - ((int128)tmp_q[4] * 45711843815181L) + ((-((int128)tmp_q[5] * 22904268989603L) - ((int128)tmp_q[6] * 9553601434968L) + ((int128)tmp_q[7] * 36420883907103L) + ((int128)tmp_q[8] * 66654373226494L) - ((int128)tmp_q[9] * 67758429304040L) + ((int128)tmp_q[10] * 29294333443620L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 29294333443620L) + ((int128)tmp_q[1] * 59747763969240L) + ((int128)tmp_q[2] * 8945956931177L) - ((int128)tmp_q[3] * 62314421255818L) + ((int128)tmp_q[4] * 77605883461296L) - ((int128)tmp_q[5] * 45711843815181L) + ((-((int128)tmp_q[6] * 22904268989603L) - ((int128)tmp_q[7] * 9553601434968L) + ((int128)tmp_q[8] * 36420883907103L) + ((int128)tmp_q[9] * 66654373226494L) - ((int128)tmp_q[10] * 67758429304040L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 67758429304040L) + ((int128)tmp_q[1] * 29294333443620L) + ((int128)tmp_q[2] * 59747763969240L) + ((int128)tmp_q[3] * 8945956931177L) - ((int128)tmp_q[4] * 62314421255818L) + ((int128)tmp_q[5] * 77605883461296L) - ((int128)tmp_q[6] * 45711843815181L) + ((-((int128)tmp_q[7] * 22904268989603L) - ((int128)tmp_q[8] * 9553601434968L) + ((int128)tmp_q[9] * 36420883907103L) + ((int128)tmp_q[10] * 66654373226494L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 66654373226494L) - ((int128)tmp_q[1] * 67758429304040L) + ((int128)tmp_q[2] * 29294333443620L) + ((int128)tmp_q[3] * 59747763969240L) + ((int128)tmp_q[4] * 8945956931177L) - ((int128)tmp_q[5] * 62314421255818L) + ((int128)tmp_q[6] * 77605883461296L) - ((int128)tmp_q[7] * 45711843815181L) + ((-((int128)tmp_q[8] * 22904268989603L) - ((int128)tmp_q[9] * 9553601434968L) + ((int128)tmp_q[10] * 36420883907103L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 36420883907103L) + ((int128)tmp_q[1] * 66654373226494L) - ((int128)tmp_q[2] * 67758429304040L) + ((int128)tmp_q[3] * 29294333443620L) + ((int128)tmp_q[4] * 59747763969240L) + ((int128)tmp_q[5] * 8945956931177L) - ((int128)tmp_q[6] * 62314421255818L) + ((int128)tmp_q[7] * 77605883461296L) - ((int128)tmp_q[8] * 45711843815181L) + ((-((int128)tmp_q[9] * 22904268989603L) - ((int128)tmp_q[10] * 9553601434968L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 9553601434968L) + ((int128)tmp_q[1] * 36420883907103L) + ((int128)tmp_q[2] * 66654373226494L) - ((int128)tmp_q[3] * 67758429304040L) + ((int128)tmp_q[4] * 29294333443620L) + ((int128)tmp_q[5] * 59747763969240L) + ((int128)tmp_q[6] * 8945956931177L) - ((int128)tmp_q[7] * 62314421255818L) + ((int128)tmp_q[8] * 77605883461296L) - ((int128)tmp_q[9] * 45711843815181L) - ((int128)tmp_q[10] * 45808537979206L);
	tmp_zero[10] = -((int128)tmp_q[0] * 22904268989603L) - ((int128)tmp_q[1] * 9553601434968L) + ((int128)tmp_q[2] * 36420883907103L) + ((int128)tmp_q[3] * 66654373226494L) - ((int128)tmp_q[4] * 67758429304040L) + ((int128)tmp_q[5] * 29294333443620L) + ((int128)tmp_q[6] * 59747763969240L) + ((int128)tmp_q[7] * 8945956931177L) - ((int128)tmp_q[8] * 62314421255818L) + ((int128)tmp_q[9] * 77605883461296L) - ((int128)tmp_q[10] * 45711843815181L);

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

