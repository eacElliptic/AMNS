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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - ((int128)pa[6] * pb[6]);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - ((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - ((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - ((int128)pa[6] * pa[6]);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13287208854219371361UL) + ((((uint64_t)op[1] * 14701228477354331046UL) + ((uint64_t)op[2] * 4160081460244702887UL) + ((uint64_t)op[3] * 15774100135670700277UL) + ((uint64_t)op[4] * 1637499953191379549UL) + ((uint64_t)op[5] * 13620960953298965930UL) + ((uint64_t)op[6] * 2394599477215148535UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 2394599477215148535UL) + ((uint64_t)op[1] * 13287208854219371361UL) + ((((uint64_t)op[2] * 14701228477354331046UL) + ((uint64_t)op[3] * 4160081460244702887UL) + ((uint64_t)op[4] * 15774100135670700277UL) + ((uint64_t)op[5] * 1637499953191379549UL) + ((uint64_t)op[6] * 13620960953298965930UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 13620960953298965930UL) + ((uint64_t)op[1] * 2394599477215148535UL) + ((uint64_t)op[2] * 13287208854219371361UL) + ((((uint64_t)op[3] * 14701228477354331046UL) + ((uint64_t)op[4] * 4160081460244702887UL) + ((uint64_t)op[5] * 15774100135670700277UL) + ((uint64_t)op[6] * 1637499953191379549UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 1637499953191379549UL) + ((uint64_t)op[1] * 13620960953298965930UL) + ((uint64_t)op[2] * 2394599477215148535UL) + ((uint64_t)op[3] * 13287208854219371361UL) + ((((uint64_t)op[4] * 14701228477354331046UL) + ((uint64_t)op[5] * 4160081460244702887UL) + ((uint64_t)op[6] * 15774100135670700277UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 15774100135670700277UL) + ((uint64_t)op[1] * 1637499953191379549UL) + ((uint64_t)op[2] * 13620960953298965930UL) + ((uint64_t)op[3] * 2394599477215148535UL) + ((uint64_t)op[4] * 13287208854219371361UL) + ((((uint64_t)op[5] * 14701228477354331046UL) + ((uint64_t)op[6] * 4160081460244702887UL)) * 18446744073709551615);
	tmp_q[5] = ((uint64_t)op[0] * 4160081460244702887UL) + ((uint64_t)op[1] * 15774100135670700277UL) + ((uint64_t)op[2] * 1637499953191379549UL) + ((uint64_t)op[3] * 13620960953298965930UL) + ((uint64_t)op[4] * 2394599477215148535UL) + ((uint64_t)op[5] * 13287208854219371361UL) + ((uint64_t)op[6] * 3745515596355220570UL);
	tmp_q[6] = ((uint64_t)op[0] * 14701228477354331046UL) + ((uint64_t)op[1] * 4160081460244702887UL) + ((uint64_t)op[2] * 15774100135670700277UL) + ((uint64_t)op[3] * 1637499953191379549UL) + ((uint64_t)op[4] * 13620960953298965930UL) + ((uint64_t)op[5] * 2394599477215148535UL) + ((uint64_t)op[6] * 13287208854219371361UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1245108078502L) - (((int128)tmp_q[1] * 63643841648L) - ((int128)tmp_q[2] * 591425512536L) + ((int128)tmp_q[3] * 3500575134611L) + ((int128)tmp_q[4] * 3348970396451L) - ((int128)tmp_q[5] * 3292244314285L) - ((int128)tmp_q[6] * 3730678300440L));
	tmp_zero[1] = -((int128)tmp_q[0] * 3730678300440L) - ((int128)tmp_q[1] * 1245108078502L) - (((int128)tmp_q[2] * 63643841648L) - ((int128)tmp_q[3] * 591425512536L) + ((int128)tmp_q[4] * 3500575134611L) + ((int128)tmp_q[5] * 3348970396451L) - ((int128)tmp_q[6] * 3292244314285L));
	tmp_zero[2] = -((int128)tmp_q[0] * 3292244314285L) - ((int128)tmp_q[1] * 3730678300440L) - ((int128)tmp_q[2] * 1245108078502L) - (((int128)tmp_q[3] * 63643841648L) - ((int128)tmp_q[4] * 591425512536L) + ((int128)tmp_q[5] * 3500575134611L) + ((int128)tmp_q[6] * 3348970396451L));
	tmp_zero[3] = ((int128)tmp_q[0] * 3348970396451L) - ((int128)tmp_q[1] * 3292244314285L) - ((int128)tmp_q[2] * 3730678300440L) - ((int128)tmp_q[3] * 1245108078502L) - (((int128)tmp_q[4] * 63643841648L) - ((int128)tmp_q[5] * 591425512536L) + ((int128)tmp_q[6] * 3500575134611L));
	tmp_zero[4] = ((int128)tmp_q[0] * 3500575134611L) + ((int128)tmp_q[1] * 3348970396451L) - ((int128)tmp_q[2] * 3292244314285L) - ((int128)tmp_q[3] * 3730678300440L) - ((int128)tmp_q[4] * 1245108078502L) - (((int128)tmp_q[5] * 63643841648L) - ((int128)tmp_q[6] * 591425512536L));
	tmp_zero[5] = -((int128)tmp_q[0] * 591425512536L) + ((int128)tmp_q[1] * 3500575134611L) + ((int128)tmp_q[2] * 3348970396451L) - ((int128)tmp_q[3] * 3292244314285L) - ((int128)tmp_q[4] * 3730678300440L) - ((int128)tmp_q[5] * 1245108078502L) - ((int128)tmp_q[6] * 63643841648L);
	tmp_zero[6] = ((int128)tmp_q[0] * 63643841648L) - ((int128)tmp_q[1] * 591425512536L) + ((int128)tmp_q[2] * 3500575134611L) + ((int128)tmp_q[3] * 3348970396451L) - ((int128)tmp_q[4] * 3292244314285L) - ((int128)tmp_q[5] * 3730678300440L) - ((int128)tmp_q[6] * 1245108078502L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

