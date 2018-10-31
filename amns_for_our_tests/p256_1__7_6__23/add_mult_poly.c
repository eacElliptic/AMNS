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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 72580153127500891UL) + ((((uint64_t)op[1] * 8468214800066082966UL) + ((uint64_t)op[2] * 55936323793264674UL) + ((uint64_t)op[3] * 12115706106394598320UL) + ((uint64_t)op[4] * 2899075337491808615UL) + ((uint64_t)op[5] * 13275765700778003002UL) + ((uint64_t)op[6] * 7212814401813249662UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 7212814401813249662UL) + ((uint64_t)op[1] * 72580153127500891UL) + ((((uint64_t)op[2] * 8468214800066082966UL) + ((uint64_t)op[3] * 55936323793264674UL) + ((uint64_t)op[4] * 12115706106394598320UL) + ((uint64_t)op[5] * 2899075337491808615UL) + ((uint64_t)op[6] * 13275765700778003002UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 13275765700778003002UL) + ((uint64_t)op[1] * 7212814401813249662UL) + ((uint64_t)op[2] * 72580153127500891UL) + ((((uint64_t)op[3] * 8468214800066082966UL) + ((uint64_t)op[4] * 55936323793264674UL) + ((uint64_t)op[5] * 12115706106394598320UL) + ((uint64_t)op[6] * 2899075337491808615UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 2899075337491808615UL) + ((uint64_t)op[1] * 13275765700778003002UL) + ((uint64_t)op[2] * 7212814401813249662UL) + ((uint64_t)op[3] * 72580153127500891UL) + ((((uint64_t)op[4] * 8468214800066082966UL) + ((uint64_t)op[5] * 55936323793264674UL) + ((uint64_t)op[6] * 12115706106394598320UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 12115706106394598320UL) + ((uint64_t)op[1] * 2899075337491808615UL) + ((uint64_t)op[2] * 13275765700778003002UL) + ((uint64_t)op[3] * 7212814401813249662UL) + ((uint64_t)op[4] * 72580153127500891UL) + ((((uint64_t)op[5] * 8468214800066082966UL) + ((uint64_t)op[6] * 55936323793264674UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 55936323793264674UL) + ((uint64_t)op[1] * 12115706106394598320UL) + ((uint64_t)op[2] * 2899075337491808615UL) + ((uint64_t)op[3] * 13275765700778003002UL) + ((uint64_t)op[4] * 7212814401813249662UL) + ((uint64_t)op[5] * 72580153127500891UL) + ((uint64_t)op[6] * 13915800652977394564UL);
	tmp_q[6] = ((uint64_t)op[0] * 8468214800066082966UL) + ((uint64_t)op[1] * 55936323793264674UL) + ((uint64_t)op[2] * 12115706106394598320UL) + ((uint64_t)op[3] * 2899075337491808615UL) + ((uint64_t)op[4] * 13275765700778003002UL) + ((uint64_t)op[5] * 7212814401813249662UL) + ((uint64_t)op[6] * 72580153127500891UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 78996798769L) + ((-((int128)tmp_q[1] * 10446709073L) + ((int128)tmp_q[2] * 16874725412L) + ((int128)tmp_q[3] * 58124376880L) - ((int128)tmp_q[4] * 54054467737L) - ((int128)tmp_q[5] * 44111315296L) + ((int128)tmp_q[6] * 8100322598L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 8100322598L) + ((int128)tmp_q[1] * 78996798769L) + ((-((int128)tmp_q[2] * 10446709073L) + ((int128)tmp_q[3] * 16874725412L) + ((int128)tmp_q[4] * 58124376880L) - ((int128)tmp_q[5] * 54054467737L) - ((int128)tmp_q[6] * 44111315296L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 44111315296L) + ((int128)tmp_q[1] * 8100322598L) + ((int128)tmp_q[2] * 78996798769L) + ((-((int128)tmp_q[3] * 10446709073L) + ((int128)tmp_q[4] * 16874725412L) + ((int128)tmp_q[5] * 58124376880L) - ((int128)tmp_q[6] * 54054467737L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 54054467737L) - ((int128)tmp_q[1] * 44111315296L) + ((int128)tmp_q[2] * 8100322598L) + ((int128)tmp_q[3] * 78996798769L) + ((-((int128)tmp_q[4] * 10446709073L) + ((int128)tmp_q[5] * 16874725412L) + ((int128)tmp_q[6] * 58124376880L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 58124376880L) - ((int128)tmp_q[1] * 54054467737L) - ((int128)tmp_q[2] * 44111315296L) + ((int128)tmp_q[3] * 8100322598L) + ((int128)tmp_q[4] * 78996798769L) + ((-((int128)tmp_q[5] * 10446709073L) + ((int128)tmp_q[6] * 16874725412L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 16874725412L) + ((int128)tmp_q[1] * 58124376880L) - ((int128)tmp_q[2] * 54054467737L) - ((int128)tmp_q[3] * 44111315296L) + ((int128)tmp_q[4] * 8100322598L) + ((int128)tmp_q[5] * 78996798769L) - ((int128)tmp_q[6] * 62680254438L);
	tmp_zero[6] = -((int128)tmp_q[0] * 10446709073L) + ((int128)tmp_q[1] * 16874725412L) + ((int128)tmp_q[2] * 58124376880L) - ((int128)tmp_q[3] * 54054467737L) - ((int128)tmp_q[4] * 44111315296L) + ((int128)tmp_q[5] * 8100322598L) + ((int128)tmp_q[6] * 78996798769L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

