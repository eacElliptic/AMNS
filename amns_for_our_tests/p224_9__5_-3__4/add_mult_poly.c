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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8336592734456880092UL) + ((((uint64_t)op[1] * 5181162617861188834UL) + ((uint64_t)op[2] * 8500272620505367746UL) + ((uint64_t)op[3] * 16269859939925057841UL) + ((uint64_t)op[4] * 16445353490661335824UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16445353490661335824UL) + ((uint64_t)op[1] * 8336592734456880092UL) + ((((uint64_t)op[2] * 5181162617861188834UL) + ((uint64_t)op[3] * 8500272620505367746UL) + ((uint64_t)op[4] * 16269859939925057841UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 16269859939925057841UL) + ((uint64_t)op[1] * 16445353490661335824UL) + ((uint64_t)op[2] * 8336592734456880092UL) + ((((uint64_t)op[3] * 5181162617861188834UL) + ((uint64_t)op[4] * 8500272620505367746UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 8500272620505367746UL) + ((uint64_t)op[1] * 16269859939925057841UL) + ((uint64_t)op[2] * 16445353490661335824UL) + ((uint64_t)op[3] * 8336592734456880092UL) + ((uint64_t)op[4] * 2903256220125985114UL);
	tmp_q[4] = ((uint64_t)op[0] * 5181162617861188834UL) + ((uint64_t)op[1] * 8500272620505367746UL) + ((uint64_t)op[2] * 16269859939925057841UL) + ((uint64_t)op[3] * 16445353490661335824UL) + ((uint64_t)op[4] * 8336592734456880092UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 9188833394674L) - ((-((int128)tmp_q[1] * 19453427850878L) - ((int128)tmp_q[2] * 21806029569133L) - ((int128)tmp_q[3] * 2980959200188L) - ((int128)tmp_q[4] * 3970953282852L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 3970953282852L) - ((int128)tmp_q[1] * 9188833394674L) - ((-((int128)tmp_q[2] * 19453427850878L) - ((int128)tmp_q[3] * 21806029569133L) - ((int128)tmp_q[4] * 2980959200188L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 2980959200188L) - ((int128)tmp_q[1] * 3970953282852L) - ((int128)tmp_q[2] * 9188833394674L) - ((-((int128)tmp_q[3] * 19453427850878L) - ((int128)tmp_q[4] * 21806029569133L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 21806029569133L) - ((int128)tmp_q[1] * 2980959200188L) - ((int128)tmp_q[2] * 3970953282852L) - ((int128)tmp_q[3] * 9188833394674L) + ((int128)tmp_q[4] * 58360283552634L);
	tmp_zero[4] = -((int128)tmp_q[0] * 19453427850878L) - ((int128)tmp_q[1] * 21806029569133L) - ((int128)tmp_q[2] * 2980959200188L) - ((int128)tmp_q[3] * 3970953282852L) - ((int128)tmp_q[4] * 9188833394674L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

