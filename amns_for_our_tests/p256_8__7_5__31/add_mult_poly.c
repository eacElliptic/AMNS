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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6936742821827891967UL) + ((((uint64_t)op[1] * 3522174349041649439UL) + ((uint64_t)op[2] * 9059561562959172719UL) + ((uint64_t)op[3] * 14289489603886369903UL) + ((uint64_t)op[4] * 17297078295116976569UL) + ((uint64_t)op[5] * 6710362729397594830UL) + ((uint64_t)op[6] * 6878654827637845538UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 6878654827637845538UL) + ((uint64_t)op[1] * 6936742821827891967UL) + ((((uint64_t)op[2] * 3522174349041649439UL) + ((uint64_t)op[3] * 9059561562959172719UL) + ((uint64_t)op[4] * 14289489603886369903UL) + ((uint64_t)op[5] * 17297078295116976569UL) + ((uint64_t)op[6] * 6710362729397594830UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 6710362729397594830UL) + ((uint64_t)op[1] * 6878654827637845538UL) + ((uint64_t)op[2] * 6936742821827891967UL) + ((((uint64_t)op[3] * 3522174349041649439UL) + ((uint64_t)op[4] * 9059561562959172719UL) + ((uint64_t)op[5] * 14289489603886369903UL) + ((uint64_t)op[6] * 17297078295116976569UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 17297078295116976569UL) + ((uint64_t)op[1] * 6710362729397594830UL) + ((uint64_t)op[2] * 6878654827637845538UL) + ((uint64_t)op[3] * 6936742821827891967UL) + ((((uint64_t)op[4] * 3522174349041649439UL) + ((uint64_t)op[5] * 9059561562959172719UL) + ((uint64_t)op[6] * 14289489603886369903UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 14289489603886369903UL) + ((uint64_t)op[1] * 17297078295116976569UL) + ((uint64_t)op[2] * 6710362729397594830UL) + ((uint64_t)op[3] * 6878654827637845538UL) + ((uint64_t)op[4] * 6936742821827891967UL) + ((((uint64_t)op[5] * 3522174349041649439UL) + ((uint64_t)op[6] * 9059561562959172719UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 9059561562959172719UL) + ((uint64_t)op[1] * 14289489603886369903UL) + ((uint64_t)op[2] * 17297078295116976569UL) + ((uint64_t)op[3] * 6710362729397594830UL) + ((uint64_t)op[4] * 6878654827637845538UL) + ((uint64_t)op[5] * 6936742821827891967UL) + ((uint64_t)op[6] * 17610871745208247195UL);
	tmp_q[6] = ((uint64_t)op[0] * 3522174349041649439UL) + ((uint64_t)op[1] * 9059561562959172719UL) + ((uint64_t)op[2] * 14289489603886369903UL) + ((uint64_t)op[3] * 17297078295116976569UL) + ((uint64_t)op[4] * 6710362729397594830UL) + ((uint64_t)op[5] * 6878654827637845538UL) + ((uint64_t)op[6] * 6936742821827891967UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3022471213L) + ((-((int128)tmp_q[1] * 39727018660L) + ((int128)tmp_q[2] * 6086297196L) + ((int128)tmp_q[3] * 54888947469L) - ((int128)tmp_q[4] * 51323590998L) - ((int128)tmp_q[5] * 61451126101L) - ((int128)tmp_q[6] * 40070525276L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 40070525276L) + ((int128)tmp_q[1] * 3022471213L) + ((-((int128)tmp_q[2] * 39727018660L) + ((int128)tmp_q[3] * 6086297196L) + ((int128)tmp_q[4] * 54888947469L) - ((int128)tmp_q[5] * 51323590998L) - ((int128)tmp_q[6] * 61451126101L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 61451126101L) - ((int128)tmp_q[1] * 40070525276L) + ((int128)tmp_q[2] * 3022471213L) + ((-((int128)tmp_q[3] * 39727018660L) + ((int128)tmp_q[4] * 6086297196L) + ((int128)tmp_q[5] * 54888947469L) - ((int128)tmp_q[6] * 51323590998L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 51323590998L) - ((int128)tmp_q[1] * 61451126101L) - ((int128)tmp_q[2] * 40070525276L) + ((int128)tmp_q[3] * 3022471213L) + ((-((int128)tmp_q[4] * 39727018660L) + ((int128)tmp_q[5] * 6086297196L) + ((int128)tmp_q[6] * 54888947469L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 54888947469L) - ((int128)tmp_q[1] * 51323590998L) - ((int128)tmp_q[2] * 61451126101L) - ((int128)tmp_q[3] * 40070525276L) + ((int128)tmp_q[4] * 3022471213L) + ((-((int128)tmp_q[5] * 39727018660L) + ((int128)tmp_q[6] * 6086297196L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 6086297196L) + ((int128)tmp_q[1] * 54888947469L) - ((int128)tmp_q[2] * 51323590998L) - ((int128)tmp_q[3] * 61451126101L) - ((int128)tmp_q[4] * 40070525276L) + ((int128)tmp_q[5] * 3022471213L) - ((int128)tmp_q[6] * 198635093300L);
	tmp_zero[6] = -((int128)tmp_q[0] * 39727018660L) + ((int128)tmp_q[1] * 6086297196L) + ((int128)tmp_q[2] * 54888947469L) - ((int128)tmp_q[3] * 51323590998L) - ((int128)tmp_q[4] * 61451126101L) - ((int128)tmp_q[5] * 40070525276L) + ((int128)tmp_q[6] * 3022471213L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

