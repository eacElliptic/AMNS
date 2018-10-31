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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13677045251406390527UL) + ((((uint64_t)op[1] * 10811934133603809114UL) + ((uint64_t)op[2] * 751928066400947453UL) + ((uint64_t)op[3] * 7065825107985564858UL) + ((uint64_t)op[4] * 12002460066399956143UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 12002460066399956143UL) + ((uint64_t)op[1] * 13677045251406390527UL) + ((((uint64_t)op[2] * 10811934133603809114UL) + ((uint64_t)op[3] * 751928066400947453UL) + ((uint64_t)op[4] * 7065825107985564858UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 7065825107985564858UL) + ((uint64_t)op[1] * 12002460066399956143UL) + ((uint64_t)op[2] * 13677045251406390527UL) + ((((uint64_t)op[3] * 10811934133603809114UL) + ((uint64_t)op[4] * 751928066400947453UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 751928066400947453UL) + ((uint64_t)op[1] * 7065825107985564858UL) + ((uint64_t)op[2] * 12002460066399956143UL) + ((uint64_t)op[3] * 13677045251406390527UL) + ((uint64_t)op[4] * 16550181433321094282UL);
	tmp_q[4] = ((uint64_t)op[0] * 10811934133603809114UL) + ((uint64_t)op[1] * 751928066400947453UL) + ((uint64_t)op[2] * 7065825107985564858UL) + ((uint64_t)op[3] * 12002460066399956143UL) + ((uint64_t)op[4] * 13677045251406390527UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 9585966997738L) - ((-((int128)tmp_q[1] * 1683582114914L) - ((int128)tmp_q[2] * 8602267392223L) - ((int128)tmp_q[3] * 4040377121949L) - ((int128)tmp_q[4] * 18916163212919L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 18916163212919L) - ((int128)tmp_q[1] * 9585966997738L) - ((-((int128)tmp_q[2] * 1683582114914L) - ((int128)tmp_q[3] * 8602267392223L) - ((int128)tmp_q[4] * 4040377121949L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 4040377121949L) - ((int128)tmp_q[1] * 18916163212919L) - ((int128)tmp_q[2] * 9585966997738L) - ((-((int128)tmp_q[3] * 1683582114914L) - ((int128)tmp_q[4] * 8602267392223L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 8602267392223L) - ((int128)tmp_q[1] * 4040377121949L) - ((int128)tmp_q[2] * 18916163212919L) - ((int128)tmp_q[3] * 9585966997738L) + ((int128)tmp_q[4] * 11785074804398L);
	tmp_zero[4] = -((int128)tmp_q[0] * 1683582114914L) - ((int128)tmp_q[1] * 8602267392223L) - ((int128)tmp_q[2] * 4040377121949L) - ((int128)tmp_q[3] * 18916163212919L) - ((int128)tmp_q[4] * 9585966997738L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

