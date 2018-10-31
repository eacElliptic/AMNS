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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15429312352470756473UL) + ((((uint64_t)op[1] * 13571914704634587412UL) + ((uint64_t)op[2] * 7493178428156621871UL) + ((uint64_t)op[3] * 13690863825093704465UL) + ((uint64_t)op[4] * 8005826515000700652UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 8005826515000700652UL) + ((uint64_t)op[1] * 15429312352470756473UL) + ((((uint64_t)op[2] * 13571914704634587412UL) + ((uint64_t)op[3] * 7493178428156621871UL) + ((uint64_t)op[4] * 13690863825093704465UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 13690863825093704465UL) + ((uint64_t)op[1] * 8005826515000700652UL) + ((uint64_t)op[2] * 15429312352470756473UL) + ((((uint64_t)op[3] * 13571914704634587412UL) + ((uint64_t)op[4] * 7493178428156621871UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 7493178428156621871UL) + ((uint64_t)op[1] * 13690863825093704465UL) + ((uint64_t)op[2] * 8005826515000700652UL) + ((uint64_t)op[3] * 15429312352470756473UL) + ((uint64_t)op[4] * 16341597268528941216UL);
	tmp_q[4] = ((uint64_t)op[0] * 13571914704634587412UL) + ((uint64_t)op[1] * 7493178428156621871UL) + ((uint64_t)op[2] * 13690863825093704465UL) + ((uint64_t)op[3] * 8005826515000700652UL) + ((uint64_t)op[4] * 15429312352470756473UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 19535472087335L) + ((-((int128)tmp_q[1] * 19937002142885L) - ((int128)tmp_q[2] * 11200928273737L) - ((int128)tmp_q[3] * 13967116652263L) - ((int128)tmp_q[4] * 10557391427412L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 10557391427412L) + ((int128)tmp_q[1] * 19535472087335L) + ((-((int128)tmp_q[2] * 19937002142885L) - ((int128)tmp_q[3] * 11200928273737L) - ((int128)tmp_q[4] * 13967116652263L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 13967116652263L) - ((int128)tmp_q[1] * 10557391427412L) + ((int128)tmp_q[2] * 19535472087335L) + ((-((int128)tmp_q[3] * 19937002142885L) - ((int128)tmp_q[4] * 11200928273737L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 11200928273737L) - ((int128)tmp_q[1] * 13967116652263L) - ((int128)tmp_q[2] * 10557391427412L) + ((int128)tmp_q[3] * 19535472087335L) - ((int128)tmp_q[4] * 159496017143080L);
	tmp_zero[4] = -((int128)tmp_q[0] * 19937002142885L) - ((int128)tmp_q[1] * 11200928273737L) - ((int128)tmp_q[2] * 13967116652263L) - ((int128)tmp_q[3] * 10557391427412L) + ((int128)tmp_q[4] * 19535472087335L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

