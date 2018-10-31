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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18244465044316407780UL) + ((((uint64_t)op[1] * 10684728608295967750UL) + ((uint64_t)op[2] * 7458549314044237819UL) + ((uint64_t)op[3] * 6691589402014324214UL) + ((uint64_t)op[4] * 586249977934931298UL) + ((uint64_t)op[5] * 14812327051859616540UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 14812327051859616540UL) + ((uint64_t)op[1] * 18244465044316407780UL) + ((((uint64_t)op[2] * 10684728608295967750UL) + ((uint64_t)op[3] * 7458549314044237819UL) + ((uint64_t)op[4] * 6691589402014324214UL) + ((uint64_t)op[5] * 586249977934931298UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 586249977934931298UL) + ((uint64_t)op[1] * 14812327051859616540UL) + ((uint64_t)op[2] * 18244465044316407780UL) + ((((uint64_t)op[3] * 10684728608295967750UL) + ((uint64_t)op[4] * 7458549314044237819UL) + ((uint64_t)op[5] * 6691589402014324214UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 6691589402014324214UL) + ((uint64_t)op[1] * 586249977934931298UL) + ((uint64_t)op[2] * 14812327051859616540UL) + ((uint64_t)op[3] * 18244465044316407780UL) + ((((uint64_t)op[4] * 10684728608295967750UL) + ((uint64_t)op[5] * 7458549314044237819UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 7458549314044237819UL) + ((uint64_t)op[1] * 6691589402014324214UL) + ((uint64_t)op[2] * 586249977934931298UL) + ((uint64_t)op[3] * 14812327051859616540UL) + ((uint64_t)op[4] * 18244465044316407780UL) + ((uint64_t)op[5] * 1916589179648816098UL);
	tmp_q[5] = ((uint64_t)op[0] * 10684728608295967750UL) + ((uint64_t)op[1] * 7458549314044237819UL) + ((uint64_t)op[2] * 6691589402014324214UL) + ((uint64_t)op[3] * 586249977934931298UL) + ((uint64_t)op[4] * 14812327051859616540UL) + ((uint64_t)op[5] * 18244465044316407780UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2653337919358L) - ((-((int128)tmp_q[1] * 2272996357300L) + ((int128)tmp_q[2] * 2922468602188L) - ((int128)tmp_q[3] * 3022415529166L) + ((int128)tmp_q[4] * 1490114234671L) + ((int128)tmp_q[5] * 589484361746L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 589484361746L) - ((int128)tmp_q[1] * 2653337919358L) - ((-((int128)tmp_q[2] * 2272996357300L) + ((int128)tmp_q[3] * 2922468602188L) - ((int128)tmp_q[4] * 3022415529166L) + ((int128)tmp_q[5] * 1490114234671L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 1490114234671L) + ((int128)tmp_q[1] * 589484361746L) - ((int128)tmp_q[2] * 2653337919358L) - ((-((int128)tmp_q[3] * 2272996357300L) + ((int128)tmp_q[4] * 2922468602188L) - ((int128)tmp_q[5] * 3022415529166L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 3022415529166L) + ((int128)tmp_q[1] * 1490114234671L) + ((int128)tmp_q[2] * 589484361746L) - ((int128)tmp_q[3] * 2653337919358L) - ((-((int128)tmp_q[4] * 2272996357300L) + ((int128)tmp_q[5] * 2922468602188L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 2922468602188L) - ((int128)tmp_q[1] * 3022415529166L) + ((int128)tmp_q[2] * 1490114234671L) + ((int128)tmp_q[3] * 589484361746L) - ((int128)tmp_q[4] * 2653337919358L) + ((int128)tmp_q[5] * 11364981786500L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2272996357300L) + ((int128)tmp_q[1] * 2922468602188L) - ((int128)tmp_q[2] * 3022415529166L) + ((int128)tmp_q[3] * 1490114234671L) + ((int128)tmp_q[4] * 589484361746L) - ((int128)tmp_q[5] * 2653337919358L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

