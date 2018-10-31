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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - ((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - ((int128)pa[4] * pa[4]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3506504934515119952UL) + ((((uint64_t)op[1] * 9293525270395398787UL) + ((uint64_t)op[2] * 5612186825336067021UL) + ((uint64_t)op[3] * 10185096323456504186UL) + ((uint64_t)op[4] * 17372939703030955905UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 17372939703030955905UL) + ((uint64_t)op[1] * 3506504934515119952UL) + ((((uint64_t)op[2] * 9293525270395398787UL) + ((uint64_t)op[3] * 5612186825336067021UL) + ((uint64_t)op[4] * 10185096323456504186UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 10185096323456504186UL) + ((uint64_t)op[1] * 17372939703030955905UL) + ((uint64_t)op[2] * 3506504934515119952UL) + ((((uint64_t)op[3] * 9293525270395398787UL) + ((uint64_t)op[4] * 5612186825336067021UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 5612186825336067021UL) + ((uint64_t)op[1] * 10185096323456504186UL) + ((uint64_t)op[2] * 17372939703030955905UL) + ((uint64_t)op[3] * 3506504934515119952UL) + ((uint64_t)op[4] * 9153218803314152829UL);
	tmp_q[4] = ((uint64_t)op[0] * 9293525270395398787UL) + ((uint64_t)op[1] * 5612186825336067021UL) + ((uint64_t)op[2] * 10185096323456504186UL) + ((uint64_t)op[3] * 17372939703030955905UL) + ((uint64_t)op[4] * 3506504934515119952UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 171633210360383L) - (((int128)tmp_q[1] * 15905447963001L) - ((int128)tmp_q[2] * 17870365844859L) - ((int128)tmp_q[3] * 103037080434948L) + ((int128)tmp_q[4] * 102371943733294L));
	tmp_zero[1] = ((int128)tmp_q[0] * 102371943733294L) + ((int128)tmp_q[1] * 171633210360383L) - (((int128)tmp_q[2] * 15905447963001L) - ((int128)tmp_q[3] * 17870365844859L) - ((int128)tmp_q[4] * 103037080434948L));
	tmp_zero[2] = -((int128)tmp_q[0] * 103037080434948L) + ((int128)tmp_q[1] * 102371943733294L) + ((int128)tmp_q[2] * 171633210360383L) - (((int128)tmp_q[3] * 15905447963001L) - ((int128)tmp_q[4] * 17870365844859L));
	tmp_zero[3] = -((int128)tmp_q[0] * 17870365844859L) - ((int128)tmp_q[1] * 103037080434948L) + ((int128)tmp_q[2] * 102371943733294L) + ((int128)tmp_q[3] * 171633210360383L) - ((int128)tmp_q[4] * 15905447963001L);
	tmp_zero[4] = ((int128)tmp_q[0] * 15905447963001L) - ((int128)tmp_q[1] * 17870365844859L) - ((int128)tmp_q[2] * 103037080434948L) + ((int128)tmp_q[3] * 102371943733294L) + ((int128)tmp_q[4] * 171633210360383L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

