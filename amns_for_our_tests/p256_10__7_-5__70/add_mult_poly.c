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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2384351891331928560UL) + ((((uint64_t)op[1] * 15482461012875313355UL) + ((uint64_t)op[2] * 12598932139699932991UL) + ((uint64_t)op[3] * 5618396267759257589UL) + ((uint64_t)op[4] * 17945256329873930615UL) + ((uint64_t)op[5] * 11353848646088835951UL) + ((uint64_t)op[6] * 8466655766032542604UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 8466655766032542604UL) + ((uint64_t)op[1] * 2384351891331928560UL) + ((((uint64_t)op[2] * 15482461012875313355UL) + ((uint64_t)op[3] * 12598932139699932991UL) + ((uint64_t)op[4] * 5618396267759257589UL) + ((uint64_t)op[5] * 17945256329873930615UL) + ((uint64_t)op[6] * 11353848646088835951UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 11353848646088835951UL) + ((uint64_t)op[1] * 8466655766032542604UL) + ((uint64_t)op[2] * 2384351891331928560UL) + ((((uint64_t)op[3] * 15482461012875313355UL) + ((uint64_t)op[4] * 12598932139699932991UL) + ((uint64_t)op[5] * 5618396267759257589UL) + ((uint64_t)op[6] * 17945256329873930615UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 17945256329873930615UL) + ((uint64_t)op[1] * 11353848646088835951UL) + ((uint64_t)op[2] * 8466655766032542604UL) + ((uint64_t)op[3] * 2384351891331928560UL) + ((((uint64_t)op[4] * 15482461012875313355UL) + ((uint64_t)op[5] * 12598932139699932991UL) + ((uint64_t)op[6] * 5618396267759257589UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 5618396267759257589UL) + ((uint64_t)op[1] * 17945256329873930615UL) + ((uint64_t)op[2] * 11353848646088835951UL) + ((uint64_t)op[3] * 8466655766032542604UL) + ((uint64_t)op[4] * 2384351891331928560UL) + ((((uint64_t)op[5] * 15482461012875313355UL) + ((uint64_t)op[6] * 12598932139699932991UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 12598932139699932991UL) + ((uint64_t)op[1] * 5618396267759257589UL) + ((uint64_t)op[2] * 17945256329873930615UL) + ((uint64_t)op[3] * 11353848646088835951UL) + ((uint64_t)op[4] * 8466655766032542604UL) + ((uint64_t)op[5] * 2384351891331928560UL) + ((uint64_t)op[6] * 14821415304171191305UL);
	tmp_q[6] = ((uint64_t)op[0] * 15482461012875313355UL) + ((uint64_t)op[1] * 12598932139699932991UL) + ((uint64_t)op[2] * 5618396267759257589UL) + ((uint64_t)op[3] * 17945256329873930615UL) + ((uint64_t)op[4] * 11353848646088835951UL) + ((uint64_t)op[5] * 8466655766032542604UL) + ((uint64_t)op[6] * 2384351891331928560UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 63332634334L) - ((-((int128)tmp_q[1] * 26532221058L) - ((int128)tmp_q[2] * 30925788505L) - ((int128)tmp_q[3] * 56250692492L) - ((int128)tmp_q[4] * 29447646197L) - ((int128)tmp_q[5] * 56909327452L) - ((int128)tmp_q[6] * 19569273977L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 19569273977L) - ((int128)tmp_q[1] * 63332634334L) - ((-((int128)tmp_q[2] * 26532221058L) - ((int128)tmp_q[3] * 30925788505L) - ((int128)tmp_q[4] * 56250692492L) - ((int128)tmp_q[5] * 29447646197L) - ((int128)tmp_q[6] * 56909327452L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 56909327452L) - ((int128)tmp_q[1] * 19569273977L) - ((int128)tmp_q[2] * 63332634334L) - ((-((int128)tmp_q[3] * 26532221058L) - ((int128)tmp_q[4] * 30925788505L) - ((int128)tmp_q[5] * 56250692492L) - ((int128)tmp_q[6] * 29447646197L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 29447646197L) - ((int128)tmp_q[1] * 56909327452L) - ((int128)tmp_q[2] * 19569273977L) - ((int128)tmp_q[3] * 63332634334L) - ((-((int128)tmp_q[4] * 26532221058L) - ((int128)tmp_q[5] * 30925788505L) - ((int128)tmp_q[6] * 56250692492L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 56250692492L) - ((int128)tmp_q[1] * 29447646197L) - ((int128)tmp_q[2] * 56909327452L) - ((int128)tmp_q[3] * 19569273977L) - ((int128)tmp_q[4] * 63332634334L) - ((-((int128)tmp_q[5] * 26532221058L) - ((int128)tmp_q[6] * 30925788505L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 30925788505L) - ((int128)tmp_q[1] * 56250692492L) - ((int128)tmp_q[2] * 29447646197L) - ((int128)tmp_q[3] * 56909327452L) - ((int128)tmp_q[4] * 19569273977L) - ((int128)tmp_q[5] * 63332634334L) + ((int128)tmp_q[6] * 132661105290L);
	tmp_zero[6] = -((int128)tmp_q[0] * 26532221058L) - ((int128)tmp_q[1] * 30925788505L) - ((int128)tmp_q[2] * 56250692492L) - ((int128)tmp_q[3] * 29447646197L) - ((int128)tmp_q[4] * 56909327452L) - ((int128)tmp_q[5] * 19569273977L) - ((int128)tmp_q[6] * 63332634334L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

