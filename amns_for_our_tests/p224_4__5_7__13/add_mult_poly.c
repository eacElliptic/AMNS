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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12412702324121363716UL) + ((((uint64_t)op[1] * 128679142214506023UL) + ((uint64_t)op[2] * 5413062697709583607UL) + ((uint64_t)op[3] * 8438783534870542003UL) + ((uint64_t)op[4] * 4274045165542911316UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 4274045165542911316UL) + ((uint64_t)op[1] * 12412702324121363716UL) + ((((uint64_t)op[2] * 128679142214506023UL) + ((uint64_t)op[3] * 5413062697709583607UL) + ((uint64_t)op[4] * 8438783534870542003UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 8438783534870542003UL) + ((uint64_t)op[1] * 4274045165542911316UL) + ((uint64_t)op[2] * 12412702324121363716UL) + ((((uint64_t)op[3] * 128679142214506023UL) + ((uint64_t)op[4] * 5413062697709583607UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 5413062697709583607UL) + ((uint64_t)op[1] * 8438783534870542003UL) + ((uint64_t)op[2] * 4274045165542911316UL) + ((uint64_t)op[3] * 12412702324121363716UL) + ((uint64_t)op[4] * 900753995501542161UL);
	tmp_q[4] = ((uint64_t)op[0] * 128679142214506023UL) + ((uint64_t)op[1] * 5413062697709583607UL) + ((uint64_t)op[2] * 8438783534870542003UL) + ((uint64_t)op[3] * 4274045165542911316UL) + ((uint64_t)op[4] * 12412702324121363716UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6539267916275L) + ((((int128)tmp_q[1] * 10257280454425L) + ((int128)tmp_q[2] * 2215747244068L) - ((int128)tmp_q[3] * 21010045099253L) + ((int128)tmp_q[4] * 11514511345376L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 11514511345376L) - ((int128)tmp_q[1] * 6539267916275L) + ((((int128)tmp_q[2] * 10257280454425L) + ((int128)tmp_q[3] * 2215747244068L) - ((int128)tmp_q[4] * 21010045099253L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 21010045099253L) + ((int128)tmp_q[1] * 11514511345376L) - ((int128)tmp_q[2] * 6539267916275L) + ((((int128)tmp_q[3] * 10257280454425L) + ((int128)tmp_q[4] * 2215747244068L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 2215747244068L) - ((int128)tmp_q[1] * 21010045099253L) + ((int128)tmp_q[2] * 11514511345376L) - ((int128)tmp_q[3] * 6539267916275L) + ((int128)tmp_q[4] * 71800963180975L);
	tmp_zero[4] = ((int128)tmp_q[0] * 10257280454425L) + ((int128)tmp_q[1] * 2215747244068L) - ((int128)tmp_q[2] * 21010045099253L) + ((int128)tmp_q[3] * 11514511345376L) - ((int128)tmp_q[4] * 6539267916275L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

