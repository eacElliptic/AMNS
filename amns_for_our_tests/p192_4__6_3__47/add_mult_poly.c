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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5377222152849561391UL) + ((((uint64_t)op[1] * 7376731871974918006UL) + ((uint64_t)op[2] * 9291952324899614025UL) + ((uint64_t)op[3] * 9972167510862812653UL) + ((uint64_t)op[4] * 423295322607257698UL) + ((uint64_t)op[5] * 7172333639940692598UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 7172333639940692598UL) + ((uint64_t)op[1] * 5377222152849561391UL) + ((((uint64_t)op[2] * 7376731871974918006UL) + ((uint64_t)op[3] * 9291952324899614025UL) + ((uint64_t)op[4] * 9972167510862812653UL) + ((uint64_t)op[5] * 423295322607257698UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 423295322607257698UL) + ((uint64_t)op[1] * 7172333639940692598UL) + ((uint64_t)op[2] * 5377222152849561391UL) + ((((uint64_t)op[3] * 7376731871974918006UL) + ((uint64_t)op[4] * 9291952324899614025UL) + ((uint64_t)op[5] * 9972167510862812653UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 9972167510862812653UL) + ((uint64_t)op[1] * 423295322607257698UL) + ((uint64_t)op[2] * 7172333639940692598UL) + ((uint64_t)op[3] * 5377222152849561391UL) + ((((uint64_t)op[4] * 7376731871974918006UL) + ((uint64_t)op[5] * 9291952324899614025UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 9291952324899614025UL) + ((uint64_t)op[1] * 9972167510862812653UL) + ((uint64_t)op[2] * 423295322607257698UL) + ((uint64_t)op[3] * 7172333639940692598UL) + ((uint64_t)op[4] * 5377222152849561391UL) + ((uint64_t)op[5] * 3683451542215202402UL);
	tmp_q[5] = ((uint64_t)op[0] * 7376731871974918006UL) + ((uint64_t)op[1] * 9291952324899614025UL) + ((uint64_t)op[2] * 9972167510862812653UL) + ((uint64_t)op[3] * 423295322607257698UL) + ((uint64_t)op[4] * 7172333639940692598UL) + ((uint64_t)op[5] * 5377222152849561391UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1893760190L) + ((-((int128)tmp_q[1] * 138196408L) + ((int128)tmp_q[2] * 1381509111L) - ((int128)tmp_q[3] * 2370008568L) - ((int128)tmp_q[4] * 12227173L) - ((int128)tmp_q[5] * 2461692581L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 2461692581L) - ((int128)tmp_q[1] * 1893760190L) + ((-((int128)tmp_q[2] * 138196408L) + ((int128)tmp_q[3] * 1381509111L) - ((int128)tmp_q[4] * 2370008568L) - ((int128)tmp_q[5] * 12227173L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 12227173L) - ((int128)tmp_q[1] * 2461692581L) - ((int128)tmp_q[2] * 1893760190L) + ((-((int128)tmp_q[3] * 138196408L) + ((int128)tmp_q[4] * 1381509111L) - ((int128)tmp_q[5] * 2370008568L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 2370008568L) - ((int128)tmp_q[1] * 12227173L) - ((int128)tmp_q[2] * 2461692581L) - ((int128)tmp_q[3] * 1893760190L) + ((-((int128)tmp_q[4] * 138196408L) + ((int128)tmp_q[5] * 1381509111L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1381509111L) - ((int128)tmp_q[1] * 2370008568L) - ((int128)tmp_q[2] * 12227173L) - ((int128)tmp_q[3] * 2461692581L) - ((int128)tmp_q[4] * 1893760190L) - ((int128)tmp_q[5] * 414589224L);
	tmp_zero[5] = -((int128)tmp_q[0] * 138196408L) + ((int128)tmp_q[1] * 1381509111L) - ((int128)tmp_q[2] * 2370008568L) - ((int128)tmp_q[3] * 12227173L) - ((int128)tmp_q[4] * 2461692581L) - ((int128)tmp_q[5] * 1893760190L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

