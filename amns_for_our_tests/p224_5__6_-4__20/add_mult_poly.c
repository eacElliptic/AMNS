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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5122616443260675607UL) + ((((uint64_t)op[1] * 17452623355435216710UL) + ((uint64_t)op[2] * 11521920147059425938UL) + ((uint64_t)op[3] * 12159503677132982685UL) + ((uint64_t)op[4] * 8030552450479092808UL) + ((uint64_t)op[5] * 2517006111093959832UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 2517006111093959832UL) + ((uint64_t)op[1] * 5122616443260675607UL) + ((((uint64_t)op[2] * 17452623355435216710UL) + ((uint64_t)op[3] * 11521920147059425938UL) + ((uint64_t)op[4] * 12159503677132982685UL) + ((uint64_t)op[5] * 8030552450479092808UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 8030552450479092808UL) + ((uint64_t)op[1] * 2517006111093959832UL) + ((uint64_t)op[2] * 5122616443260675607UL) + ((((uint64_t)op[3] * 17452623355435216710UL) + ((uint64_t)op[4] * 11521920147059425938UL) + ((uint64_t)op[5] * 12159503677132982685UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 12159503677132982685UL) + ((uint64_t)op[1] * 8030552450479092808UL) + ((uint64_t)op[2] * 2517006111093959832UL) + ((uint64_t)op[3] * 5122616443260675607UL) + ((((uint64_t)op[4] * 17452623355435216710UL) + ((uint64_t)op[5] * 11521920147059425938UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 11521920147059425938UL) + ((uint64_t)op[1] * 12159503677132982685UL) + ((uint64_t)op[2] * 8030552450479092808UL) + ((uint64_t)op[3] * 2517006111093959832UL) + ((uint64_t)op[4] * 5122616443260675607UL) + ((uint64_t)op[5] * 3976482873097339624UL);
	tmp_q[5] = ((uint64_t)op[0] * 17452623355435216710UL) + ((uint64_t)op[1] * 11521920147059425938UL) + ((uint64_t)op[2] * 12159503677132982685UL) + ((uint64_t)op[3] * 8030552450479092808UL) + ((uint64_t)op[4] * 2517006111093959832UL) + ((uint64_t)op[5] * 5122616443260675607UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 92168746597L) - ((-((int128)tmp_q[1] * 9069659714L) + ((int128)tmp_q[2] * 75577818650L) - ((int128)tmp_q[3] * 62982214679L) + ((int128)tmp_q[4] * 93664826792L) + ((int128)tmp_q[5] * 23180817640L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 23180817640L) + ((int128)tmp_q[1] * 92168746597L) - ((-((int128)tmp_q[2] * 9069659714L) + ((int128)tmp_q[3] * 75577818650L) - ((int128)tmp_q[4] * 62982214679L) + ((int128)tmp_q[5] * 93664826792L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 93664826792L) + ((int128)tmp_q[1] * 23180817640L) + ((int128)tmp_q[2] * 92168746597L) - ((-((int128)tmp_q[3] * 9069659714L) + ((int128)tmp_q[4] * 75577818650L) - ((int128)tmp_q[5] * 62982214679L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 62982214679L) + ((int128)tmp_q[1] * 93664826792L) + ((int128)tmp_q[2] * 23180817640L) + ((int128)tmp_q[3] * 92168746597L) - ((-((int128)tmp_q[4] * 9069659714L) + ((int128)tmp_q[5] * 75577818650L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 75577818650L) - ((int128)tmp_q[1] * 62982214679L) + ((int128)tmp_q[2] * 93664826792L) + ((int128)tmp_q[3] * 23180817640L) + ((int128)tmp_q[4] * 92168746597L) + ((int128)tmp_q[5] * 36278638856L);
	tmp_zero[5] = -((int128)tmp_q[0] * 9069659714L) + ((int128)tmp_q[1] * 75577818650L) - ((int128)tmp_q[2] * 62982214679L) + ((int128)tmp_q[3] * 93664826792L) + ((int128)tmp_q[4] * 23180817640L) + ((int128)tmp_q[5] * 92168746597L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

