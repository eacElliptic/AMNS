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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10051760207667066366UL) + ((((uint64_t)op[1] * 668954469913936017UL) + ((uint64_t)op[2] * 12717675867016249724UL) + ((uint64_t)op[3] * 5358081389780012818UL) + ((uint64_t)op[4] * 1412937687757721800UL) + ((uint64_t)op[5] * 7111436683047240294UL) + ((uint64_t)op[6] * 1032857142503917066UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 1032857142503917066UL) + ((uint64_t)op[1] * 10051760207667066366UL) + ((((uint64_t)op[2] * 668954469913936017UL) + ((uint64_t)op[3] * 12717675867016249724UL) + ((uint64_t)op[4] * 5358081389780012818UL) + ((uint64_t)op[5] * 1412937687757721800UL) + ((uint64_t)op[6] * 7111436683047240294UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 7111436683047240294UL) + ((uint64_t)op[1] * 1032857142503917066UL) + ((uint64_t)op[2] * 10051760207667066366UL) + ((((uint64_t)op[3] * 668954469913936017UL) + ((uint64_t)op[4] * 12717675867016249724UL) + ((uint64_t)op[5] * 5358081389780012818UL) + ((uint64_t)op[6] * 1412937687757721800UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 1412937687757721800UL) + ((uint64_t)op[1] * 7111436683047240294UL) + ((uint64_t)op[2] * 1032857142503917066UL) + ((uint64_t)op[3] * 10051760207667066366UL) + ((((uint64_t)op[4] * 668954469913936017UL) + ((uint64_t)op[5] * 12717675867016249724UL) + ((uint64_t)op[6] * 5358081389780012818UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 5358081389780012818UL) + ((uint64_t)op[1] * 1412937687757721800UL) + ((uint64_t)op[2] * 7111436683047240294UL) + ((uint64_t)op[3] * 1032857142503917066UL) + ((uint64_t)op[4] * 10051760207667066366UL) + ((((uint64_t)op[5] * 668954469913936017UL) + ((uint64_t)op[6] * 12717675867016249724UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 12717675867016249724UL) + ((uint64_t)op[1] * 5358081389780012818UL) + ((uint64_t)op[2] * 1412937687757721800UL) + ((uint64_t)op[3] * 7111436683047240294UL) + ((uint64_t)op[4] * 1032857142503917066UL) + ((uint64_t)op[5] * 10051760207667066366UL) + ((uint64_t)op[6] * 16439880663967743565UL);
	tmp_q[6] = ((uint64_t)op[0] * 668954469913936017UL) + ((uint64_t)op[1] * 12717675867016249724UL) + ((uint64_t)op[2] * 5358081389780012818UL) + ((uint64_t)op[3] * 1412937687757721800UL) + ((uint64_t)op[4] * 7111436683047240294UL) + ((uint64_t)op[5] * 1032857142503917066UL) + ((uint64_t)op[6] * 10051760207667066366UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 23167563488L) - ((((int128)tmp_q[1] * 32271461154L) + ((int128)tmp_q[2] * 12412863108L) - ((int128)tmp_q[3] * 1305995406L) + ((int128)tmp_q[4] * 25688291750L) - ((int128)tmp_q[5] * 47006974122L) + ((int128)tmp_q[6] * 3792317803L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 3792317803L) - ((int128)tmp_q[1] * 23167563488L) - ((((int128)tmp_q[2] * 32271461154L) + ((int128)tmp_q[3] * 12412863108L) - ((int128)tmp_q[4] * 1305995406L) + ((int128)tmp_q[5] * 25688291750L) - ((int128)tmp_q[6] * 47006974122L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 47006974122L) + ((int128)tmp_q[1] * 3792317803L) - ((int128)tmp_q[2] * 23167563488L) - ((((int128)tmp_q[3] * 32271461154L) + ((int128)tmp_q[4] * 12412863108L) - ((int128)tmp_q[5] * 1305995406L) + ((int128)tmp_q[6] * 25688291750L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 25688291750L) - ((int128)tmp_q[1] * 47006974122L) + ((int128)tmp_q[2] * 3792317803L) - ((int128)tmp_q[3] * 23167563488L) - ((((int128)tmp_q[4] * 32271461154L) + ((int128)tmp_q[5] * 12412863108L) - ((int128)tmp_q[6] * 1305995406L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 1305995406L) + ((int128)tmp_q[1] * 25688291750L) - ((int128)tmp_q[2] * 47006974122L) + ((int128)tmp_q[3] * 3792317803L) - ((int128)tmp_q[4] * 23167563488L) - ((((int128)tmp_q[5] * 32271461154L) + ((int128)tmp_q[6] * 12412863108L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 12412863108L) - ((int128)tmp_q[1] * 1305995406L) + ((int128)tmp_q[2] * 25688291750L) - ((int128)tmp_q[3] * 47006974122L) + ((int128)tmp_q[4] * 3792317803L) - ((int128)tmp_q[5] * 23167563488L) - ((int128)tmp_q[6] * 96814383462L);
	tmp_zero[6] = ((int128)tmp_q[0] * 32271461154L) + ((int128)tmp_q[1] * 12412863108L) - ((int128)tmp_q[2] * 1305995406L) + ((int128)tmp_q[3] * 25688291750L) - ((int128)tmp_q[4] * 47006974122L) + ((int128)tmp_q[5] * 3792317803L) - ((int128)tmp_q[6] * 23167563488L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

