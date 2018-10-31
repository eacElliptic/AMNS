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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11964839494229500147UL) + ((((uint64_t)op[1] * 17623898660730507570UL) + ((uint64_t)op[2] * 12866215149553484756UL) + ((uint64_t)op[3] * 4131397615548142776UL) + ((uint64_t)op[4] * 1693389385148690543UL) + ((uint64_t)op[5] * 15451408093375436818UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 15451408093375436818UL) + ((uint64_t)op[1] * 11964839494229500147UL) + ((((uint64_t)op[2] * 17623898660730507570UL) + ((uint64_t)op[3] * 12866215149553484756UL) + ((uint64_t)op[4] * 4131397615548142776UL) + ((uint64_t)op[5] * 1693389385148690543UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 1693389385148690543UL) + ((uint64_t)op[1] * 15451408093375436818UL) + ((uint64_t)op[2] * 11964839494229500147UL) + ((((uint64_t)op[3] * 17623898660730507570UL) + ((uint64_t)op[4] * 12866215149553484756UL) + ((uint64_t)op[5] * 4131397615548142776UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 4131397615548142776UL) + ((uint64_t)op[1] * 1693389385148690543UL) + ((uint64_t)op[2] * 15451408093375436818UL) + ((uint64_t)op[3] * 11964839494229500147UL) + ((((uint64_t)op[4] * 17623898660730507570UL) + ((uint64_t)op[5] * 12866215149553484756UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 12866215149553484756UL) + ((uint64_t)op[1] * 4131397615548142776UL) + ((uint64_t)op[2] * 1693389385148690543UL) + ((uint64_t)op[3] * 15451408093375436818UL) + ((uint64_t)op[4] * 11964839494229500147UL) + ((uint64_t)op[5] * 11863980769877199248UL);
	tmp_q[5] = ((uint64_t)op[0] * 17623898660730507570UL) + ((uint64_t)op[1] * 12866215149553484756UL) + ((uint64_t)op[2] * 4131397615548142776UL) + ((uint64_t)op[3] * 1693389385148690543UL) + ((uint64_t)op[4] * 15451408093375436818UL) + ((uint64_t)op[5] * 11964839494229500147UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 20472959415229L) + ((((int128)tmp_q[1] * 17801344128152L) + ((int128)tmp_q[2] * 33038828639085L) + ((int128)tmp_q[3] * 56327611140588L) - ((int128)tmp_q[4] * 57465333149997L) - ((int128)tmp_q[5] * 32614141602334L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 32614141602334L) + ((int128)tmp_q[1] * 20472959415229L) + ((((int128)tmp_q[2] * 17801344128152L) + ((int128)tmp_q[3] * 33038828639085L) + ((int128)tmp_q[4] * 56327611140588L) - ((int128)tmp_q[5] * 57465333149997L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 57465333149997L) - ((int128)tmp_q[1] * 32614141602334L) + ((int128)tmp_q[2] * 20472959415229L) + ((((int128)tmp_q[3] * 17801344128152L) + ((int128)tmp_q[4] * 33038828639085L) + ((int128)tmp_q[5] * 56327611140588L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 56327611140588L) - ((int128)tmp_q[1] * 57465333149997L) - ((int128)tmp_q[2] * 32614141602334L) + ((int128)tmp_q[3] * 20472959415229L) + ((((int128)tmp_q[4] * 17801344128152L) + ((int128)tmp_q[5] * 33038828639085L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 33038828639085L) + ((int128)tmp_q[1] * 56327611140588L) - ((int128)tmp_q[2] * 57465333149997L) - ((int128)tmp_q[3] * 32614141602334L) + ((int128)tmp_q[4] * 20472959415229L) + ((int128)tmp_q[5] * 142410753025216L);
	tmp_zero[5] = ((int128)tmp_q[0] * 17801344128152L) + ((int128)tmp_q[1] * 33038828639085L) + ((int128)tmp_q[2] * 56327611140588L) - ((int128)tmp_q[3] * 57465333149997L) - ((int128)tmp_q[4] * 32614141602334L) + ((int128)tmp_q[5] * 20472959415229L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

