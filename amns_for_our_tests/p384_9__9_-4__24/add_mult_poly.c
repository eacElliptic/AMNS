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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4945503768996424777UL) + ((((uint64_t)op[1] * 8075837279414788818UL) + ((uint64_t)op[2] * 15189432189208986333UL) + ((uint64_t)op[3] * 10526958009101403359UL) + ((uint64_t)op[4] * 13192951823810783017UL) + ((uint64_t)op[5] * 12669724166897788449UL) + ((uint64_t)op[6] * 14287769930603504455UL) + ((uint64_t)op[7] * 13389380990140826448UL) + ((uint64_t)op[8] * 5883581066734317561UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 5883581066734317561UL) + ((uint64_t)op[1] * 4945503768996424777UL) + ((((uint64_t)op[2] * 8075837279414788818UL) + ((uint64_t)op[3] * 15189432189208986333UL) + ((uint64_t)op[4] * 10526958009101403359UL) + ((uint64_t)op[5] * 13192951823810783017UL) + ((uint64_t)op[6] * 12669724166897788449UL) + ((uint64_t)op[7] * 14287769930603504455UL) + ((uint64_t)op[8] * 13389380990140826448UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 13389380990140826448UL) + ((uint64_t)op[1] * 5883581066734317561UL) + ((uint64_t)op[2] * 4945503768996424777UL) + ((((uint64_t)op[3] * 8075837279414788818UL) + ((uint64_t)op[4] * 15189432189208986333UL) + ((uint64_t)op[5] * 10526958009101403359UL) + ((uint64_t)op[6] * 13192951823810783017UL) + ((uint64_t)op[7] * 12669724166897788449UL) + ((uint64_t)op[8] * 14287769930603504455UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 14287769930603504455UL) + ((uint64_t)op[1] * 13389380990140826448UL) + ((uint64_t)op[2] * 5883581066734317561UL) + ((uint64_t)op[3] * 4945503768996424777UL) + ((((uint64_t)op[4] * 8075837279414788818UL) + ((uint64_t)op[5] * 15189432189208986333UL) + ((uint64_t)op[6] * 10526958009101403359UL) + ((uint64_t)op[7] * 13192951823810783017UL) + ((uint64_t)op[8] * 12669724166897788449UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 12669724166897788449UL) + ((uint64_t)op[1] * 14287769930603504455UL) + ((uint64_t)op[2] * 13389380990140826448UL) + ((uint64_t)op[3] * 5883581066734317561UL) + ((uint64_t)op[4] * 4945503768996424777UL) + ((((uint64_t)op[5] * 8075837279414788818UL) + ((uint64_t)op[6] * 15189432189208986333UL) + ((uint64_t)op[7] * 10526958009101403359UL) + ((uint64_t)op[8] * 13192951823810783017UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 13192951823810783017UL) + ((uint64_t)op[1] * 12669724166897788449UL) + ((uint64_t)op[2] * 14287769930603504455UL) + ((uint64_t)op[3] * 13389380990140826448UL) + ((uint64_t)op[4] * 5883581066734317561UL) + ((uint64_t)op[5] * 4945503768996424777UL) + ((((uint64_t)op[6] * 8075837279414788818UL) + ((uint64_t)op[7] * 15189432189208986333UL) + ((uint64_t)op[8] * 10526958009101403359UL)) * 18446744073709551612);
	tmp_q[6] = ((uint64_t)op[0] * 10526958009101403359UL) + ((uint64_t)op[1] * 13192951823810783017UL) + ((uint64_t)op[2] * 12669724166897788449UL) + ((uint64_t)op[3] * 14287769930603504455UL) + ((uint64_t)op[4] * 13389380990140826448UL) + ((uint64_t)op[5] * 5883581066734317561UL) + ((uint64_t)op[6] * 4945503768996424777UL) + ((((uint64_t)op[7] * 8075837279414788818UL) + ((uint64_t)op[8] * 15189432189208986333UL)) * 18446744073709551612);
	tmp_q[7] = ((uint64_t)op[0] * 15189432189208986333UL) + ((uint64_t)op[1] * 10526958009101403359UL) + ((uint64_t)op[2] * 13192951823810783017UL) + ((uint64_t)op[3] * 12669724166897788449UL) + ((uint64_t)op[4] * 14287769930603504455UL) + ((uint64_t)op[5] * 13389380990140826448UL) + ((uint64_t)op[6] * 5883581066734317561UL) + ((uint64_t)op[7] * 4945503768996424777UL) + ((uint64_t)op[8] * 4590139029759947960UL);
	tmp_q[8] = ((uint64_t)op[0] * 8075837279414788818UL) + ((uint64_t)op[1] * 15189432189208986333UL) + ((uint64_t)op[2] * 10526958009101403359UL) + ((uint64_t)op[3] * 13192951823810783017UL) + ((uint64_t)op[4] * 12669724166897788449UL) + ((uint64_t)op[5] * 14287769930603504455UL) + ((uint64_t)op[6] * 13389380990140826448UL) + ((uint64_t)op[7] * 5883581066734317561UL) + ((uint64_t)op[8] * 4945503768996424777UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 85838361729L) - ((((int128)tmp_q[1] * 497171087448L) - ((int128)tmp_q[2] * 2318628355821L) - ((int128)tmp_q[3] * 957770443086L) - ((int128)tmp_q[4] * 1001866831091L) + ((int128)tmp_q[5] * 1201610799054L) + ((int128)tmp_q[6] * 4306080210736L) + ((int128)tmp_q[7] * 125700476951L) + ((int128)tmp_q[8] * 1031204949857L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 1031204949857L) - ((int128)tmp_q[1] * 85838361729L) - ((((int128)tmp_q[2] * 497171087448L) - ((int128)tmp_q[3] * 2318628355821L) - ((int128)tmp_q[4] * 957770443086L) - ((int128)tmp_q[5] * 1001866831091L) + ((int128)tmp_q[6] * 1201610799054L) + ((int128)tmp_q[7] * 4306080210736L) + ((int128)tmp_q[8] * 125700476951L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 125700476951L) + ((int128)tmp_q[1] * 1031204949857L) - ((int128)tmp_q[2] * 85838361729L) - ((((int128)tmp_q[3] * 497171087448L) - ((int128)tmp_q[4] * 2318628355821L) - ((int128)tmp_q[5] * 957770443086L) - ((int128)tmp_q[6] * 1001866831091L) + ((int128)tmp_q[7] * 1201610799054L) + ((int128)tmp_q[8] * 4306080210736L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 4306080210736L) + ((int128)tmp_q[1] * 125700476951L) + ((int128)tmp_q[2] * 1031204949857L) - ((int128)tmp_q[3] * 85838361729L) - ((((int128)tmp_q[4] * 497171087448L) - ((int128)tmp_q[5] * 2318628355821L) - ((int128)tmp_q[6] * 957770443086L) - ((int128)tmp_q[7] * 1001866831091L) + ((int128)tmp_q[8] * 1201610799054L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 1201610799054L) + ((int128)tmp_q[1] * 4306080210736L) + ((int128)tmp_q[2] * 125700476951L) + ((int128)tmp_q[3] * 1031204949857L) - ((int128)tmp_q[4] * 85838361729L) - ((((int128)tmp_q[5] * 497171087448L) - ((int128)tmp_q[6] * 2318628355821L) - ((int128)tmp_q[7] * 957770443086L) - ((int128)tmp_q[8] * 1001866831091L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 1001866831091L) + ((int128)tmp_q[1] * 1201610799054L) + ((int128)tmp_q[2] * 4306080210736L) + ((int128)tmp_q[3] * 125700476951L) + ((int128)tmp_q[4] * 1031204949857L) - ((int128)tmp_q[5] * 85838361729L) - ((((int128)tmp_q[6] * 497171087448L) - ((int128)tmp_q[7] * 2318628355821L) - ((int128)tmp_q[8] * 957770443086L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 957770443086L) - ((int128)tmp_q[1] * 1001866831091L) + ((int128)tmp_q[2] * 1201610799054L) + ((int128)tmp_q[3] * 4306080210736L) + ((int128)tmp_q[4] * 125700476951L) + ((int128)tmp_q[5] * 1031204949857L) - ((int128)tmp_q[6] * 85838361729L) - ((((int128)tmp_q[7] * 497171087448L) - ((int128)tmp_q[8] * 2318628355821L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 2318628355821L) - ((int128)tmp_q[1] * 957770443086L) - ((int128)tmp_q[2] * 1001866831091L) + ((int128)tmp_q[3] * 1201610799054L) + ((int128)tmp_q[4] * 4306080210736L) + ((int128)tmp_q[5] * 125700476951L) + ((int128)tmp_q[6] * 1031204949857L) - ((int128)tmp_q[7] * 85838361729L) - ((int128)tmp_q[8] * 1988684349792L);
	tmp_zero[8] = ((int128)tmp_q[0] * 497171087448L) - ((int128)tmp_q[1] * 2318628355821L) - ((int128)tmp_q[2] * 957770443086L) - ((int128)tmp_q[3] * 1001866831091L) + ((int128)tmp_q[4] * 1201610799054L) + ((int128)tmp_q[5] * 4306080210736L) + ((int128)tmp_q[6] * 125700476951L) + ((int128)tmp_q[7] * 1031204949857L) - ((int128)tmp_q[8] * 85838361729L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
}

