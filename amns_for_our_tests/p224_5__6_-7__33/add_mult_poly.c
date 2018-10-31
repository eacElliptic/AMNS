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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15975786809541715381UL) + ((((uint64_t)op[1] * 14696727487331019519UL) + ((uint64_t)op[2] * 3946568951690776482UL) + ((uint64_t)op[3] * 18143254459981919277UL) + ((uint64_t)op[4] * 12288905067766293697UL) + ((uint64_t)op[5] * 15381323095617184895UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 15381323095617184895UL) + ((uint64_t)op[1] * 15975786809541715381UL) + ((((uint64_t)op[2] * 14696727487331019519UL) + ((uint64_t)op[3] * 3946568951690776482UL) + ((uint64_t)op[4] * 18143254459981919277UL) + ((uint64_t)op[5] * 12288905067766293697UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 12288905067766293697UL) + ((uint64_t)op[1] * 15381323095617184895UL) + ((uint64_t)op[2] * 15975786809541715381UL) + ((((uint64_t)op[3] * 14696727487331019519UL) + ((uint64_t)op[4] * 3946568951690776482UL) + ((uint64_t)op[5] * 18143254459981919277UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 18143254459981919277UL) + ((uint64_t)op[1] * 12288905067766293697UL) + ((uint64_t)op[2] * 15381323095617184895UL) + ((uint64_t)op[3] * 15975786809541715381UL) + ((((uint64_t)op[4] * 14696727487331019519UL) + ((uint64_t)op[5] * 3946568951690776482UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 3946568951690776482UL) + ((uint64_t)op[1] * 18143254459981919277UL) + ((uint64_t)op[2] * 12288905067766293697UL) + ((uint64_t)op[3] * 15381323095617184895UL) + ((uint64_t)op[4] * 15975786809541715381UL) + ((uint64_t)op[5] * 7803372030940173063UL);
	tmp_q[5] = ((uint64_t)op[0] * 14696727487331019519UL) + ((uint64_t)op[1] * 3946568951690776482UL) + ((uint64_t)op[2] * 18143254459981919277UL) + ((uint64_t)op[3] * 12288905067766293697UL) + ((uint64_t)op[4] * 15381323095617184895UL) + ((uint64_t)op[5] * 15975786809541715381UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 58290688697L) - ((((int128)tmp_q[1] * 86710037469L) - ((int128)tmp_q[2] * 21305411001L) + ((int128)tmp_q[3] * 2354731821L) - ((int128)tmp_q[4] * 10038939054L) + ((int128)tmp_q[5] * 109314627787L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 109314627787L) - ((int128)tmp_q[1] * 58290688697L) - ((((int128)tmp_q[2] * 86710037469L) - ((int128)tmp_q[3] * 21305411001L) + ((int128)tmp_q[4] * 2354731821L) - ((int128)tmp_q[5] * 10038939054L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 10038939054L) + ((int128)tmp_q[1] * 109314627787L) - ((int128)tmp_q[2] * 58290688697L) - ((((int128)tmp_q[3] * 86710037469L) - ((int128)tmp_q[4] * 21305411001L) + ((int128)tmp_q[5] * 2354731821L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 2354731821L) - ((int128)tmp_q[1] * 10038939054L) + ((int128)tmp_q[2] * 109314627787L) - ((int128)tmp_q[3] * 58290688697L) - ((((int128)tmp_q[4] * 86710037469L) - ((int128)tmp_q[5] * 21305411001L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 21305411001L) + ((int128)tmp_q[1] * 2354731821L) - ((int128)tmp_q[2] * 10038939054L) + ((int128)tmp_q[3] * 109314627787L) - ((int128)tmp_q[4] * 58290688697L) - ((int128)tmp_q[5] * 606970262283L);
	tmp_zero[5] = ((int128)tmp_q[0] * 86710037469L) - ((int128)tmp_q[1] * 21305411001L) + ((int128)tmp_q[2] * 2354731821L) - ((int128)tmp_q[3] * 10038939054L) + ((int128)tmp_q[4] * 109314627787L) - ((int128)tmp_q[5] * 58290688697L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

