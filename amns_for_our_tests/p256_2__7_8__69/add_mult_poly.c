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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14590047459663261049UL) + ((((uint64_t)op[1] * 14045452291257871415UL) + ((uint64_t)op[2] * 10123493666044419409UL) + ((uint64_t)op[3] * 5364523226126592619UL) + ((uint64_t)op[4] * 2907855036599398910UL) + ((uint64_t)op[5] * 8082187280395474432UL) + ((uint64_t)op[6] * 17871144740696341532UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 17871144740696341532UL) + ((uint64_t)op[1] * 14590047459663261049UL) + ((((uint64_t)op[2] * 14045452291257871415UL) + ((uint64_t)op[3] * 10123493666044419409UL) + ((uint64_t)op[4] * 5364523226126592619UL) + ((uint64_t)op[5] * 2907855036599398910UL) + ((uint64_t)op[6] * 8082187280395474432UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 8082187280395474432UL) + ((uint64_t)op[1] * 17871144740696341532UL) + ((uint64_t)op[2] * 14590047459663261049UL) + ((((uint64_t)op[3] * 14045452291257871415UL) + ((uint64_t)op[4] * 10123493666044419409UL) + ((uint64_t)op[5] * 5364523226126592619UL) + ((uint64_t)op[6] * 2907855036599398910UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 2907855036599398910UL) + ((uint64_t)op[1] * 8082187280395474432UL) + ((uint64_t)op[2] * 17871144740696341532UL) + ((uint64_t)op[3] * 14590047459663261049UL) + ((((uint64_t)op[4] * 14045452291257871415UL) + ((uint64_t)op[5] * 10123493666044419409UL) + ((uint64_t)op[6] * 5364523226126592619UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 5364523226126592619UL) + ((uint64_t)op[1] * 2907855036599398910UL) + ((uint64_t)op[2] * 8082187280395474432UL) + ((uint64_t)op[3] * 17871144740696341532UL) + ((uint64_t)op[4] * 14590047459663261049UL) + ((((uint64_t)op[5] * 14045452291257871415UL) + ((uint64_t)op[6] * 10123493666044419409UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 10123493666044419409UL) + ((uint64_t)op[1] * 5364523226126592619UL) + ((uint64_t)op[2] * 2907855036599398910UL) + ((uint64_t)op[3] * 8082187280395474432UL) + ((uint64_t)op[4] * 17871144740696341532UL) + ((uint64_t)op[5] * 14590047459663261049UL) + ((uint64_t)op[6] * 1683153887805661624UL);
	tmp_q[6] = ((uint64_t)op[0] * 14045452291257871415UL) + ((uint64_t)op[1] * 10123493666044419409UL) + ((uint64_t)op[2] * 5364523226126592619UL) + ((uint64_t)op[3] * 2907855036599398910UL) + ((uint64_t)op[4] * 8082187280395474432UL) + ((uint64_t)op[5] * 17871144740696341532UL) + ((uint64_t)op[6] * 14590047459663261049UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 10865520169L) + ((((int128)tmp_q[1] * 42982560643L) - ((int128)tmp_q[2] * 41868922583L) - ((int128)tmp_q[3] * 23591280341L) - ((int128)tmp_q[4] * 11835020506L) + ((int128)tmp_q[5] * 8862285120L) - ((int128)tmp_q[6] * 4678042444L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 4678042444L) - ((int128)tmp_q[1] * 10865520169L) + ((((int128)tmp_q[2] * 42982560643L) - ((int128)tmp_q[3] * 41868922583L) - ((int128)tmp_q[4] * 23591280341L) - ((int128)tmp_q[5] * 11835020506L) + ((int128)tmp_q[6] * 8862285120L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 8862285120L) - ((int128)tmp_q[1] * 4678042444L) - ((int128)tmp_q[2] * 10865520169L) + ((((int128)tmp_q[3] * 42982560643L) - ((int128)tmp_q[4] * 41868922583L) - ((int128)tmp_q[5] * 23591280341L) - ((int128)tmp_q[6] * 11835020506L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 11835020506L) + ((int128)tmp_q[1] * 8862285120L) - ((int128)tmp_q[2] * 4678042444L) - ((int128)tmp_q[3] * 10865520169L) + ((((int128)tmp_q[4] * 42982560643L) - ((int128)tmp_q[5] * 41868922583L) - ((int128)tmp_q[6] * 23591280341L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 23591280341L) - ((int128)tmp_q[1] * 11835020506L) + ((int128)tmp_q[2] * 8862285120L) - ((int128)tmp_q[3] * 4678042444L) - ((int128)tmp_q[4] * 10865520169L) + ((((int128)tmp_q[5] * 42982560643L) - ((int128)tmp_q[6] * 41868922583L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 41868922583L) - ((int128)tmp_q[1] * 23591280341L) - ((int128)tmp_q[2] * 11835020506L) + ((int128)tmp_q[3] * 8862285120L) - ((int128)tmp_q[4] * 4678042444L) - ((int128)tmp_q[5] * 10865520169L) + ((int128)tmp_q[6] * 343860485144L);
	tmp_zero[6] = ((int128)tmp_q[0] * 42982560643L) - ((int128)tmp_q[1] * 41868922583L) - ((int128)tmp_q[2] * 23591280341L) - ((int128)tmp_q[3] * 11835020506L) + ((int128)tmp_q[4] * 8862285120L) - ((int128)tmp_q[5] * 4678042444L) - ((int128)tmp_q[6] * 10865520169L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

