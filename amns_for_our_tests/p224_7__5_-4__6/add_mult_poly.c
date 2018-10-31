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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7898607226691565149UL) + ((((uint64_t)op[1] * 15159032300605020280UL) + ((uint64_t)op[2] * 7852265130802090293UL) + ((uint64_t)op[3] * 6406503754976613627UL) + ((uint64_t)op[4] * 5555820543469565491UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 5555820543469565491UL) + ((uint64_t)op[1] * 7898607226691565149UL) + ((((uint64_t)op[2] * 15159032300605020280UL) + ((uint64_t)op[3] * 7852265130802090293UL) + ((uint64_t)op[4] * 6406503754976613627UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 6406503754976613627UL) + ((uint64_t)op[1] * 5555820543469565491UL) + ((uint64_t)op[2] * 7898607226691565149UL) + ((((uint64_t)op[3] * 15159032300605020280UL) + ((uint64_t)op[4] * 7852265130802090293UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 7852265130802090293UL) + ((uint64_t)op[1] * 6406503754976613627UL) + ((uint64_t)op[2] * 5555820543469565491UL) + ((uint64_t)op[3] * 7898607226691565149UL) + ((uint64_t)op[4] * 13150847092418125344UL);
	tmp_q[4] = ((uint64_t)op[0] * 15159032300605020280UL) + ((uint64_t)op[1] * 7852265130802090293UL) + ((uint64_t)op[2] * 6406503754976613627UL) + ((uint64_t)op[3] * 5555820543469565491UL) + ((uint64_t)op[4] * 7898607226691565149UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 8466612363527L) - ((((int128)tmp_q[1] * 1179149378149L) - ((int128)tmp_q[2] * 17015445850718L) + ((int128)tmp_q[3] * 680894402030L) - ((int128)tmp_q[4] * 4397556523349L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 4397556523349L) + ((int128)tmp_q[1] * 8466612363527L) - ((((int128)tmp_q[2] * 1179149378149L) - ((int128)tmp_q[3] * 17015445850718L) + ((int128)tmp_q[4] * 680894402030L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 680894402030L) - ((int128)tmp_q[1] * 4397556523349L) + ((int128)tmp_q[2] * 8466612363527L) - ((((int128)tmp_q[3] * 1179149378149L) - ((int128)tmp_q[4] * 17015445850718L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 17015445850718L) + ((int128)tmp_q[1] * 680894402030L) - ((int128)tmp_q[2] * 4397556523349L) + ((int128)tmp_q[3] * 8466612363527L) - ((int128)tmp_q[4] * 4716597512596L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1179149378149L) - ((int128)tmp_q[1] * 17015445850718L) + ((int128)tmp_q[2] * 680894402030L) - ((int128)tmp_q[3] * 4397556523349L) + ((int128)tmp_q[4] * 8466612363527L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

