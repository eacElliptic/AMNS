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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11521848772154745860UL) + ((((uint64_t)op[1] * 363292294635378755UL) + ((uint64_t)op[2] * 11357973236365493591UL) + ((uint64_t)op[3] * 3573627728262613876UL) + ((uint64_t)op[4] * 7488666866062775584UL) + ((uint64_t)op[5] * 10605293442572625219UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 10605293442572625219UL) + ((uint64_t)op[1] * 11521848772154745860UL) + ((((uint64_t)op[2] * 363292294635378755UL) + ((uint64_t)op[3] * 11357973236365493591UL) + ((uint64_t)op[4] * 3573627728262613876UL) + ((uint64_t)op[5] * 7488666866062775584UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 7488666866062775584UL) + ((uint64_t)op[1] * 10605293442572625219UL) + ((uint64_t)op[2] * 11521848772154745860UL) + ((((uint64_t)op[3] * 363292294635378755UL) + ((uint64_t)op[4] * 11357973236365493591UL) + ((uint64_t)op[5] * 3573627728262613876UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 3573627728262613876UL) + ((uint64_t)op[1] * 7488666866062775584UL) + ((uint64_t)op[2] * 10605293442572625219UL) + ((uint64_t)op[3] * 11521848772154745860UL) + ((((uint64_t)op[4] * 363292294635378755UL) + ((uint64_t)op[5] * 11357973236365493591UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 11357973236365493591UL) + ((uint64_t)op[1] * 3573627728262613876UL) + ((uint64_t)op[2] * 7488666866062775584UL) + ((uint64_t)op[3] * 10605293442572625219UL) + ((uint64_t)op[4] * 11521848772154745860UL) + ((uint64_t)op[5] * 17356867189803415351UL);
	tmp_q[5] = ((uint64_t)op[0] * 363292294635378755UL) + ((uint64_t)op[1] * 11357973236365493591UL) + ((uint64_t)op[2] * 3573627728262613876UL) + ((uint64_t)op[3] * 7488666866062775584UL) + ((uint64_t)op[4] * 10605293442572625219UL) + ((uint64_t)op[5] * 11521848772154745860UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4537927349597L) - ((((int128)tmp_q[1] * 1280715661706L) + ((int128)tmp_q[2] * 909092519784L) - ((int128)tmp_q[3] * 358228643669L) + ((int128)tmp_q[4] * 2716673816166L) + ((int128)tmp_q[5] * 1108392432697L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 1108392432697L) - ((int128)tmp_q[1] * 4537927349597L) - ((((int128)tmp_q[2] * 1280715661706L) + ((int128)tmp_q[3] * 909092519784L) - ((int128)tmp_q[4] * 358228643669L) + ((int128)tmp_q[5] * 2716673816166L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2716673816166L) + ((int128)tmp_q[1] * 1108392432697L) - ((int128)tmp_q[2] * 4537927349597L) - ((((int128)tmp_q[3] * 1280715661706L) + ((int128)tmp_q[4] * 909092519784L) - ((int128)tmp_q[5] * 358228643669L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 358228643669L) + ((int128)tmp_q[1] * 2716673816166L) + ((int128)tmp_q[2] * 1108392432697L) - ((int128)tmp_q[3] * 4537927349597L) - ((((int128)tmp_q[4] * 1280715661706L) + ((int128)tmp_q[5] * 909092519784L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 909092519784L) - ((int128)tmp_q[1] * 358228643669L) + ((int128)tmp_q[2] * 2716673816166L) + ((int128)tmp_q[3] * 1108392432697L) - ((int128)tmp_q[4] * 4537927349597L) - ((int128)tmp_q[5] * 3842146985118L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1280715661706L) + ((int128)tmp_q[1] * 909092519784L) - ((int128)tmp_q[2] * 358228643669L) + ((int128)tmp_q[3] * 2716673816166L) + ((int128)tmp_q[4] * 1108392432697L) - ((int128)tmp_q[5] * 4537927349597L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

