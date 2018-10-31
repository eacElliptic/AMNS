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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10824601389993097013UL) + ((((uint64_t)op[1] * 13313266261704002868UL) + ((uint64_t)op[2] * 2624163691593157155UL) + ((uint64_t)op[3] * 800543273142175200UL) + ((uint64_t)op[4] * 11939769737105245625UL) + ((uint64_t)op[5] * 3846732329237910775UL) + ((uint64_t)op[6] * 5801362673895653526UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 5801362673895653526UL) + ((uint64_t)op[1] * 10824601389993097013UL) + ((((uint64_t)op[2] * 13313266261704002868UL) + ((uint64_t)op[3] * 2624163691593157155UL) + ((uint64_t)op[4] * 800543273142175200UL) + ((uint64_t)op[5] * 11939769737105245625UL) + ((uint64_t)op[6] * 3846732329237910775UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 3846732329237910775UL) + ((uint64_t)op[1] * 5801362673895653526UL) + ((uint64_t)op[2] * 10824601389993097013UL) + ((((uint64_t)op[3] * 13313266261704002868UL) + ((uint64_t)op[4] * 2624163691593157155UL) + ((uint64_t)op[5] * 800543273142175200UL) + ((uint64_t)op[6] * 11939769737105245625UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 11939769737105245625UL) + ((uint64_t)op[1] * 3846732329237910775UL) + ((uint64_t)op[2] * 5801362673895653526UL) + ((uint64_t)op[3] * 10824601389993097013UL) + ((((uint64_t)op[4] * 13313266261704002868UL) + ((uint64_t)op[5] * 2624163691593157155UL) + ((uint64_t)op[6] * 800543273142175200UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 800543273142175200UL) + ((uint64_t)op[1] * 11939769737105245625UL) + ((uint64_t)op[2] * 3846732329237910775UL) + ((uint64_t)op[3] * 5801362673895653526UL) + ((uint64_t)op[4] * 10824601389993097013UL) + ((((uint64_t)op[5] * 13313266261704002868UL) + ((uint64_t)op[6] * 2624163691593157155UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 2624163691593157155UL) + ((uint64_t)op[1] * 800543273142175200UL) + ((uint64_t)op[2] * 11939769737105245625UL) + ((uint64_t)op[3] * 3846732329237910775UL) + ((uint64_t)op[4] * 5801362673895653526UL) + ((uint64_t)op[5] * 10824601389993097013UL) + ((uint64_t)op[6] * 10266955624011097496UL);
	tmp_q[6] = ((uint64_t)op[0] * 13313266261704002868UL) + ((uint64_t)op[1] * 2624163691593157155UL) + ((uint64_t)op[2] * 800543273142175200UL) + ((uint64_t)op[3] * 11939769737105245625UL) + ((uint64_t)op[4] * 3846732329237910775UL) + ((uint64_t)op[5] * 5801362673895653526UL) + ((uint64_t)op[6] * 10824601389993097013UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 12253194395345975L) - ((-((int128)tmp_q[1] * 2786724442451478L) - ((int128)tmp_q[2] * 7642652496425787L) + ((int128)tmp_q[3] * 14518417682501411L) + ((int128)tmp_q[4] * 3081232561050645L) + ((int128)tmp_q[5] * 12473233215350655L) + ((int128)tmp_q[6] * 6429671150180222L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 6429671150180222L) - ((int128)tmp_q[1] * 12253194395345975L) - ((-((int128)tmp_q[2] * 2786724442451478L) - ((int128)tmp_q[3] * 7642652496425787L) + ((int128)tmp_q[4] * 14518417682501411L) + ((int128)tmp_q[5] * 3081232561050645L) + ((int128)tmp_q[6] * 12473233215350655L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 12473233215350655L) + ((int128)tmp_q[1] * 6429671150180222L) - ((int128)tmp_q[2] * 12253194395345975L) - ((-((int128)tmp_q[3] * 2786724442451478L) - ((int128)tmp_q[4] * 7642652496425787L) + ((int128)tmp_q[5] * 14518417682501411L) + ((int128)tmp_q[6] * 3081232561050645L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 3081232561050645L) + ((int128)tmp_q[1] * 12473233215350655L) + ((int128)tmp_q[2] * 6429671150180222L) - ((int128)tmp_q[3] * 12253194395345975L) - ((-((int128)tmp_q[4] * 2786724442451478L) - ((int128)tmp_q[5] * 7642652496425787L) + ((int128)tmp_q[6] * 14518417682501411L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 14518417682501411L) + ((int128)tmp_q[1] * 3081232561050645L) + ((int128)tmp_q[2] * 12473233215350655L) + ((int128)tmp_q[3] * 6429671150180222L) - ((int128)tmp_q[4] * 12253194395345975L) - ((-((int128)tmp_q[5] * 2786724442451478L) - ((int128)tmp_q[6] * 7642652496425787L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 7642652496425787L) + ((int128)tmp_q[1] * 14518417682501411L) + ((int128)tmp_q[2] * 3081232561050645L) + ((int128)tmp_q[3] * 12473233215350655L) + ((int128)tmp_q[4] * 6429671150180222L) - ((int128)tmp_q[5] * 12253194395345975L) + ((int128)tmp_q[6] * 5573448884902956L);
	tmp_zero[6] = -((int128)tmp_q[0] * 2786724442451478L) - ((int128)tmp_q[1] * 7642652496425787L) + ((int128)tmp_q[2] * 14518417682501411L) + ((int128)tmp_q[3] * 3081232561050645L) + ((int128)tmp_q[4] * 12473233215350655L) + ((int128)tmp_q[5] * 6429671150180222L) - ((int128)tmp_q[6] * 12253194395345975L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

