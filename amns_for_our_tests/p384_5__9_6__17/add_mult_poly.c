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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[8] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13584838056613211991UL) + ((((uint64_t)op[1] * 11716116551999882841UL) + ((uint64_t)op[2] * 12331427776448161388UL) + ((uint64_t)op[3] * 921401310736177333UL) + ((uint64_t)op[4] * 16849785241750233064UL) + ((uint64_t)op[5] * 12237024956014123486UL) + ((uint64_t)op[6] * 4942303750685077511UL) + ((uint64_t)op[7] * 1212380383683010735UL) + ((uint64_t)op[8] * 17491284684779935631UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 17491284684779935631UL) + ((uint64_t)op[1] * 13584838056613211991UL) + ((((uint64_t)op[2] * 11716116551999882841UL) + ((uint64_t)op[3] * 12331427776448161388UL) + ((uint64_t)op[4] * 921401310736177333UL) + ((uint64_t)op[5] * 16849785241750233064UL) + ((uint64_t)op[6] * 12237024956014123486UL) + ((uint64_t)op[7] * 4942303750685077511UL) + ((uint64_t)op[8] * 1212380383683010735UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 1212380383683010735UL) + ((uint64_t)op[1] * 17491284684779935631UL) + ((uint64_t)op[2] * 13584838056613211991UL) + ((((uint64_t)op[3] * 11716116551999882841UL) + ((uint64_t)op[4] * 12331427776448161388UL) + ((uint64_t)op[5] * 921401310736177333UL) + ((uint64_t)op[6] * 16849785241750233064UL) + ((uint64_t)op[7] * 12237024956014123486UL) + ((uint64_t)op[8] * 4942303750685077511UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 4942303750685077511UL) + ((uint64_t)op[1] * 1212380383683010735UL) + ((uint64_t)op[2] * 17491284684779935631UL) + ((uint64_t)op[3] * 13584838056613211991UL) + ((((uint64_t)op[4] * 11716116551999882841UL) + ((uint64_t)op[5] * 12331427776448161388UL) + ((uint64_t)op[6] * 921401310736177333UL) + ((uint64_t)op[7] * 16849785241750233064UL) + ((uint64_t)op[8] * 12237024956014123486UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 12237024956014123486UL) + ((uint64_t)op[1] * 4942303750685077511UL) + ((uint64_t)op[2] * 1212380383683010735UL) + ((uint64_t)op[3] * 17491284684779935631UL) + ((uint64_t)op[4] * 13584838056613211991UL) + ((((uint64_t)op[5] * 11716116551999882841UL) + ((uint64_t)op[6] * 12331427776448161388UL) + ((uint64_t)op[7] * 921401310736177333UL) + ((uint64_t)op[8] * 16849785241750233064UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 16849785241750233064UL) + ((uint64_t)op[1] * 12237024956014123486UL) + ((uint64_t)op[2] * 4942303750685077511UL) + ((uint64_t)op[3] * 1212380383683010735UL) + ((uint64_t)op[4] * 17491284684779935631UL) + ((uint64_t)op[5] * 13584838056613211991UL) + ((((uint64_t)op[6] * 11716116551999882841UL) + ((uint64_t)op[7] * 12331427776448161388UL) + ((uint64_t)op[8] * 921401310736177333UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 921401310736177333UL) + ((uint64_t)op[1] * 16849785241750233064UL) + ((uint64_t)op[2] * 12237024956014123486UL) + ((uint64_t)op[3] * 4942303750685077511UL) + ((uint64_t)op[4] * 1212380383683010735UL) + ((uint64_t)op[5] * 17491284684779935631UL) + ((uint64_t)op[6] * 13584838056613211991UL) + ((((uint64_t)op[7] * 11716116551999882841UL) + ((uint64_t)op[8] * 12331427776448161388UL)) * 6);
	tmp_q[7] = ((uint64_t)op[0] * 12331427776448161388UL) + ((uint64_t)op[1] * 921401310736177333UL) + ((uint64_t)op[2] * 16849785241750233064UL) + ((uint64_t)op[3] * 12237024956014123486UL) + ((uint64_t)op[4] * 4942303750685077511UL) + ((uint64_t)op[5] * 1212380383683010735UL) + ((uint64_t)op[6] * 17491284684779935631UL) + ((uint64_t)op[7] * 13584838056613211991UL) + ((uint64_t)op[8] * 14956467090870642198UL);
	tmp_q[8] = ((uint64_t)op[0] * 11716116551999882841UL) + ((uint64_t)op[1] * 12331427776448161388UL) + ((uint64_t)op[2] * 921401310736177333UL) + ((uint64_t)op[3] * 16849785241750233064UL) + ((uint64_t)op[4] * 12237024956014123486UL) + ((uint64_t)op[5] * 4942303750685077511UL) + ((uint64_t)op[6] * 1212380383683010735UL) + ((uint64_t)op[7] * 17491284684779935631UL) + ((uint64_t)op[8] * 13584838056613211991UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 930027645025L) + ((-((int128)tmp_q[1] * 4307084963511L) + ((int128)tmp_q[2] * 242481408252L) + ((int128)tmp_q[3] * 15980362805L) + ((int128)tmp_q[4] * 3809802179143L) - ((int128)tmp_q[5] * 2921349544185L) + ((int128)tmp_q[6] * 3661298228640L) - ((int128)tmp_q[7] * 2195268924692L) - ((int128)tmp_q[8] * 357516344071L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 357516344071L) - ((int128)tmp_q[1] * 930027645025L) + ((-((int128)tmp_q[2] * 4307084963511L) + ((int128)tmp_q[3] * 242481408252L) + ((int128)tmp_q[4] * 15980362805L) + ((int128)tmp_q[5] * 3809802179143L) - ((int128)tmp_q[6] * 2921349544185L) + ((int128)tmp_q[7] * 3661298228640L) - ((int128)tmp_q[8] * 2195268924692L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 2195268924692L) - ((int128)tmp_q[1] * 357516344071L) - ((int128)tmp_q[2] * 930027645025L) + ((-((int128)tmp_q[3] * 4307084963511L) + ((int128)tmp_q[4] * 242481408252L) + ((int128)tmp_q[5] * 15980362805L) + ((int128)tmp_q[6] * 3809802179143L) - ((int128)tmp_q[7] * 2921349544185L) + ((int128)tmp_q[8] * 3661298228640L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 3661298228640L) - ((int128)tmp_q[1] * 2195268924692L) - ((int128)tmp_q[2] * 357516344071L) - ((int128)tmp_q[3] * 930027645025L) + ((-((int128)tmp_q[4] * 4307084963511L) + ((int128)tmp_q[5] * 242481408252L) + ((int128)tmp_q[6] * 15980362805L) + ((int128)tmp_q[7] * 3809802179143L) - ((int128)tmp_q[8] * 2921349544185L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 2921349544185L) + ((int128)tmp_q[1] * 3661298228640L) - ((int128)tmp_q[2] * 2195268924692L) - ((int128)tmp_q[3] * 357516344071L) - ((int128)tmp_q[4] * 930027645025L) + ((-((int128)tmp_q[5] * 4307084963511L) + ((int128)tmp_q[6] * 242481408252L) + ((int128)tmp_q[7] * 15980362805L) + ((int128)tmp_q[8] * 3809802179143L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 3809802179143L) - ((int128)tmp_q[1] * 2921349544185L) + ((int128)tmp_q[2] * 3661298228640L) - ((int128)tmp_q[3] * 2195268924692L) - ((int128)tmp_q[4] * 357516344071L) - ((int128)tmp_q[5] * 930027645025L) + ((-((int128)tmp_q[6] * 4307084963511L) + ((int128)tmp_q[7] * 242481408252L) + ((int128)tmp_q[8] * 15980362805L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 15980362805L) + ((int128)tmp_q[1] * 3809802179143L) - ((int128)tmp_q[2] * 2921349544185L) + ((int128)tmp_q[3] * 3661298228640L) - ((int128)tmp_q[4] * 2195268924692L) - ((int128)tmp_q[5] * 357516344071L) - ((int128)tmp_q[6] * 930027645025L) + ((-((int128)tmp_q[7] * 4307084963511L) + ((int128)tmp_q[8] * 242481408252L)) * 6);
	tmp_zero[7] = ((int128)tmp_q[0] * 242481408252L) + ((int128)tmp_q[1] * 15980362805L) + ((int128)tmp_q[2] * 3809802179143L) - ((int128)tmp_q[3] * 2921349544185L) + ((int128)tmp_q[4] * 3661298228640L) - ((int128)tmp_q[5] * 2195268924692L) - ((int128)tmp_q[6] * 357516344071L) - ((int128)tmp_q[7] * 930027645025L) - ((int128)tmp_q[8] * 25842509781066L);
	tmp_zero[8] = -((int128)tmp_q[0] * 4307084963511L) + ((int128)tmp_q[1] * 242481408252L) + ((int128)tmp_q[2] * 15980362805L) + ((int128)tmp_q[3] * 3809802179143L) - ((int128)tmp_q[4] * 2921349544185L) + ((int128)tmp_q[5] * 3661298228640L) - ((int128)tmp_q[6] * 2195268924692L) - ((int128)tmp_q[7] * 357516344071L) - ((int128)tmp_q[8] * 930027645025L);

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

