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
	tmp_q[0] = ((uint64_t)op[0] * 15330133497183739967UL) + ((((uint64_t)op[1] * 11232131588037336444UL) + ((uint64_t)op[2] * 1873456372851225408UL) + ((uint64_t)op[3] * 11428430107267506624UL) + ((uint64_t)op[4] * 3465389611299559872UL) + ((uint64_t)op[5] * 14098255663379004048UL) + ((uint64_t)op[6] * 13383517588130894646UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 13383517588130894646UL) + ((uint64_t)op[1] * 15330133497183739967UL) + ((((uint64_t)op[2] * 11232131588037336444UL) + ((uint64_t)op[3] * 1873456372851225408UL) + ((uint64_t)op[4] * 11428430107267506624UL) + ((uint64_t)op[5] * 3465389611299559872UL) + ((uint64_t)op[6] * 14098255663379004048UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 14098255663379004048UL) + ((uint64_t)op[1] * 13383517588130894646UL) + ((uint64_t)op[2] * 15330133497183739967UL) + ((((uint64_t)op[3] * 11232131588037336444UL) + ((uint64_t)op[4] * 1873456372851225408UL) + ((uint64_t)op[5] * 11428430107267506624UL) + ((uint64_t)op[6] * 3465389611299559872UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 3465389611299559872UL) + ((uint64_t)op[1] * 14098255663379004048UL) + ((uint64_t)op[2] * 13383517588130894646UL) + ((uint64_t)op[3] * 15330133497183739967UL) + ((((uint64_t)op[4] * 11232131588037336444UL) + ((uint64_t)op[5] * 1873456372851225408UL) + ((uint64_t)op[6] * 11428430107267506624UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 11428430107267506624UL) + ((uint64_t)op[1] * 3465389611299559872UL) + ((uint64_t)op[2] * 14098255663379004048UL) + ((uint64_t)op[3] * 13383517588130894646UL) + ((uint64_t)op[4] * 15330133497183739967UL) + ((((uint64_t)op[5] * 11232131588037336444UL) + ((uint64_t)op[6] * 1873456372851225408UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 1873456372851225408UL) + ((uint64_t)op[1] * 11428430107267506624UL) + ((uint64_t)op[2] * 3465389611299559872UL) + ((uint64_t)op[3] * 14098255663379004048UL) + ((uint64_t)op[4] * 13383517588130894646UL) + ((uint64_t)op[5] * 15330133497183739967UL) + ((uint64_t)op[6] * 16070076409460485088UL);
	tmp_q[6] = ((uint64_t)op[0] * 11232131588037336444UL) + ((uint64_t)op[1] * 1873456372851225408UL) + ((uint64_t)op[2] * 11428430107267506624UL) + ((uint64_t)op[3] * 3465389611299559872UL) + ((uint64_t)op[4] * 14098255663379004048UL) + ((uint64_t)op[5] * 13383517588130894646UL) + ((uint64_t)op[6] * 15330133497183739967UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 91616327361L) + ((-((int128)tmp_q[1] * 75076034628L) + ((int128)tmp_q[2] * 9124750112L) - ((int128)tmp_q[3] * 57774699632L) - ((int128)tmp_q[4] * 18478659688L) + ((int128)tmp_q[5] * 11484948724L) - ((int128)tmp_q[6] * 14050902858L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 14050902858L) + ((int128)tmp_q[1] * 91616327361L) + ((-((int128)tmp_q[2] * 75076034628L) + ((int128)tmp_q[3] * 9124750112L) - ((int128)tmp_q[4] * 57774699632L) - ((int128)tmp_q[5] * 18478659688L) + ((int128)tmp_q[6] * 11484948724L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 11484948724L) - ((int128)tmp_q[1] * 14050902858L) + ((int128)tmp_q[2] * 91616327361L) + ((-((int128)tmp_q[3] * 75076034628L) + ((int128)tmp_q[4] * 9124750112L) - ((int128)tmp_q[5] * 57774699632L) - ((int128)tmp_q[6] * 18478659688L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 18478659688L) + ((int128)tmp_q[1] * 11484948724L) - ((int128)tmp_q[2] * 14050902858L) + ((int128)tmp_q[3] * 91616327361L) + ((-((int128)tmp_q[4] * 75076034628L) + ((int128)tmp_q[5] * 9124750112L) - ((int128)tmp_q[6] * 57774699632L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 57774699632L) - ((int128)tmp_q[1] * 18478659688L) + ((int128)tmp_q[2] * 11484948724L) - ((int128)tmp_q[3] * 14050902858L) + ((int128)tmp_q[4] * 91616327361L) + ((-((int128)tmp_q[5] * 75076034628L) + ((int128)tmp_q[6] * 9124750112L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 9124750112L) - ((int128)tmp_q[1] * 57774699632L) - ((int128)tmp_q[2] * 18478659688L) + ((int128)tmp_q[3] * 11484948724L) - ((int128)tmp_q[4] * 14050902858L) + ((int128)tmp_q[5] * 91616327361L) - ((int128)tmp_q[6] * 600608277024L);
	tmp_zero[6] = -((int128)tmp_q[0] * 75076034628L) + ((int128)tmp_q[1] * 9124750112L) - ((int128)tmp_q[2] * 57774699632L) - ((int128)tmp_q[3] * 18478659688L) + ((int128)tmp_q[4] * 11484948724L) - ((int128)tmp_q[5] * 14050902858L) + ((int128)tmp_q[6] * 91616327361L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

