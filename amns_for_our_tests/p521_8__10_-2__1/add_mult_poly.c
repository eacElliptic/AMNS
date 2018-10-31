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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[9] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15647682621430997225UL) + ((((uint64_t)op[1] * 6415123949656556873UL) + ((uint64_t)op[2] * 6710810586351003951UL) + ((uint64_t)op[3] * 3397519781539329897UL) + ((uint64_t)op[4] * 1737439529138179359UL) + ((uint64_t)op[5] * 4358242460326989289UL) + ((uint64_t)op[6] * 18241668298754037588UL) + ((uint64_t)op[7] * 9832664814121185089UL) + ((uint64_t)op[8] * 9801534740758993140UL) + ((uint64_t)op[9] * 5955317286210582802UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 5955317286210582802UL) + ((uint64_t)op[1] * 15647682621430997225UL) + ((((uint64_t)op[2] * 6415123949656556873UL) + ((uint64_t)op[3] * 6710810586351003951UL) + ((uint64_t)op[4] * 3397519781539329897UL) + ((uint64_t)op[5] * 1737439529138179359UL) + ((uint64_t)op[6] * 4358242460326989289UL) + ((uint64_t)op[7] * 18241668298754037588UL) + ((uint64_t)op[8] * 9832664814121185089UL) + ((uint64_t)op[9] * 9801534740758993140UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 9801534740758993140UL) + ((uint64_t)op[1] * 5955317286210582802UL) + ((uint64_t)op[2] * 15647682621430997225UL) + ((((uint64_t)op[3] * 6415123949656556873UL) + ((uint64_t)op[4] * 6710810586351003951UL) + ((uint64_t)op[5] * 3397519781539329897UL) + ((uint64_t)op[6] * 1737439529138179359UL) + ((uint64_t)op[7] * 4358242460326989289UL) + ((uint64_t)op[8] * 18241668298754037588UL) + ((uint64_t)op[9] * 9832664814121185089UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 9832664814121185089UL) + ((uint64_t)op[1] * 9801534740758993140UL) + ((uint64_t)op[2] * 5955317286210582802UL) + ((uint64_t)op[3] * 15647682621430997225UL) + ((((uint64_t)op[4] * 6415123949656556873UL) + ((uint64_t)op[5] * 6710810586351003951UL) + ((uint64_t)op[6] * 3397519781539329897UL) + ((uint64_t)op[7] * 1737439529138179359UL) + ((uint64_t)op[8] * 4358242460326989289UL) + ((uint64_t)op[9] * 18241668298754037588UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 18241668298754037588UL) + ((uint64_t)op[1] * 9832664814121185089UL) + ((uint64_t)op[2] * 9801534740758993140UL) + ((uint64_t)op[3] * 5955317286210582802UL) + ((uint64_t)op[4] * 15647682621430997225UL) + ((((uint64_t)op[5] * 6415123949656556873UL) + ((uint64_t)op[6] * 6710810586351003951UL) + ((uint64_t)op[7] * 3397519781539329897UL) + ((uint64_t)op[8] * 1737439529138179359UL) + ((uint64_t)op[9] * 4358242460326989289UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 4358242460326989289UL) + ((uint64_t)op[1] * 18241668298754037588UL) + ((uint64_t)op[2] * 9832664814121185089UL) + ((uint64_t)op[3] * 9801534740758993140UL) + ((uint64_t)op[4] * 5955317286210582802UL) + ((uint64_t)op[5] * 15647682621430997225UL) + ((((uint64_t)op[6] * 6415123949656556873UL) + ((uint64_t)op[7] * 6710810586351003951UL) + ((uint64_t)op[8] * 3397519781539329897UL) + ((uint64_t)op[9] * 1737439529138179359UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 1737439529138179359UL) + ((uint64_t)op[1] * 4358242460326989289UL) + ((uint64_t)op[2] * 18241668298754037588UL) + ((uint64_t)op[3] * 9832664814121185089UL) + ((uint64_t)op[4] * 9801534740758993140UL) + ((uint64_t)op[5] * 5955317286210582802UL) + ((uint64_t)op[6] * 15647682621430997225UL) + ((((uint64_t)op[7] * 6415123949656556873UL) + ((uint64_t)op[8] * 6710810586351003951UL) + ((uint64_t)op[9] * 3397519781539329897UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 3397519781539329897UL) + ((uint64_t)op[1] * 1737439529138179359UL) + ((uint64_t)op[2] * 4358242460326989289UL) + ((uint64_t)op[3] * 18241668298754037588UL) + ((uint64_t)op[4] * 9832664814121185089UL) + ((uint64_t)op[5] * 9801534740758993140UL) + ((uint64_t)op[6] * 5955317286210582802UL) + ((uint64_t)op[7] * 15647682621430997225UL) + ((((uint64_t)op[8] * 6415123949656556873UL) + ((uint64_t)op[9] * 6710810586351003951UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 6710810586351003951UL) + ((uint64_t)op[1] * 3397519781539329897UL) + ((uint64_t)op[2] * 1737439529138179359UL) + ((uint64_t)op[3] * 4358242460326989289UL) + ((uint64_t)op[4] * 18241668298754037588UL) + ((uint64_t)op[5] * 9832664814121185089UL) + ((uint64_t)op[6] * 9801534740758993140UL) + ((uint64_t)op[7] * 5955317286210582802UL) + ((uint64_t)op[8] * 15647682621430997225UL) + ((uint64_t)op[9] * 5616496174396437870UL);
	tmp_q[9] = ((uint64_t)op[0] * 6415123949656556873UL) + ((uint64_t)op[1] * 6710810586351003951UL) + ((uint64_t)op[2] * 3397519781539329897UL) + ((uint64_t)op[3] * 1737439529138179359UL) + ((uint64_t)op[4] * 4358242460326989289UL) + ((uint64_t)op[5] * 18241668298754037588UL) + ((uint64_t)op[6] * 9832664814121185089UL) + ((uint64_t)op[7] * 9801534740758993140UL) + ((uint64_t)op[8] * 5955317286210582802UL) + ((uint64_t)op[9] * 15647682621430997225UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 824342836518485L) - ((-((int128)tmp_q[1] * 2850347629408298L) + ((int128)tmp_q[2] * 1099898425657175L) + ((int128)tmp_q[3] * 2690538089052247L) - ((int128)tmp_q[4] * 1751218754954890L) + ((int128)tmp_q[5] * 4498638859554481L) - ((int128)tmp_q[6] * 4597861214799460L) - ((int128)tmp_q[7] * 395552995479479L) - ((int128)tmp_q[8] * 822730704834354L) - ((int128)tmp_q[9] * 2269045611191320L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 2269045611191320L) + ((int128)tmp_q[1] * 824342836518485L) - ((-((int128)tmp_q[2] * 2850347629408298L) + ((int128)tmp_q[3] * 1099898425657175L) + ((int128)tmp_q[4] * 2690538089052247L) - ((int128)tmp_q[5] * 1751218754954890L) + ((int128)tmp_q[6] * 4498638859554481L) - ((int128)tmp_q[7] * 4597861214799460L) - ((int128)tmp_q[8] * 395552995479479L) - ((int128)tmp_q[9] * 822730704834354L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 822730704834354L) - ((int128)tmp_q[1] * 2269045611191320L) + ((int128)tmp_q[2] * 824342836518485L) - ((-((int128)tmp_q[3] * 2850347629408298L) + ((int128)tmp_q[4] * 1099898425657175L) + ((int128)tmp_q[5] * 2690538089052247L) - ((int128)tmp_q[6] * 1751218754954890L) + ((int128)tmp_q[7] * 4498638859554481L) - ((int128)tmp_q[8] * 4597861214799460L) - ((int128)tmp_q[9] * 395552995479479L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 395552995479479L) - ((int128)tmp_q[1] * 822730704834354L) - ((int128)tmp_q[2] * 2269045611191320L) + ((int128)tmp_q[3] * 824342836518485L) - ((-((int128)tmp_q[4] * 2850347629408298L) + ((int128)tmp_q[5] * 1099898425657175L) + ((int128)tmp_q[6] * 2690538089052247L) - ((int128)tmp_q[7] * 1751218754954890L) + ((int128)tmp_q[8] * 4498638859554481L) - ((int128)tmp_q[9] * 4597861214799460L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 4597861214799460L) - ((int128)tmp_q[1] * 395552995479479L) - ((int128)tmp_q[2] * 822730704834354L) - ((int128)tmp_q[3] * 2269045611191320L) + ((int128)tmp_q[4] * 824342836518485L) - ((-((int128)tmp_q[5] * 2850347629408298L) + ((int128)tmp_q[6] * 1099898425657175L) + ((int128)tmp_q[7] * 2690538089052247L) - ((int128)tmp_q[8] * 1751218754954890L) + ((int128)tmp_q[9] * 4498638859554481L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 4498638859554481L) - ((int128)tmp_q[1] * 4597861214799460L) - ((int128)tmp_q[2] * 395552995479479L) - ((int128)tmp_q[3] * 822730704834354L) - ((int128)tmp_q[4] * 2269045611191320L) + ((int128)tmp_q[5] * 824342836518485L) - ((-((int128)tmp_q[6] * 2850347629408298L) + ((int128)tmp_q[7] * 1099898425657175L) + ((int128)tmp_q[8] * 2690538089052247L) - ((int128)tmp_q[9] * 1751218754954890L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 1751218754954890L) + ((int128)tmp_q[1] * 4498638859554481L) - ((int128)tmp_q[2] * 4597861214799460L) - ((int128)tmp_q[3] * 395552995479479L) - ((int128)tmp_q[4] * 822730704834354L) - ((int128)tmp_q[5] * 2269045611191320L) + ((int128)tmp_q[6] * 824342836518485L) - ((-((int128)tmp_q[7] * 2850347629408298L) + ((int128)tmp_q[8] * 1099898425657175L) + ((int128)tmp_q[9] * 2690538089052247L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 2690538089052247L) - ((int128)tmp_q[1] * 1751218754954890L) + ((int128)tmp_q[2] * 4498638859554481L) - ((int128)tmp_q[3] * 4597861214799460L) - ((int128)tmp_q[4] * 395552995479479L) - ((int128)tmp_q[5] * 822730704834354L) - ((int128)tmp_q[6] * 2269045611191320L) + ((int128)tmp_q[7] * 824342836518485L) - ((-((int128)tmp_q[8] * 2850347629408298L) + ((int128)tmp_q[9] * 1099898425657175L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 1099898425657175L) + ((int128)tmp_q[1] * 2690538089052247L) - ((int128)tmp_q[2] * 1751218754954890L) + ((int128)tmp_q[3] * 4498638859554481L) - ((int128)tmp_q[4] * 4597861214799460L) - ((int128)tmp_q[5] * 395552995479479L) - ((int128)tmp_q[6] * 822730704834354L) - ((int128)tmp_q[7] * 2269045611191320L) + ((int128)tmp_q[8] * 824342836518485L) + ((int128)tmp_q[9] * 5700695258816596L);
	tmp_zero[9] = -((int128)tmp_q[0] * 2850347629408298L) + ((int128)tmp_q[1] * 1099898425657175L) + ((int128)tmp_q[2] * 2690538089052247L) - ((int128)tmp_q[3] * 1751218754954890L) + ((int128)tmp_q[4] * 4498638859554481L) - ((int128)tmp_q[5] * 4597861214799460L) - ((int128)tmp_q[6] * 395552995479479L) - ((int128)tmp_q[7] * 822730704834354L) - ((int128)tmp_q[8] * 2269045611191320L) + ((int128)tmp_q[9] * 824342836518485L);

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
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
}

