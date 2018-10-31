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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16315988510657202202UL) + ((((uint64_t)op[1] * 6773698739804261198UL) + ((uint64_t)op[2] * 10466725930857149775UL) + ((uint64_t)op[3] * 1229437860093050029UL) + ((uint64_t)op[4] * 7501489525299061949UL) + ((uint64_t)op[5] * 17374028719461085057UL) + ((uint64_t)op[6] * 6259457223450618179UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 6259457223450618179UL) + ((uint64_t)op[1] * 16315988510657202202UL) + ((((uint64_t)op[2] * 6773698739804261198UL) + ((uint64_t)op[3] * 10466725930857149775UL) + ((uint64_t)op[4] * 1229437860093050029UL) + ((uint64_t)op[5] * 7501489525299061949UL) + ((uint64_t)op[6] * 17374028719461085057UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 17374028719461085057UL) + ((uint64_t)op[1] * 6259457223450618179UL) + ((uint64_t)op[2] * 16315988510657202202UL) + ((((uint64_t)op[3] * 6773698739804261198UL) + ((uint64_t)op[4] * 10466725930857149775UL) + ((uint64_t)op[5] * 1229437860093050029UL) + ((uint64_t)op[6] * 7501489525299061949UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7501489525299061949UL) + ((uint64_t)op[1] * 17374028719461085057UL) + ((uint64_t)op[2] * 6259457223450618179UL) + ((uint64_t)op[3] * 16315988510657202202UL) + ((((uint64_t)op[4] * 6773698739804261198UL) + ((uint64_t)op[5] * 10466725930857149775UL) + ((uint64_t)op[6] * 1229437860093050029UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 1229437860093050029UL) + ((uint64_t)op[1] * 7501489525299061949UL) + ((uint64_t)op[2] * 17374028719461085057UL) + ((uint64_t)op[3] * 6259457223450618179UL) + ((uint64_t)op[4] * 16315988510657202202UL) + ((((uint64_t)op[5] * 6773698739804261198UL) + ((uint64_t)op[6] * 10466725930857149775UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 10466725930857149775UL) + ((uint64_t)op[1] * 1229437860093050029UL) + ((uint64_t)op[2] * 7501489525299061949UL) + ((uint64_t)op[3] * 17374028719461085057UL) + ((uint64_t)op[4] * 6259457223450618179UL) + ((uint64_t)op[5] * 16315988510657202202UL) + ((uint64_t)op[6] * 16572391928006319638UL);
	tmp_q[6] = ((uint64_t)op[0] * 6773698739804261198UL) + ((uint64_t)op[1] * 10466725930857149775UL) + ((uint64_t)op[2] * 1229437860093050029UL) + ((uint64_t)op[3] * 7501489525299061949UL) + ((uint64_t)op[4] * 17374028719461085057UL) + ((uint64_t)op[5] * 6259457223450618179UL) + ((uint64_t)op[6] * 16315988510657202202UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 46472808260L) - ((((int128)tmp_q[1] * 52727916693L) + ((int128)tmp_q[2] * 996925026L) + ((int128)tmp_q[3] * 21637903033L) + ((int128)tmp_q[4] * 22328270150L) - ((int128)tmp_q[5] * 49451802855L) + ((int128)tmp_q[6] * 34573781516L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 34573781516L) + ((int128)tmp_q[1] * 46472808260L) - ((((int128)tmp_q[2] * 52727916693L) + ((int128)tmp_q[3] * 996925026L) + ((int128)tmp_q[4] * 21637903033L) + ((int128)tmp_q[5] * 22328270150L) - ((int128)tmp_q[6] * 49451802855L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 49451802855L) + ((int128)tmp_q[1] * 34573781516L) + ((int128)tmp_q[2] * 46472808260L) - ((((int128)tmp_q[3] * 52727916693L) + ((int128)tmp_q[4] * 996925026L) + ((int128)tmp_q[5] * 21637903033L) + ((int128)tmp_q[6] * 22328270150L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 22328270150L) - ((int128)tmp_q[1] * 49451802855L) + ((int128)tmp_q[2] * 34573781516L) + ((int128)tmp_q[3] * 46472808260L) - ((((int128)tmp_q[4] * 52727916693L) + ((int128)tmp_q[5] * 996925026L) + ((int128)tmp_q[6] * 21637903033L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 21637903033L) + ((int128)tmp_q[1] * 22328270150L) - ((int128)tmp_q[2] * 49451802855L) + ((int128)tmp_q[3] * 34573781516L) + ((int128)tmp_q[4] * 46472808260L) - ((((int128)tmp_q[5] * 52727916693L) + ((int128)tmp_q[6] * 996925026L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 996925026L) + ((int128)tmp_q[1] * 21637903033L) + ((int128)tmp_q[2] * 22328270150L) - ((int128)tmp_q[3] * 49451802855L) + ((int128)tmp_q[4] * 34573781516L) + ((int128)tmp_q[5] * 46472808260L) - ((int128)tmp_q[6] * 158183750079L);
	tmp_zero[6] = ((int128)tmp_q[0] * 52727916693L) + ((int128)tmp_q[1] * 996925026L) + ((int128)tmp_q[2] * 21637903033L) + ((int128)tmp_q[3] * 22328270150L) - ((int128)tmp_q[4] * 49451802855L) + ((int128)tmp_q[5] * 34573781516L) + ((int128)tmp_q[6] * 46472808260L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

