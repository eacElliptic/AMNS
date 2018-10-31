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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[9] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6716707418836620577UL) + ((((uint64_t)op[1] * 6859754077524849273UL) + ((uint64_t)op[2] * 12773765759568315446UL) + ((uint64_t)op[3] * 14130529983452883530UL) + ((uint64_t)op[4] * 16472754894996165140UL) + ((uint64_t)op[5] * 17478959635137798189UL) + ((uint64_t)op[6] * 4264835536011942951UL) + ((uint64_t)op[7] * 1880711957005820535UL) + ((uint64_t)op[8] * 1823513099801667250UL) + ((uint64_t)op[9] * 1091380140581603821UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 1091380140581603821UL) + ((uint64_t)op[1] * 6716707418836620577UL) + ((((uint64_t)op[2] * 6859754077524849273UL) + ((uint64_t)op[3] * 12773765759568315446UL) + ((uint64_t)op[4] * 14130529983452883530UL) + ((uint64_t)op[5] * 16472754894996165140UL) + ((uint64_t)op[6] * 17478959635137798189UL) + ((uint64_t)op[7] * 4264835536011942951UL) + ((uint64_t)op[8] * 1880711957005820535UL) + ((uint64_t)op[9] * 1823513099801667250UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 1823513099801667250UL) + ((uint64_t)op[1] * 1091380140581603821UL) + ((uint64_t)op[2] * 6716707418836620577UL) + ((((uint64_t)op[3] * 6859754077524849273UL) + ((uint64_t)op[4] * 12773765759568315446UL) + ((uint64_t)op[5] * 14130529983452883530UL) + ((uint64_t)op[6] * 16472754894996165140UL) + ((uint64_t)op[7] * 17478959635137798189UL) + ((uint64_t)op[8] * 4264835536011942951UL) + ((uint64_t)op[9] * 1880711957005820535UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 1880711957005820535UL) + ((uint64_t)op[1] * 1823513099801667250UL) + ((uint64_t)op[2] * 1091380140581603821UL) + ((uint64_t)op[3] * 6716707418836620577UL) + ((((uint64_t)op[4] * 6859754077524849273UL) + ((uint64_t)op[5] * 12773765759568315446UL) + ((uint64_t)op[6] * 14130529983452883530UL) + ((uint64_t)op[7] * 16472754894996165140UL) + ((uint64_t)op[8] * 17478959635137798189UL) + ((uint64_t)op[9] * 4264835536011942951UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 4264835536011942951UL) + ((uint64_t)op[1] * 1880711957005820535UL) + ((uint64_t)op[2] * 1823513099801667250UL) + ((uint64_t)op[3] * 1091380140581603821UL) + ((uint64_t)op[4] * 6716707418836620577UL) + ((((uint64_t)op[5] * 6859754077524849273UL) + ((uint64_t)op[6] * 12773765759568315446UL) + ((uint64_t)op[7] * 14130529983452883530UL) + ((uint64_t)op[8] * 16472754894996165140UL) + ((uint64_t)op[9] * 17478959635137798189UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 17478959635137798189UL) + ((uint64_t)op[1] * 4264835536011942951UL) + ((uint64_t)op[2] * 1880711957005820535UL) + ((uint64_t)op[3] * 1823513099801667250UL) + ((uint64_t)op[4] * 1091380140581603821UL) + ((uint64_t)op[5] * 6716707418836620577UL) + ((((uint64_t)op[6] * 6859754077524849273UL) + ((uint64_t)op[7] * 12773765759568315446UL) + ((uint64_t)op[8] * 14130529983452883530UL) + ((uint64_t)op[9] * 16472754894996165140UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 16472754894996165140UL) + ((uint64_t)op[1] * 17478959635137798189UL) + ((uint64_t)op[2] * 4264835536011942951UL) + ((uint64_t)op[3] * 1880711957005820535UL) + ((uint64_t)op[4] * 1823513099801667250UL) + ((uint64_t)op[5] * 1091380140581603821UL) + ((uint64_t)op[6] * 6716707418836620577UL) + ((((uint64_t)op[7] * 6859754077524849273UL) + ((uint64_t)op[8] * 12773765759568315446UL) + ((uint64_t)op[9] * 14130529983452883530UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 14130529983452883530UL) + ((uint64_t)op[1] * 16472754894996165140UL) + ((uint64_t)op[2] * 17478959635137798189UL) + ((uint64_t)op[3] * 4264835536011942951UL) + ((uint64_t)op[4] * 1880711957005820535UL) + ((uint64_t)op[5] * 1823513099801667250UL) + ((uint64_t)op[6] * 1091380140581603821UL) + ((uint64_t)op[7] * 6716707418836620577UL) + ((((uint64_t)op[8] * 6859754077524849273UL) + ((uint64_t)op[9] * 12773765759568315446UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 12773765759568315446UL) + ((uint64_t)op[1] * 14130529983452883530UL) + ((uint64_t)op[2] * 16472754894996165140UL) + ((uint64_t)op[3] * 17478959635137798189UL) + ((uint64_t)op[4] * 4264835536011942951UL) + ((uint64_t)op[5] * 1880711957005820535UL) + ((uint64_t)op[6] * 1823513099801667250UL) + ((uint64_t)op[7] * 1091380140581603821UL) + ((uint64_t)op[8] * 6716707418836620577UL) + ((uint64_t)op[9] * 13719508155049698546UL);
	tmp_q[9] = ((uint64_t)op[0] * 6859754077524849273UL) + ((uint64_t)op[1] * 12773765759568315446UL) + ((uint64_t)op[2] * 14130529983452883530UL) + ((uint64_t)op[3] * 16472754894996165140UL) + ((uint64_t)op[4] * 17478959635137798189UL) + ((uint64_t)op[5] * 4264835536011942951UL) + ((uint64_t)op[6] * 1880711957005820535UL) + ((uint64_t)op[7] * 1823513099801667250UL) + ((uint64_t)op[8] * 1091380140581603821UL) + ((uint64_t)op[9] * 6716707418836620577UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 350214212492787L) + ((-((int128)tmp_q[1] * 2691227340386188L) + ((int128)tmp_q[2] * 714007026433863L) + ((int128)tmp_q[3] * 69388539114458L) - ((int128)tmp_q[4] * 346824066764993L) - ((int128)tmp_q[5] * 15257917953135L) - ((int128)tmp_q[6] * 1654908634011636L) + ((int128)tmp_q[7] * 906441219524482L) + ((int128)tmp_q[8] * 577944121758135L) + ((int128)tmp_q[9] * 1283986867552773L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1283986867552773L) - ((int128)tmp_q[1] * 350214212492787L) + ((-((int128)tmp_q[2] * 2691227340386188L) + ((int128)tmp_q[3] * 714007026433863L) + ((int128)tmp_q[4] * 69388539114458L) - ((int128)tmp_q[5] * 346824066764993L) - ((int128)tmp_q[6] * 15257917953135L) - ((int128)tmp_q[7] * 1654908634011636L) + ((int128)tmp_q[8] * 906441219524482L) + ((int128)tmp_q[9] * 577944121758135L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 577944121758135L) + ((int128)tmp_q[1] * 1283986867552773L) - ((int128)tmp_q[2] * 350214212492787L) + ((-((int128)tmp_q[3] * 2691227340386188L) + ((int128)tmp_q[4] * 714007026433863L) + ((int128)tmp_q[5] * 69388539114458L) - ((int128)tmp_q[6] * 346824066764993L) - ((int128)tmp_q[7] * 15257917953135L) - ((int128)tmp_q[8] * 1654908634011636L) + ((int128)tmp_q[9] * 906441219524482L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 906441219524482L) + ((int128)tmp_q[1] * 577944121758135L) + ((int128)tmp_q[2] * 1283986867552773L) - ((int128)tmp_q[3] * 350214212492787L) + ((-((int128)tmp_q[4] * 2691227340386188L) + ((int128)tmp_q[5] * 714007026433863L) + ((int128)tmp_q[6] * 69388539114458L) - ((int128)tmp_q[7] * 346824066764993L) - ((int128)tmp_q[8] * 15257917953135L) - ((int128)tmp_q[9] * 1654908634011636L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1654908634011636L) + ((int128)tmp_q[1] * 906441219524482L) + ((int128)tmp_q[2] * 577944121758135L) + ((int128)tmp_q[3] * 1283986867552773L) - ((int128)tmp_q[4] * 350214212492787L) + ((-((int128)tmp_q[5] * 2691227340386188L) + ((int128)tmp_q[6] * 714007026433863L) + ((int128)tmp_q[7] * 69388539114458L) - ((int128)tmp_q[8] * 346824066764993L) - ((int128)tmp_q[9] * 15257917953135L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 15257917953135L) - ((int128)tmp_q[1] * 1654908634011636L) + ((int128)tmp_q[2] * 906441219524482L) + ((int128)tmp_q[3] * 577944121758135L) + ((int128)tmp_q[4] * 1283986867552773L) - ((int128)tmp_q[5] * 350214212492787L) + ((-((int128)tmp_q[6] * 2691227340386188L) + ((int128)tmp_q[7] * 714007026433863L) + ((int128)tmp_q[8] * 69388539114458L) - ((int128)tmp_q[9] * 346824066764993L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 346824066764993L) - ((int128)tmp_q[1] * 15257917953135L) - ((int128)tmp_q[2] * 1654908634011636L) + ((int128)tmp_q[3] * 906441219524482L) + ((int128)tmp_q[4] * 577944121758135L) + ((int128)tmp_q[5] * 1283986867552773L) - ((int128)tmp_q[6] * 350214212492787L) + ((-((int128)tmp_q[7] * 2691227340386188L) + ((int128)tmp_q[8] * 714007026433863L) + ((int128)tmp_q[9] * 69388539114458L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 69388539114458L) - ((int128)tmp_q[1] * 346824066764993L) - ((int128)tmp_q[2] * 15257917953135L) - ((int128)tmp_q[3] * 1654908634011636L) + ((int128)tmp_q[4] * 906441219524482L) + ((int128)tmp_q[5] * 577944121758135L) + ((int128)tmp_q[6] * 1283986867552773L) - ((int128)tmp_q[7] * 350214212492787L) + ((-((int128)tmp_q[8] * 2691227340386188L) + ((int128)tmp_q[9] * 714007026433863L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 714007026433863L) + ((int128)tmp_q[1] * 69388539114458L) - ((int128)tmp_q[2] * 346824066764993L) - ((int128)tmp_q[3] * 15257917953135L) - ((int128)tmp_q[4] * 1654908634011636L) + ((int128)tmp_q[5] * 906441219524482L) + ((int128)tmp_q[6] * 577944121758135L) + ((int128)tmp_q[7] * 1283986867552773L) - ((int128)tmp_q[8] * 350214212492787L) - ((int128)tmp_q[9] * 5382454680772376L);
	tmp_zero[9] = -((int128)tmp_q[0] * 2691227340386188L) + ((int128)tmp_q[1] * 714007026433863L) + ((int128)tmp_q[2] * 69388539114458L) - ((int128)tmp_q[3] * 346824066764993L) - ((int128)tmp_q[4] * 15257917953135L) - ((int128)tmp_q[5] * 1654908634011636L) + ((int128)tmp_q[6] * 906441219524482L) + ((int128)tmp_q[7] * 577944121758135L) + ((int128)tmp_q[8] * 1283986867552773L) - ((int128)tmp_q[9] * 350214212492787L);

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

