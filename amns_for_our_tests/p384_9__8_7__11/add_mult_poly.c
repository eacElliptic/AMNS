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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8100713476708046688UL) + ((((uint64_t)op[1] * 14940174662585222709UL) + ((uint64_t)op[2] * 2061517899167447900UL) + ((uint64_t)op[3] * 203515623537422144UL) + ((uint64_t)op[4] * 17168660275846083140UL) + ((uint64_t)op[5] * 8798780635299717374UL) + ((uint64_t)op[6] * 10763982422043019833UL) + ((uint64_t)op[7] * 6464869141138656513UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 6464869141138656513UL) + ((uint64_t)op[1] * 8100713476708046688UL) + ((((uint64_t)op[2] * 14940174662585222709UL) + ((uint64_t)op[3] * 2061517899167447900UL) + ((uint64_t)op[4] * 203515623537422144UL) + ((uint64_t)op[5] * 17168660275846083140UL) + ((uint64_t)op[6] * 8798780635299717374UL) + ((uint64_t)op[7] * 10763982422043019833UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 10763982422043019833UL) + ((uint64_t)op[1] * 6464869141138656513UL) + ((uint64_t)op[2] * 8100713476708046688UL) + ((((uint64_t)op[3] * 14940174662585222709UL) + ((uint64_t)op[4] * 2061517899167447900UL) + ((uint64_t)op[5] * 203515623537422144UL) + ((uint64_t)op[6] * 17168660275846083140UL) + ((uint64_t)op[7] * 8798780635299717374UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 8798780635299717374UL) + ((uint64_t)op[1] * 10763982422043019833UL) + ((uint64_t)op[2] * 6464869141138656513UL) + ((uint64_t)op[3] * 8100713476708046688UL) + ((((uint64_t)op[4] * 14940174662585222709UL) + ((uint64_t)op[5] * 2061517899167447900UL) + ((uint64_t)op[6] * 203515623537422144UL) + ((uint64_t)op[7] * 17168660275846083140UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 17168660275846083140UL) + ((uint64_t)op[1] * 8798780635299717374UL) + ((uint64_t)op[2] * 10763982422043019833UL) + ((uint64_t)op[3] * 6464869141138656513UL) + ((uint64_t)op[4] * 8100713476708046688UL) + ((((uint64_t)op[5] * 14940174662585222709UL) + ((uint64_t)op[6] * 2061517899167447900UL) + ((uint64_t)op[7] * 203515623537422144UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 203515623537422144UL) + ((uint64_t)op[1] * 17168660275846083140UL) + ((uint64_t)op[2] * 8798780635299717374UL) + ((uint64_t)op[3] * 10763982422043019833UL) + ((uint64_t)op[4] * 6464869141138656513UL) + ((uint64_t)op[5] * 8100713476708046688UL) + ((((uint64_t)op[6] * 14940174662585222709UL) + ((uint64_t)op[7] * 2061517899167447900UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 2061517899167447900UL) + ((uint64_t)op[1] * 203515623537422144UL) + ((uint64_t)op[2] * 17168660275846083140UL) + ((uint64_t)op[3] * 8798780635299717374UL) + ((uint64_t)op[4] * 10763982422043019833UL) + ((uint64_t)op[5] * 6464869141138656513UL) + ((uint64_t)op[6] * 8100713476708046688UL) + ((uint64_t)op[7] * 12347502269548800883UL);
	tmp_q[7] = ((uint64_t)op[0] * 14940174662585222709UL) + ((uint64_t)op[1] * 2061517899167447900UL) + ((uint64_t)op[2] * 203515623537422144UL) + ((uint64_t)op[3] * 17168660275846083140UL) + ((uint64_t)op[4] * 8798780635299717374UL) + ((uint64_t)op[5] * 10763982422043019833UL) + ((uint64_t)op[6] * 6464869141138656513UL) + ((uint64_t)op[7] * 8100713476708046688UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 18123098192655L) + ((-((int128)tmp_q[1] * 150010046635941L) - ((int128)tmp_q[2] * 42547259597713L) + ((int128)tmp_q[3] * 42323415066872L) + ((int128)tmp_q[4] * 137246696852475L) - ((int128)tmp_q[5] * 64787152733356L) - ((int128)tmp_q[6] * 46236027879894L) - ((int128)tmp_q[7] * 175065970193577L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 175065970193577L) + ((int128)tmp_q[1] * 18123098192655L) + ((-((int128)tmp_q[2] * 150010046635941L) - ((int128)tmp_q[3] * 42547259597713L) + ((int128)tmp_q[4] * 42323415066872L) + ((int128)tmp_q[5] * 137246696852475L) - ((int128)tmp_q[6] * 64787152733356L) - ((int128)tmp_q[7] * 46236027879894L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 46236027879894L) - ((int128)tmp_q[1] * 175065970193577L) + ((int128)tmp_q[2] * 18123098192655L) + ((-((int128)tmp_q[3] * 150010046635941L) - ((int128)tmp_q[4] * 42547259597713L) + ((int128)tmp_q[5] * 42323415066872L) + ((int128)tmp_q[6] * 137246696852475L) - ((int128)tmp_q[7] * 64787152733356L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 64787152733356L) - ((int128)tmp_q[1] * 46236027879894L) - ((int128)tmp_q[2] * 175065970193577L) + ((int128)tmp_q[3] * 18123098192655L) + ((-((int128)tmp_q[4] * 150010046635941L) - ((int128)tmp_q[5] * 42547259597713L) + ((int128)tmp_q[6] * 42323415066872L) + ((int128)tmp_q[7] * 137246696852475L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 137246696852475L) - ((int128)tmp_q[1] * 64787152733356L) - ((int128)tmp_q[2] * 46236027879894L) - ((int128)tmp_q[3] * 175065970193577L) + ((int128)tmp_q[4] * 18123098192655L) + ((-((int128)tmp_q[5] * 150010046635941L) - ((int128)tmp_q[6] * 42547259597713L) + ((int128)tmp_q[7] * 42323415066872L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 42323415066872L) + ((int128)tmp_q[1] * 137246696852475L) - ((int128)tmp_q[2] * 64787152733356L) - ((int128)tmp_q[3] * 46236027879894L) - ((int128)tmp_q[4] * 175065970193577L) + ((int128)tmp_q[5] * 18123098192655L) + ((-((int128)tmp_q[6] * 150010046635941L) - ((int128)tmp_q[7] * 42547259597713L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 42547259597713L) + ((int128)tmp_q[1] * 42323415066872L) + ((int128)tmp_q[2] * 137246696852475L) - ((int128)tmp_q[3] * 64787152733356L) - ((int128)tmp_q[4] * 46236027879894L) - ((int128)tmp_q[5] * 175065970193577L) + ((int128)tmp_q[6] * 18123098192655L) - ((int128)tmp_q[7] * 1050070326451587L);
	tmp_zero[7] = -((int128)tmp_q[0] * 150010046635941L) - ((int128)tmp_q[1] * 42547259597713L) + ((int128)tmp_q[2] * 42323415066872L) + ((int128)tmp_q[3] * 137246696852475L) - ((int128)tmp_q[4] * 64787152733356L) - ((int128)tmp_q[5] * 46236027879894L) - ((int128)tmp_q[6] * 175065970193577L) + ((int128)tmp_q[7] * 18123098192655L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

