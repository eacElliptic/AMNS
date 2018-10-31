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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3733878089892917969UL) + ((((uint64_t)op[1] * 3795566916124336252UL) + ((uint64_t)op[2] * 1247264899790482450UL) + ((uint64_t)op[3] * 15501723649346778429UL) + ((uint64_t)op[4] * 5576811590702170914UL) + ((uint64_t)op[5] * 11551824129326818103UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 11551824129326818103UL) + ((uint64_t)op[1] * 3733878089892917969UL) + ((((uint64_t)op[2] * 3795566916124336252UL) + ((uint64_t)op[3] * 1247264899790482450UL) + ((uint64_t)op[4] * 15501723649346778429UL) + ((uint64_t)op[5] * 5576811590702170914UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 5576811590702170914UL) + ((uint64_t)op[1] * 11551824129326818103UL) + ((uint64_t)op[2] * 3733878089892917969UL) + ((((uint64_t)op[3] * 3795566916124336252UL) + ((uint64_t)op[4] * 1247264899790482450UL) + ((uint64_t)op[5] * 15501723649346778429UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 15501723649346778429UL) + ((uint64_t)op[1] * 5576811590702170914UL) + ((uint64_t)op[2] * 11551824129326818103UL) + ((uint64_t)op[3] * 3733878089892917969UL) + ((((uint64_t)op[4] * 3795566916124336252UL) + ((uint64_t)op[5] * 1247264899790482450UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 1247264899790482450UL) + ((uint64_t)op[1] * 15501723649346778429UL) + ((uint64_t)op[2] * 5576811590702170914UL) + ((uint64_t)op[3] * 11551824129326818103UL) + ((uint64_t)op[4] * 3733878089892917969UL) + ((uint64_t)op[5] * 11386700748373008756UL);
	tmp_q[5] = ((uint64_t)op[0] * 3795566916124336252UL) + ((uint64_t)op[1] * 1247264899790482450UL) + ((uint64_t)op[2] * 15501723649346778429UL) + ((uint64_t)op[3] * 5576811590702170914UL) + ((uint64_t)op[4] * 11551824129326818103UL) + ((uint64_t)op[5] * 3733878089892917969UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1843716382364L) + ((-((int128)tmp_q[1] * 627880366203L) + ((int128)tmp_q[2] * 1990807641139L) - ((int128)tmp_q[3] * 1626440705436L) - ((int128)tmp_q[4] * 1303445534436L) + ((int128)tmp_q[5] * 4048865897329L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 4048865897329L) + ((int128)tmp_q[1] * 1843716382364L) + ((-((int128)tmp_q[2] * 627880366203L) + ((int128)tmp_q[3] * 1990807641139L) - ((int128)tmp_q[4] * 1626440705436L) - ((int128)tmp_q[5] * 1303445534436L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1303445534436L) + ((int128)tmp_q[1] * 4048865897329L) + ((int128)tmp_q[2] * 1843716382364L) + ((-((int128)tmp_q[3] * 627880366203L) + ((int128)tmp_q[4] * 1990807641139L) - ((int128)tmp_q[5] * 1626440705436L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 1626440705436L) - ((int128)tmp_q[1] * 1303445534436L) + ((int128)tmp_q[2] * 4048865897329L) + ((int128)tmp_q[3] * 1843716382364L) + ((-((int128)tmp_q[4] * 627880366203L) + ((int128)tmp_q[5] * 1990807641139L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1990807641139L) - ((int128)tmp_q[1] * 1626440705436L) - ((int128)tmp_q[2] * 1303445534436L) + ((int128)tmp_q[3] * 4048865897329L) + ((int128)tmp_q[4] * 1843716382364L) - ((int128)tmp_q[5] * 1883641098609L);
	tmp_zero[5] = -((int128)tmp_q[0] * 627880366203L) + ((int128)tmp_q[1] * 1990807641139L) - ((int128)tmp_q[2] * 1626440705436L) - ((int128)tmp_q[3] * 1303445534436L) + ((int128)tmp_q[4] * 4048865897329L) + ((int128)tmp_q[5] * 1843716382364L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

