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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10339561890339207591UL) + ((((uint64_t)op[1] * 10775152400154290475UL) + ((uint64_t)op[2] * 547712738508328561UL) + ((uint64_t)op[3] * 1351787146618180400UL) + ((uint64_t)op[4] * 7516400335553397953UL) + ((uint64_t)op[5] * 304134797001482277UL) + ((uint64_t)op[6] * 14499674114735526096UL) + ((uint64_t)op[7] * 17293731447877993510UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 17293731447877993510UL) + ((uint64_t)op[1] * 10339561890339207591UL) + ((((uint64_t)op[2] * 10775152400154290475UL) + ((uint64_t)op[3] * 547712738508328561UL) + ((uint64_t)op[4] * 1351787146618180400UL) + ((uint64_t)op[5] * 7516400335553397953UL) + ((uint64_t)op[6] * 304134797001482277UL) + ((uint64_t)op[7] * 14499674114735526096UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 14499674114735526096UL) + ((uint64_t)op[1] * 17293731447877993510UL) + ((uint64_t)op[2] * 10339561890339207591UL) + ((((uint64_t)op[3] * 10775152400154290475UL) + ((uint64_t)op[4] * 547712738508328561UL) + ((uint64_t)op[5] * 1351787146618180400UL) + ((uint64_t)op[6] * 7516400335553397953UL) + ((uint64_t)op[7] * 304134797001482277UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 304134797001482277UL) + ((uint64_t)op[1] * 14499674114735526096UL) + ((uint64_t)op[2] * 17293731447877993510UL) + ((uint64_t)op[3] * 10339561890339207591UL) + ((((uint64_t)op[4] * 10775152400154290475UL) + ((uint64_t)op[5] * 547712738508328561UL) + ((uint64_t)op[6] * 1351787146618180400UL) + ((uint64_t)op[7] * 7516400335553397953UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 7516400335553397953UL) + ((uint64_t)op[1] * 304134797001482277UL) + ((uint64_t)op[2] * 14499674114735526096UL) + ((uint64_t)op[3] * 17293731447877993510UL) + ((uint64_t)op[4] * 10339561890339207591UL) + ((((uint64_t)op[5] * 10775152400154290475UL) + ((uint64_t)op[6] * 547712738508328561UL) + ((uint64_t)op[7] * 1351787146618180400UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 1351787146618180400UL) + ((uint64_t)op[1] * 7516400335553397953UL) + ((uint64_t)op[2] * 304134797001482277UL) + ((uint64_t)op[3] * 14499674114735526096UL) + ((uint64_t)op[4] * 17293731447877993510UL) + ((uint64_t)op[5] * 10339561890339207591UL) + ((((uint64_t)op[6] * 10775152400154290475UL) + ((uint64_t)op[7] * 547712738508328561UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 547712738508328561UL) + ((uint64_t)op[1] * 1351787146618180400UL) + ((uint64_t)op[2] * 7516400335553397953UL) + ((uint64_t)op[3] * 304134797001482277UL) + ((uint64_t)op[4] * 14499674114735526096UL) + ((uint64_t)op[5] * 17293731447877993510UL) + ((uint64_t)op[6] * 10339561890339207591UL) + ((uint64_t)op[7] * 16982273853352349143UL);
	tmp_q[7] = ((uint64_t)op[0] * 10775152400154290475UL) + ((uint64_t)op[1] * 547712738508328561UL) + ((uint64_t)op[2] * 1351787146618180400UL) + ((uint64_t)op[3] * 7516400335553397953UL) + ((uint64_t)op[4] * 304134797001482277UL) + ((uint64_t)op[5] * 14499674114735526096UL) + ((uint64_t)op[6] * 17293731447877993510UL) + ((uint64_t)op[7] * 10339561890339207591UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 163465289960341L) + ((((int128)tmp_q[1] * 56173068320253L) + ((int128)tmp_q[2] * 26501689764026L) + ((int128)tmp_q[3] * 15036676722530L) + ((int128)tmp_q[4] * 136436369817801L) - ((int128)tmp_q[5] * 80798381729161L) - ((int128)tmp_q[6] * 196206109786651L) + ((int128)tmp_q[7] * 131238450327244L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 131238450327244L) - ((int128)tmp_q[1] * 163465289960341L) + ((((int128)tmp_q[2] * 56173068320253L) + ((int128)tmp_q[3] * 26501689764026L) + ((int128)tmp_q[4] * 15036676722530L) + ((int128)tmp_q[5] * 136436369817801L) - ((int128)tmp_q[6] * 80798381729161L) - ((int128)tmp_q[7] * 196206109786651L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 196206109786651L) + ((int128)tmp_q[1] * 131238450327244L) - ((int128)tmp_q[2] * 163465289960341L) + ((((int128)tmp_q[3] * 56173068320253L) + ((int128)tmp_q[4] * 26501689764026L) + ((int128)tmp_q[5] * 15036676722530L) + ((int128)tmp_q[6] * 136436369817801L) - ((int128)tmp_q[7] * 80798381729161L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 80798381729161L) - ((int128)tmp_q[1] * 196206109786651L) + ((int128)tmp_q[2] * 131238450327244L) - ((int128)tmp_q[3] * 163465289960341L) + ((((int128)tmp_q[4] * 56173068320253L) + ((int128)tmp_q[5] * 26501689764026L) + ((int128)tmp_q[6] * 15036676722530L) + ((int128)tmp_q[7] * 136436369817801L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 136436369817801L) - ((int128)tmp_q[1] * 80798381729161L) - ((int128)tmp_q[2] * 196206109786651L) + ((int128)tmp_q[3] * 131238450327244L) - ((int128)tmp_q[4] * 163465289960341L) + ((((int128)tmp_q[5] * 56173068320253L) + ((int128)tmp_q[6] * 26501689764026L) + ((int128)tmp_q[7] * 15036676722530L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 15036676722530L) + ((int128)tmp_q[1] * 136436369817801L) - ((int128)tmp_q[2] * 80798381729161L) - ((int128)tmp_q[3] * 196206109786651L) + ((int128)tmp_q[4] * 131238450327244L) - ((int128)tmp_q[5] * 163465289960341L) + ((((int128)tmp_q[6] * 56173068320253L) + ((int128)tmp_q[7] * 26501689764026L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 26501689764026L) + ((int128)tmp_q[1] * 15036676722530L) + ((int128)tmp_q[2] * 136436369817801L) - ((int128)tmp_q[3] * 80798381729161L) - ((int128)tmp_q[4] * 196206109786651L) + ((int128)tmp_q[5] * 131238450327244L) - ((int128)tmp_q[6] * 163465289960341L) + ((int128)tmp_q[7] * 280865341601265L);
	tmp_zero[7] = ((int128)tmp_q[0] * 56173068320253L) + ((int128)tmp_q[1] * 26501689764026L) + ((int128)tmp_q[2] * 15036676722530L) + ((int128)tmp_q[3] * 136436369817801L) - ((int128)tmp_q[4] * 80798381729161L) - ((int128)tmp_q[5] * 196206109786651L) + ((int128)tmp_q[6] * 131238450327244L) - ((int128)tmp_q[7] * 163465289960341L);

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

