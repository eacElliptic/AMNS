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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3002384401159145841UL) + ((((uint64_t)op[1] * 15803187526268403198UL) + ((uint64_t)op[2] * 9379193917681268541UL) + ((uint64_t)op[3] * 14975066608983637906UL) + ((uint64_t)op[4] * 7727533014328160449UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 7727533014328160449UL) + ((uint64_t)op[1] * 3002384401159145841UL) + ((((uint64_t)op[2] * 15803187526268403198UL) + ((uint64_t)op[3] * 9379193917681268541UL) + ((uint64_t)op[4] * 14975066608983637906UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 14975066608983637906UL) + ((uint64_t)op[1] * 7727533014328160449UL) + ((uint64_t)op[2] * 3002384401159145841UL) + ((((uint64_t)op[3] * 15803187526268403198UL) + ((uint64_t)op[4] * 9379193917681268541UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 9379193917681268541UL) + ((uint64_t)op[1] * 14975066608983637906UL) + ((uint64_t)op[2] * 7727533014328160449UL) + ((uint64_t)op[3] * 3002384401159145841UL) + ((uint64_t)op[4] * 7930669642323445254UL);
	tmp_q[4] = ((uint64_t)op[0] * 15803187526268403198UL) + ((uint64_t)op[1] * 9379193917681268541UL) + ((uint64_t)op[2] * 14975066608983637906UL) + ((uint64_t)op[3] * 7727533014328160449UL) + ((uint64_t)op[4] * 3002384401159145841UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2090616818693860L) - ((((int128)tmp_q[1] * 97561858288684L) + ((int128)tmp_q[2] * 898099519903807L) + ((int128)tmp_q[3] * 410208674889303L) + ((int128)tmp_q[4] * 663051365049947L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 663051365049947L) - ((int128)tmp_q[1] * 2090616818693860L) - ((((int128)tmp_q[2] * 97561858288684L) + ((int128)tmp_q[3] * 898099519903807L) + ((int128)tmp_q[4] * 410208674889303L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 410208674889303L) + ((int128)tmp_q[1] * 663051365049947L) - ((int128)tmp_q[2] * 2090616818693860L) - ((((int128)tmp_q[3] * 97561858288684L) + ((int128)tmp_q[4] * 898099519903807L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 898099519903807L) + ((int128)tmp_q[1] * 410208674889303L) + ((int128)tmp_q[2] * 663051365049947L) - ((int128)tmp_q[3] * 2090616818693860L) - ((int128)tmp_q[4] * 292685574866052L);
	tmp_zero[4] = ((int128)tmp_q[0] * 97561858288684L) + ((int128)tmp_q[1] * 898099519903807L) + ((int128)tmp_q[2] * 410208674889303L) + ((int128)tmp_q[3] * 663051365049947L) - ((int128)tmp_q[4] * 2090616818693860L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

