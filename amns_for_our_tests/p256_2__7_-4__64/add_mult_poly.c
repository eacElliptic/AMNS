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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18049762443000244315UL) + ((((uint64_t)op[1] * 12428454108010828594UL) + ((uint64_t)op[2] * 17530306842020031415UL) + ((uint64_t)op[3] * 7880396490724724564UL) + ((uint64_t)op[4] * 4050317895966490404UL) + ((uint64_t)op[5] * 13578291119429057241UL) + ((uint64_t)op[6] * 17320778009269289596UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 17320778009269289596UL) + ((uint64_t)op[1] * 18049762443000244315UL) + ((((uint64_t)op[2] * 12428454108010828594UL) + ((uint64_t)op[3] * 17530306842020031415UL) + ((uint64_t)op[4] * 7880396490724724564UL) + ((uint64_t)op[5] * 4050317895966490404UL) + ((uint64_t)op[6] * 13578291119429057241UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 13578291119429057241UL) + ((uint64_t)op[1] * 17320778009269289596UL) + ((uint64_t)op[2] * 18049762443000244315UL) + ((((uint64_t)op[3] * 12428454108010828594UL) + ((uint64_t)op[4] * 17530306842020031415UL) + ((uint64_t)op[5] * 7880396490724724564UL) + ((uint64_t)op[6] * 4050317895966490404UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 4050317895966490404UL) + ((uint64_t)op[1] * 13578291119429057241UL) + ((uint64_t)op[2] * 17320778009269289596UL) + ((uint64_t)op[3] * 18049762443000244315UL) + ((((uint64_t)op[4] * 12428454108010828594UL) + ((uint64_t)op[5] * 17530306842020031415UL) + ((uint64_t)op[6] * 7880396490724724564UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 7880396490724724564UL) + ((uint64_t)op[1] * 4050317895966490404UL) + ((uint64_t)op[2] * 13578291119429057241UL) + ((uint64_t)op[3] * 17320778009269289596UL) + ((uint64_t)op[4] * 18049762443000244315UL) + ((((uint64_t)op[5] * 12428454108010828594UL) + ((uint64_t)op[6] * 17530306842020031415UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 17530306842020031415UL) + ((uint64_t)op[1] * 7880396490724724564UL) + ((uint64_t)op[2] * 4050317895966490404UL) + ((uint64_t)op[3] * 13578291119429057241UL) + ((uint64_t)op[4] * 17320778009269289596UL) + ((uint64_t)op[5] * 18049762443000244315UL) + ((uint64_t)op[6] * 5626415789085340472UL);
	tmp_q[6] = ((uint64_t)op[0] * 12428454108010828594UL) + ((uint64_t)op[1] * 17530306842020031415UL) + ((uint64_t)op[2] * 7880396490724724564UL) + ((uint64_t)op[3] * 4050317895966490404UL) + ((uint64_t)op[4] * 13578291119429057241UL) + ((uint64_t)op[5] * 17320778009269289596UL) + ((uint64_t)op[6] * 18049762443000244315UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 20969554997L) - ((-((int128)tmp_q[1] * 41221325009L) + ((int128)tmp_q[2] * 24893438955L) + ((int128)tmp_q[3] * 48456215561L) + ((int128)tmp_q[4] * 48909769308L) - ((int128)tmp_q[5] * 1068529539L) + ((int128)tmp_q[6] * 8829100472L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 8829100472L) + ((int128)tmp_q[1] * 20969554997L) - ((-((int128)tmp_q[2] * 41221325009L) + ((int128)tmp_q[3] * 24893438955L) + ((int128)tmp_q[4] * 48456215561L) + ((int128)tmp_q[5] * 48909769308L) - ((int128)tmp_q[6] * 1068529539L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 1068529539L) + ((int128)tmp_q[1] * 8829100472L) + ((int128)tmp_q[2] * 20969554997L) - ((-((int128)tmp_q[3] * 41221325009L) + ((int128)tmp_q[4] * 24893438955L) + ((int128)tmp_q[5] * 48456215561L) + ((int128)tmp_q[6] * 48909769308L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 48909769308L) - ((int128)tmp_q[1] * 1068529539L) + ((int128)tmp_q[2] * 8829100472L) + ((int128)tmp_q[3] * 20969554997L) - ((-((int128)tmp_q[4] * 41221325009L) + ((int128)tmp_q[5] * 24893438955L) + ((int128)tmp_q[6] * 48456215561L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 48456215561L) + ((int128)tmp_q[1] * 48909769308L) - ((int128)tmp_q[2] * 1068529539L) + ((int128)tmp_q[3] * 8829100472L) + ((int128)tmp_q[4] * 20969554997L) - ((-((int128)tmp_q[5] * 41221325009L) + ((int128)tmp_q[6] * 24893438955L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 24893438955L) + ((int128)tmp_q[1] * 48456215561L) + ((int128)tmp_q[2] * 48909769308L) - ((int128)tmp_q[3] * 1068529539L) + ((int128)tmp_q[4] * 8829100472L) + ((int128)tmp_q[5] * 20969554997L) + ((int128)tmp_q[6] * 164885300036L);
	tmp_zero[6] = -((int128)tmp_q[0] * 41221325009L) + ((int128)tmp_q[1] * 24893438955L) + ((int128)tmp_q[2] * 48456215561L) + ((int128)tmp_q[3] * 48909769308L) - ((int128)tmp_q[4] * 1068529539L) + ((int128)tmp_q[5] * 8829100472L) + ((int128)tmp_q[6] * 20969554997L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

