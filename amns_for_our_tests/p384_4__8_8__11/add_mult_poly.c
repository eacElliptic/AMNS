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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 4);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7242471128902094845UL) + ((((uint64_t)op[1] * 7644172920721650380UL) + ((uint64_t)op[2] * 14480972243987255592UL) + ((uint64_t)op[3] * 16440940025626422747UL) + ((uint64_t)op[4] * 13206231148641176531UL) + ((uint64_t)op[5] * 2050967588468694256UL) + ((uint64_t)op[6] * 10885218272097998207UL) + ((uint64_t)op[7] * 5404545498586556099UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 5404545498586556099UL) + ((uint64_t)op[1] * 7242471128902094845UL) + ((((uint64_t)op[2] * 7644172920721650380UL) + ((uint64_t)op[3] * 14480972243987255592UL) + ((uint64_t)op[4] * 16440940025626422747UL) + ((uint64_t)op[5] * 13206231148641176531UL) + ((uint64_t)op[6] * 2050967588468694256UL) + ((uint64_t)op[7] * 10885218272097998207UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 10885218272097998207UL) + ((uint64_t)op[1] * 5404545498586556099UL) + ((uint64_t)op[2] * 7242471128902094845UL) + ((((uint64_t)op[3] * 7644172920721650380UL) + ((uint64_t)op[4] * 14480972243987255592UL) + ((uint64_t)op[5] * 16440940025626422747UL) + ((uint64_t)op[6] * 13206231148641176531UL) + ((uint64_t)op[7] * 2050967588468694256UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 2050967588468694256UL) + ((uint64_t)op[1] * 10885218272097998207UL) + ((uint64_t)op[2] * 5404545498586556099UL) + ((uint64_t)op[3] * 7242471128902094845UL) + ((((uint64_t)op[4] * 7644172920721650380UL) + ((uint64_t)op[5] * 14480972243987255592UL) + ((uint64_t)op[6] * 16440940025626422747UL) + ((uint64_t)op[7] * 13206231148641176531UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 13206231148641176531UL) + ((uint64_t)op[1] * 2050967588468694256UL) + ((uint64_t)op[2] * 10885218272097998207UL) + ((uint64_t)op[3] * 5404545498586556099UL) + ((uint64_t)op[4] * 7242471128902094845UL) + ((((uint64_t)op[5] * 7644172920721650380UL) + ((uint64_t)op[6] * 14480972243987255592UL) + ((uint64_t)op[7] * 16440940025626422747UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 16440940025626422747UL) + ((uint64_t)op[1] * 13206231148641176531UL) + ((uint64_t)op[2] * 2050967588468694256UL) + ((uint64_t)op[3] * 10885218272097998207UL) + ((uint64_t)op[4] * 5404545498586556099UL) + ((uint64_t)op[5] * 7242471128902094845UL) + ((((uint64_t)op[6] * 7644172920721650380UL) + ((uint64_t)op[7] * 14480972243987255592UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 14480972243987255592UL) + ((uint64_t)op[1] * 16440940025626422747UL) + ((uint64_t)op[2] * 13206231148641176531UL) + ((uint64_t)op[3] * 2050967588468694256UL) + ((uint64_t)op[4] * 10885218272097998207UL) + ((uint64_t)op[5] * 5404545498586556099UL) + ((uint64_t)op[6] * 7242471128902094845UL) + ((uint64_t)op[7] * 5813151144644548192UL);
	tmp_q[7] = ((uint64_t)op[0] * 7644172920721650380UL) + ((uint64_t)op[1] * 14480972243987255592UL) + ((uint64_t)op[2] * 16440940025626422747UL) + ((uint64_t)op[3] * 13206231148641176531UL) + ((uint64_t)op[4] * 2050967588468694256UL) + ((uint64_t)op[5] * 10885218272097998207UL) + ((uint64_t)op[6] * 5404545498586556099UL) + ((uint64_t)op[7] * 7242471128902094845UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 115477741724829L) + ((-((int128)tmp_q[1] * 14601016889216L) - ((int128)tmp_q[2] * 135395430639564L) + ((int128)tmp_q[3] * 90984076226041L) + ((int128)tmp_q[4] * 48081346443830L) - ((int128)tmp_q[5] * 6768647266423L) - ((int128)tmp_q[6] * 71710245822206L) + ((int128)tmp_q[7] * 88845416132251L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 88845416132251L) - ((int128)tmp_q[1] * 115477741724829L) + ((-((int128)tmp_q[2] * 14601016889216L) - ((int128)tmp_q[3] * 135395430639564L) + ((int128)tmp_q[4] * 90984076226041L) + ((int128)tmp_q[5] * 48081346443830L) - ((int128)tmp_q[6] * 6768647266423L) - ((int128)tmp_q[7] * 71710245822206L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 71710245822206L) + ((int128)tmp_q[1] * 88845416132251L) - ((int128)tmp_q[2] * 115477741724829L) + ((-((int128)tmp_q[3] * 14601016889216L) - ((int128)tmp_q[4] * 135395430639564L) + ((int128)tmp_q[5] * 90984076226041L) + ((int128)tmp_q[6] * 48081346443830L) - ((int128)tmp_q[7] * 6768647266423L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 6768647266423L) - ((int128)tmp_q[1] * 71710245822206L) + ((int128)tmp_q[2] * 88845416132251L) - ((int128)tmp_q[3] * 115477741724829L) + ((-((int128)tmp_q[4] * 14601016889216L) - ((int128)tmp_q[5] * 135395430639564L) + ((int128)tmp_q[6] * 90984076226041L) + ((int128)tmp_q[7] * 48081346443830L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 48081346443830L) - ((int128)tmp_q[1] * 6768647266423L) - ((int128)tmp_q[2] * 71710245822206L) + ((int128)tmp_q[3] * 88845416132251L) - ((int128)tmp_q[4] * 115477741724829L) + ((-((int128)tmp_q[5] * 14601016889216L) - ((int128)tmp_q[6] * 135395430639564L) + ((int128)tmp_q[7] * 90984076226041L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 90984076226041L) + ((int128)tmp_q[1] * 48081346443830L) - ((int128)tmp_q[2] * 6768647266423L) - ((int128)tmp_q[3] * 71710245822206L) + ((int128)tmp_q[4] * 88845416132251L) - ((int128)tmp_q[5] * 115477741724829L) + ((-((int128)tmp_q[6] * 14601016889216L) - ((int128)tmp_q[7] * 135395430639564L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 135395430639564L) + ((int128)tmp_q[1] * 90984076226041L) + ((int128)tmp_q[2] * 48081346443830L) - ((int128)tmp_q[3] * 6768647266423L) - ((int128)tmp_q[4] * 71710245822206L) + ((int128)tmp_q[5] * 88845416132251L) - ((int128)tmp_q[6] * 115477741724829L) - ((int128)tmp_q[7] * 116808135113728L);
	tmp_zero[7] = -((int128)tmp_q[0] * 14601016889216L) - ((int128)tmp_q[1] * 135395430639564L) + ((int128)tmp_q[2] * 90984076226041L) + ((int128)tmp_q[3] * 48081346443830L) - ((int128)tmp_q[4] * 6768647266423L) - ((int128)tmp_q[5] * 71710245822206L) + ((int128)tmp_q[6] * 88845416132251L) - ((int128)tmp_q[7] * 115477741724829L);

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

