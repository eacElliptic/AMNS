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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16934668233843138616UL) + ((((uint64_t)op[1] * 2331170766855607400UL) + ((uint64_t)op[2] * 13127639771389117299UL) + ((uint64_t)op[3] * 8374529937500559224UL) + ((uint64_t)op[4] * 2358062124876237051UL) + ((uint64_t)op[5] * 854300810887776094UL) + ((uint64_t)op[6] * 7201865685962875529UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 7201865685962875529UL) + ((uint64_t)op[1] * 16934668233843138616UL) + ((((uint64_t)op[2] * 2331170766855607400UL) + ((uint64_t)op[3] * 13127639771389117299UL) + ((uint64_t)op[4] * 8374529937500559224UL) + ((uint64_t)op[5] * 2358062124876237051UL) + ((uint64_t)op[6] * 854300810887776094UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 854300810887776094UL) + ((uint64_t)op[1] * 7201865685962875529UL) + ((uint64_t)op[2] * 16934668233843138616UL) + ((((uint64_t)op[3] * 2331170766855607400UL) + ((uint64_t)op[4] * 13127639771389117299UL) + ((uint64_t)op[5] * 8374529937500559224UL) + ((uint64_t)op[6] * 2358062124876237051UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 2358062124876237051UL) + ((uint64_t)op[1] * 854300810887776094UL) + ((uint64_t)op[2] * 7201865685962875529UL) + ((uint64_t)op[3] * 16934668233843138616UL) + ((((uint64_t)op[4] * 2331170766855607400UL) + ((uint64_t)op[5] * 13127639771389117299UL) + ((uint64_t)op[6] * 8374529937500559224UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 8374529937500559224UL) + ((uint64_t)op[1] * 2358062124876237051UL) + ((uint64_t)op[2] * 854300810887776094UL) + ((uint64_t)op[3] * 7201865685962875529UL) + ((uint64_t)op[4] * 16934668233843138616UL) + ((((uint64_t)op[5] * 2331170766855607400UL) + ((uint64_t)op[6] * 13127639771389117299UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 13127639771389117299UL) + ((uint64_t)op[1] * 8374529937500559224UL) + ((uint64_t)op[2] * 2358062124876237051UL) + ((uint64_t)op[3] * 854300810887776094UL) + ((uint64_t)op[4] * 7201865685962875529UL) + ((uint64_t)op[5] * 16934668233843138616UL) + ((uint64_t)op[6] * 6993512300566822200UL);
	tmp_q[6] = ((uint64_t)op[0] * 2331170766855607400UL) + ((uint64_t)op[1] * 13127639771389117299UL) + ((uint64_t)op[2] * 8374529937500559224UL) + ((uint64_t)op[3] * 2358062124876237051UL) + ((uint64_t)op[4] * 854300810887776094UL) + ((uint64_t)op[5] * 7201865685962875529UL) + ((uint64_t)op[6] * 16934668233843138616UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 24792832462L) + ((((int128)tmp_q[1] * 4567494305L) - ((int128)tmp_q[2] * 12008407257L) + ((int128)tmp_q[3] * 25666663763L) + ((int128)tmp_q[4] * 53062800085L) + ((int128)tmp_q[5] * 17697206059L) - ((int128)tmp_q[6] * 70045970848L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 70045970848L) - ((int128)tmp_q[1] * 24792832462L) + ((((int128)tmp_q[2] * 4567494305L) - ((int128)tmp_q[3] * 12008407257L) + ((int128)tmp_q[4] * 25666663763L) + ((int128)tmp_q[5] * 53062800085L) + ((int128)tmp_q[6] * 17697206059L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 17697206059L) - ((int128)tmp_q[1] * 70045970848L) - ((int128)tmp_q[2] * 24792832462L) + ((((int128)tmp_q[3] * 4567494305L) - ((int128)tmp_q[4] * 12008407257L) + ((int128)tmp_q[5] * 25666663763L) + ((int128)tmp_q[6] * 53062800085L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 53062800085L) + ((int128)tmp_q[1] * 17697206059L) - ((int128)tmp_q[2] * 70045970848L) - ((int128)tmp_q[3] * 24792832462L) + ((((int128)tmp_q[4] * 4567494305L) - ((int128)tmp_q[5] * 12008407257L) + ((int128)tmp_q[6] * 25666663763L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 25666663763L) + ((int128)tmp_q[1] * 53062800085L) + ((int128)tmp_q[2] * 17697206059L) - ((int128)tmp_q[3] * 70045970848L) - ((int128)tmp_q[4] * 24792832462L) + ((((int128)tmp_q[5] * 4567494305L) - ((int128)tmp_q[6] * 12008407257L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 12008407257L) + ((int128)tmp_q[1] * 25666663763L) + ((int128)tmp_q[2] * 53062800085L) + ((int128)tmp_q[3] * 17697206059L) - ((int128)tmp_q[4] * 70045970848L) - ((int128)tmp_q[5] * 24792832462L) + ((int128)tmp_q[6] * 13702482915L);
	tmp_zero[6] = ((int128)tmp_q[0] * 4567494305L) - ((int128)tmp_q[1] * 12008407257L) + ((int128)tmp_q[2] * 25666663763L) + ((int128)tmp_q[3] * 53062800085L) + ((int128)tmp_q[4] * 17697206059L) - ((int128)tmp_q[5] * 70045970848L) - ((int128)tmp_q[6] * 24792832462L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

