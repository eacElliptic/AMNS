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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 605019360372903003UL) + ((((uint64_t)op[1] * 1584614376596801529UL) + ((uint64_t)op[2] * 13002154218893442434UL) + ((uint64_t)op[3] * 5782342310597378278UL) + ((uint64_t)op[4] * 997445221193231493UL) + ((uint64_t)op[5] * 13375432650701855941UL) + ((uint64_t)op[6] * 1705753999475556693UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 1705753999475556693UL) + ((uint64_t)op[1] * 605019360372903003UL) + ((((uint64_t)op[2] * 1584614376596801529UL) + ((uint64_t)op[3] * 13002154218893442434UL) + ((uint64_t)op[4] * 5782342310597378278UL) + ((uint64_t)op[5] * 997445221193231493UL) + ((uint64_t)op[6] * 13375432650701855941UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13375432650701855941UL) + ((uint64_t)op[1] * 1705753999475556693UL) + ((uint64_t)op[2] * 605019360372903003UL) + ((((uint64_t)op[3] * 1584614376596801529UL) + ((uint64_t)op[4] * 13002154218893442434UL) + ((uint64_t)op[5] * 5782342310597378278UL) + ((uint64_t)op[6] * 997445221193231493UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 997445221193231493UL) + ((uint64_t)op[1] * 13375432650701855941UL) + ((uint64_t)op[2] * 1705753999475556693UL) + ((uint64_t)op[3] * 605019360372903003UL) + ((((uint64_t)op[4] * 1584614376596801529UL) + ((uint64_t)op[5] * 13002154218893442434UL) + ((uint64_t)op[6] * 5782342310597378278UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 5782342310597378278UL) + ((uint64_t)op[1] * 997445221193231493UL) + ((uint64_t)op[2] * 13375432650701855941UL) + ((uint64_t)op[3] * 1705753999475556693UL) + ((uint64_t)op[4] * 605019360372903003UL) + ((((uint64_t)op[5] * 1584614376596801529UL) + ((uint64_t)op[6] * 13002154218893442434UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 13002154218893442434UL) + ((uint64_t)op[1] * 5782342310597378278UL) + ((uint64_t)op[2] * 997445221193231493UL) + ((uint64_t)op[3] * 13375432650701855941UL) + ((uint64_t)op[4] * 1705753999475556693UL) + ((uint64_t)op[5] * 605019360372903003UL) + ((uint64_t)op[6] * 7923071882984007645UL);
	tmp_q[6] = ((uint64_t)op[0] * 1584614376596801529UL) + ((uint64_t)op[1] * 13002154218893442434UL) + ((uint64_t)op[2] * 5782342310597378278UL) + ((uint64_t)op[3] * 997445221193231493UL) + ((uint64_t)op[4] * 13375432650701855941UL) + ((uint64_t)op[5] * 1705753999475556693UL) + ((uint64_t)op[6] * 605019360372903003UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 63232498896L) + ((-((int128)tmp_q[1] * 3162344551L) - ((int128)tmp_q[2] * 47491041974L) + ((int128)tmp_q[3] * 33972260365L) + ((int128)tmp_q[4] * 14274052752L) - ((int128)tmp_q[5] * 27487306774L) + ((int128)tmp_q[6] * 21941799399L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 21941799399L) - ((int128)tmp_q[1] * 63232498896L) + ((-((int128)tmp_q[2] * 3162344551L) - ((int128)tmp_q[3] * 47491041974L) + ((int128)tmp_q[4] * 33972260365L) + ((int128)tmp_q[5] * 14274052752L) - ((int128)tmp_q[6] * 27487306774L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 27487306774L) + ((int128)tmp_q[1] * 21941799399L) - ((int128)tmp_q[2] * 63232498896L) + ((-((int128)tmp_q[3] * 3162344551L) - ((int128)tmp_q[4] * 47491041974L) + ((int128)tmp_q[5] * 33972260365L) + ((int128)tmp_q[6] * 14274052752L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 14274052752L) - ((int128)tmp_q[1] * 27487306774L) + ((int128)tmp_q[2] * 21941799399L) - ((int128)tmp_q[3] * 63232498896L) + ((-((int128)tmp_q[4] * 3162344551L) - ((int128)tmp_q[5] * 47491041974L) + ((int128)tmp_q[6] * 33972260365L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 33972260365L) + ((int128)tmp_q[1] * 14274052752L) - ((int128)tmp_q[2] * 27487306774L) + ((int128)tmp_q[3] * 21941799399L) - ((int128)tmp_q[4] * 63232498896L) + ((-((int128)tmp_q[5] * 3162344551L) - ((int128)tmp_q[6] * 47491041974L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 47491041974L) + ((int128)tmp_q[1] * 33972260365L) + ((int128)tmp_q[2] * 14274052752L) - ((int128)tmp_q[3] * 27487306774L) + ((int128)tmp_q[4] * 21941799399L) - ((int128)tmp_q[5] * 63232498896L) - ((int128)tmp_q[6] * 15811722755L);
	tmp_zero[6] = -((int128)tmp_q[0] * 3162344551L) - ((int128)tmp_q[1] * 47491041974L) + ((int128)tmp_q[2] * 33972260365L) + ((int128)tmp_q[3] * 14274052752L) - ((int128)tmp_q[4] * 27487306774L) + ((int128)tmp_q[5] * 21941799399L) - ((int128)tmp_q[6] * 63232498896L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

