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
	tmp_q[0] = ((uint64_t)op[0] * 17066428720457971153UL) + ((((uint64_t)op[1] * 9771631849464649710UL) + ((uint64_t)op[2] * 4871386358000886865UL) + ((uint64_t)op[3] * 1119713460644651599UL) + ((uint64_t)op[4] * 11507238027146437049UL) + ((uint64_t)op[5] * 18013589229468155765UL) + ((uint64_t)op[6] * 3480434409091962808UL) + ((uint64_t)op[7] * 6442054633437066790UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 6442054633437066790UL) + ((uint64_t)op[1] * 17066428720457971153UL) + ((((uint64_t)op[2] * 9771631849464649710UL) + ((uint64_t)op[3] * 4871386358000886865UL) + ((uint64_t)op[4] * 1119713460644651599UL) + ((uint64_t)op[5] * 11507238027146437049UL) + ((uint64_t)op[6] * 18013589229468155765UL) + ((uint64_t)op[7] * 3480434409091962808UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 3480434409091962808UL) + ((uint64_t)op[1] * 6442054633437066790UL) + ((uint64_t)op[2] * 17066428720457971153UL) + ((((uint64_t)op[3] * 9771631849464649710UL) + ((uint64_t)op[4] * 4871386358000886865UL) + ((uint64_t)op[5] * 1119713460644651599UL) + ((uint64_t)op[6] * 11507238027146437049UL) + ((uint64_t)op[7] * 18013589229468155765UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 18013589229468155765UL) + ((uint64_t)op[1] * 3480434409091962808UL) + ((uint64_t)op[2] * 6442054633437066790UL) + ((uint64_t)op[3] * 17066428720457971153UL) + ((((uint64_t)op[4] * 9771631849464649710UL) + ((uint64_t)op[5] * 4871386358000886865UL) + ((uint64_t)op[6] * 1119713460644651599UL) + ((uint64_t)op[7] * 11507238027146437049UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 11507238027146437049UL) + ((uint64_t)op[1] * 18013589229468155765UL) + ((uint64_t)op[2] * 3480434409091962808UL) + ((uint64_t)op[3] * 6442054633437066790UL) + ((uint64_t)op[4] * 17066428720457971153UL) + ((((uint64_t)op[5] * 9771631849464649710UL) + ((uint64_t)op[6] * 4871386358000886865UL) + ((uint64_t)op[7] * 1119713460644651599UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 1119713460644651599UL) + ((uint64_t)op[1] * 11507238027146437049UL) + ((uint64_t)op[2] * 18013589229468155765UL) + ((uint64_t)op[3] * 3480434409091962808UL) + ((uint64_t)op[4] * 6442054633437066790UL) + ((uint64_t)op[5] * 17066428720457971153UL) + ((((uint64_t)op[6] * 9771631849464649710UL) + ((uint64_t)op[7] * 4871386358000886865UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 4871386358000886865UL) + ((uint64_t)op[1] * 1119713460644651599UL) + ((uint64_t)op[2] * 11507238027146437049UL) + ((uint64_t)op[3] * 18013589229468155765UL) + ((uint64_t)op[4] * 3480434409091962808UL) + ((uint64_t)op[5] * 6442054633437066790UL) + ((uint64_t)op[6] * 17066428720457971153UL) + ((uint64_t)op[7] * 4386078500878991216UL);
	tmp_q[7] = ((uint64_t)op[0] * 9771631849464649710UL) + ((uint64_t)op[1] * 4871386358000886865UL) + ((uint64_t)op[2] * 1119713460644651599UL) + ((uint64_t)op[3] * 11507238027146437049UL) + ((uint64_t)op[4] * 18013589229468155765UL) + ((uint64_t)op[5] * 3480434409091962808UL) + ((uint64_t)op[6] * 6442054633437066790UL) + ((uint64_t)op[7] * 17066428720457971153UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 43517420736599L) + ((-((int128)tmp_q[1] * 23176430282378L) - ((int128)tmp_q[2] * 46049970125920L) + ((int128)tmp_q[3] * 89416433352495L) + ((int128)tmp_q[4] * 41442138400653L) - ((int128)tmp_q[5] * 167635377767923L) + ((int128)tmp_q[6] * 169364453943796L) + ((int128)tmp_q[7] * 109695672309182L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 109695672309182L) + ((int128)tmp_q[1] * 43517420736599L) + ((-((int128)tmp_q[2] * 23176430282378L) - ((int128)tmp_q[3] * 46049970125920L) + ((int128)tmp_q[4] * 89416433352495L) + ((int128)tmp_q[5] * 41442138400653L) - ((int128)tmp_q[6] * 167635377767923L) + ((int128)tmp_q[7] * 169364453943796L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 169364453943796L) + ((int128)tmp_q[1] * 109695672309182L) + ((int128)tmp_q[2] * 43517420736599L) + ((-((int128)tmp_q[3] * 23176430282378L) - ((int128)tmp_q[4] * 46049970125920L) + ((int128)tmp_q[5] * 89416433352495L) + ((int128)tmp_q[6] * 41442138400653L) - ((int128)tmp_q[7] * 167635377767923L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 167635377767923L) + ((int128)tmp_q[1] * 169364453943796L) + ((int128)tmp_q[2] * 109695672309182L) + ((int128)tmp_q[3] * 43517420736599L) + ((-((int128)tmp_q[4] * 23176430282378L) - ((int128)tmp_q[5] * 46049970125920L) + ((int128)tmp_q[6] * 89416433352495L) + ((int128)tmp_q[7] * 41442138400653L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 41442138400653L) - ((int128)tmp_q[1] * 167635377767923L) + ((int128)tmp_q[2] * 169364453943796L) + ((int128)tmp_q[3] * 109695672309182L) + ((int128)tmp_q[4] * 43517420736599L) + ((-((int128)tmp_q[5] * 23176430282378L) - ((int128)tmp_q[6] * 46049970125920L) + ((int128)tmp_q[7] * 89416433352495L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 89416433352495L) + ((int128)tmp_q[1] * 41442138400653L) - ((int128)tmp_q[2] * 167635377767923L) + ((int128)tmp_q[3] * 169364453943796L) + ((int128)tmp_q[4] * 109695672309182L) + ((int128)tmp_q[5] * 43517420736599L) + ((-((int128)tmp_q[6] * 23176430282378L) - ((int128)tmp_q[7] * 46049970125920L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 46049970125920L) + ((int128)tmp_q[1] * 89416433352495L) + ((int128)tmp_q[2] * 41442138400653L) - ((int128)tmp_q[3] * 167635377767923L) + ((int128)tmp_q[4] * 169364453943796L) + ((int128)tmp_q[5] * 109695672309182L) + ((int128)tmp_q[6] * 43517420736599L) - ((int128)tmp_q[7] * 185411442259024L);
	tmp_zero[7] = -((int128)tmp_q[0] * 23176430282378L) - ((int128)tmp_q[1] * 46049970125920L) + ((int128)tmp_q[2] * 89416433352495L) + ((int128)tmp_q[3] * 41442138400653L) - ((int128)tmp_q[4] * 167635377767923L) + ((int128)tmp_q[5] * 169364453943796L) + ((int128)tmp_q[6] * 109695672309182L) + ((int128)tmp_q[7] * 43517420736599L);

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

