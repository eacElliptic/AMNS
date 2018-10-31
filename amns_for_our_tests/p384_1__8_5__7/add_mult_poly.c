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
	tmp_q[0] = ((uint64_t)op[0] * 13431428614945706983UL) + ((((uint64_t)op[1] * 6543086642243262747UL) + ((uint64_t)op[2] * 14102767249534476159UL) + ((uint64_t)op[3] * 12485043677022581432UL) + ((uint64_t)op[4] * 14382495900998386669UL) + ((uint64_t)op[5] * 10964813254738783590UL) + ((uint64_t)op[6] * 11048369016878479399UL) + ((uint64_t)op[7] * 17700457578340067660UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 17700457578340067660UL) + ((uint64_t)op[1] * 13431428614945706983UL) + ((((uint64_t)op[2] * 6543086642243262747UL) + ((uint64_t)op[3] * 14102767249534476159UL) + ((uint64_t)op[4] * 12485043677022581432UL) + ((uint64_t)op[5] * 14382495900998386669UL) + ((uint64_t)op[6] * 10964813254738783590UL) + ((uint64_t)op[7] * 11048369016878479399UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 11048369016878479399UL) + ((uint64_t)op[1] * 17700457578340067660UL) + ((uint64_t)op[2] * 13431428614945706983UL) + ((((uint64_t)op[3] * 6543086642243262747UL) + ((uint64_t)op[4] * 14102767249534476159UL) + ((uint64_t)op[5] * 12485043677022581432UL) + ((uint64_t)op[6] * 14382495900998386669UL) + ((uint64_t)op[7] * 10964813254738783590UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 10964813254738783590UL) + ((uint64_t)op[1] * 11048369016878479399UL) + ((uint64_t)op[2] * 17700457578340067660UL) + ((uint64_t)op[3] * 13431428614945706983UL) + ((((uint64_t)op[4] * 6543086642243262747UL) + ((uint64_t)op[5] * 14102767249534476159UL) + ((uint64_t)op[6] * 12485043677022581432UL) + ((uint64_t)op[7] * 14382495900998386669UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 14382495900998386669UL) + ((uint64_t)op[1] * 10964813254738783590UL) + ((uint64_t)op[2] * 11048369016878479399UL) + ((uint64_t)op[3] * 17700457578340067660UL) + ((uint64_t)op[4] * 13431428614945706983UL) + ((((uint64_t)op[5] * 6543086642243262747UL) + ((uint64_t)op[6] * 14102767249534476159UL) + ((uint64_t)op[7] * 12485043677022581432UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 12485043677022581432UL) + ((uint64_t)op[1] * 14382495900998386669UL) + ((uint64_t)op[2] * 10964813254738783590UL) + ((uint64_t)op[3] * 11048369016878479399UL) + ((uint64_t)op[4] * 17700457578340067660UL) + ((uint64_t)op[5] * 13431428614945706983UL) + ((((uint64_t)op[6] * 6543086642243262747UL) + ((uint64_t)op[7] * 14102767249534476159UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 14102767249534476159UL) + ((uint64_t)op[1] * 12485043677022581432UL) + ((uint64_t)op[2] * 14382495900998386669UL) + ((uint64_t)op[3] * 10964813254738783590UL) + ((uint64_t)op[4] * 11048369016878479399UL) + ((uint64_t)op[5] * 17700457578340067660UL) + ((uint64_t)op[6] * 13431428614945706983UL) + ((uint64_t)op[7] * 14268689137506762119UL);
	tmp_q[7] = ((uint64_t)op[0] * 6543086642243262747UL) + ((uint64_t)op[1] * 14102767249534476159UL) + ((uint64_t)op[2] * 12485043677022581432UL) + ((uint64_t)op[3] * 14382495900998386669UL) + ((uint64_t)op[4] * 10964813254738783590UL) + ((uint64_t)op[5] * 11048369016878479399UL) + ((uint64_t)op[6] * 17700457578340067660UL) + ((uint64_t)op[7] * 13431428614945706983UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 77638088259009L) + ((((int128)tmp_q[1] * 125654209935620L) + ((int128)tmp_q[2] * 79836819430313L) - ((int128)tmp_q[3] * 57148271461138L) - ((int128)tmp_q[4] * 162562855408829L) + ((int128)tmp_q[5] * 61169911749964L) + ((int128)tmp_q[6] * 41093109675043L) + ((int128)tmp_q[7] * 173270000508261L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 173270000508261L) - ((int128)tmp_q[1] * 77638088259009L) + ((((int128)tmp_q[2] * 125654209935620L) + ((int128)tmp_q[3] * 79836819430313L) - ((int128)tmp_q[4] * 57148271461138L) - ((int128)tmp_q[5] * 162562855408829L) + ((int128)tmp_q[6] * 61169911749964L) + ((int128)tmp_q[7] * 41093109675043L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 41093109675043L) + ((int128)tmp_q[1] * 173270000508261L) - ((int128)tmp_q[2] * 77638088259009L) + ((((int128)tmp_q[3] * 125654209935620L) + ((int128)tmp_q[4] * 79836819430313L) - ((int128)tmp_q[5] * 57148271461138L) - ((int128)tmp_q[6] * 162562855408829L) + ((int128)tmp_q[7] * 61169911749964L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 61169911749964L) + ((int128)tmp_q[1] * 41093109675043L) + ((int128)tmp_q[2] * 173270000508261L) - ((int128)tmp_q[3] * 77638088259009L) + ((((int128)tmp_q[4] * 125654209935620L) + ((int128)tmp_q[5] * 79836819430313L) - ((int128)tmp_q[6] * 57148271461138L) - ((int128)tmp_q[7] * 162562855408829L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 162562855408829L) + ((int128)tmp_q[1] * 61169911749964L) + ((int128)tmp_q[2] * 41093109675043L) + ((int128)tmp_q[3] * 173270000508261L) - ((int128)tmp_q[4] * 77638088259009L) + ((((int128)tmp_q[5] * 125654209935620L) + ((int128)tmp_q[6] * 79836819430313L) - ((int128)tmp_q[7] * 57148271461138L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 57148271461138L) - ((int128)tmp_q[1] * 162562855408829L) + ((int128)tmp_q[2] * 61169911749964L) + ((int128)tmp_q[3] * 41093109675043L) + ((int128)tmp_q[4] * 173270000508261L) - ((int128)tmp_q[5] * 77638088259009L) + ((((int128)tmp_q[6] * 125654209935620L) + ((int128)tmp_q[7] * 79836819430313L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 79836819430313L) - ((int128)tmp_q[1] * 57148271461138L) - ((int128)tmp_q[2] * 162562855408829L) + ((int128)tmp_q[3] * 61169911749964L) + ((int128)tmp_q[4] * 41093109675043L) + ((int128)tmp_q[5] * 173270000508261L) - ((int128)tmp_q[6] * 77638088259009L) + ((int128)tmp_q[7] * 628271049678100L);
	tmp_zero[7] = ((int128)tmp_q[0] * 125654209935620L) + ((int128)tmp_q[1] * 79836819430313L) - ((int128)tmp_q[2] * 57148271461138L) - ((int128)tmp_q[3] * 162562855408829L) + ((int128)tmp_q[4] * 61169911749964L) + ((int128)tmp_q[5] * 41093109675043L) + ((int128)tmp_q[6] * 173270000508261L) - ((int128)tmp_q[7] * 77638088259009L);

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

