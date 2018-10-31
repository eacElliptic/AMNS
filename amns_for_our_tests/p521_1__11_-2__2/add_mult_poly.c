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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11546390512710933485UL) + ((((uint64_t)op[1] * 16043893474339552361UL) + ((uint64_t)op[2] * 13380327731213645298UL) + ((uint64_t)op[3] * 3433701996778468474UL) + ((uint64_t)op[4] * 7218873199623368031UL) + ((uint64_t)op[5] * 2117895820833350699UL) + ((uint64_t)op[6] * 16813832461427639530UL) + ((uint64_t)op[7] * 15066400468617644908UL) + ((uint64_t)op[8] * 2118195139002276687UL) + ((uint64_t)op[9] * 9670333886905090383UL) + ((uint64_t)op[10] * 12589051207877381941UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 12589051207877381941UL) + ((uint64_t)op[1] * 11546390512710933485UL) + ((((uint64_t)op[2] * 16043893474339552361UL) + ((uint64_t)op[3] * 13380327731213645298UL) + ((uint64_t)op[4] * 3433701996778468474UL) + ((uint64_t)op[5] * 7218873199623368031UL) + ((uint64_t)op[6] * 2117895820833350699UL) + ((uint64_t)op[7] * 16813832461427639530UL) + ((uint64_t)op[8] * 15066400468617644908UL) + ((uint64_t)op[9] * 2118195139002276687UL) + ((uint64_t)op[10] * 9670333886905090383UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 9670333886905090383UL) + ((uint64_t)op[1] * 12589051207877381941UL) + ((uint64_t)op[2] * 11546390512710933485UL) + ((((uint64_t)op[3] * 16043893474339552361UL) + ((uint64_t)op[4] * 13380327731213645298UL) + ((uint64_t)op[5] * 3433701996778468474UL) + ((uint64_t)op[6] * 7218873199623368031UL) + ((uint64_t)op[7] * 2117895820833350699UL) + ((uint64_t)op[8] * 16813832461427639530UL) + ((uint64_t)op[9] * 15066400468617644908UL) + ((uint64_t)op[10] * 2118195139002276687UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 2118195139002276687UL) + ((uint64_t)op[1] * 9670333886905090383UL) + ((uint64_t)op[2] * 12589051207877381941UL) + ((uint64_t)op[3] * 11546390512710933485UL) + ((((uint64_t)op[4] * 16043893474339552361UL) + ((uint64_t)op[5] * 13380327731213645298UL) + ((uint64_t)op[6] * 3433701996778468474UL) + ((uint64_t)op[7] * 7218873199623368031UL) + ((uint64_t)op[8] * 2117895820833350699UL) + ((uint64_t)op[9] * 16813832461427639530UL) + ((uint64_t)op[10] * 15066400468617644908UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 15066400468617644908UL) + ((uint64_t)op[1] * 2118195139002276687UL) + ((uint64_t)op[2] * 9670333886905090383UL) + ((uint64_t)op[3] * 12589051207877381941UL) + ((uint64_t)op[4] * 11546390512710933485UL) + ((((uint64_t)op[5] * 16043893474339552361UL) + ((uint64_t)op[6] * 13380327731213645298UL) + ((uint64_t)op[7] * 3433701996778468474UL) + ((uint64_t)op[8] * 7218873199623368031UL) + ((uint64_t)op[9] * 2117895820833350699UL) + ((uint64_t)op[10] * 16813832461427639530UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 16813832461427639530UL) + ((uint64_t)op[1] * 15066400468617644908UL) + ((uint64_t)op[2] * 2118195139002276687UL) + ((uint64_t)op[3] * 9670333886905090383UL) + ((uint64_t)op[4] * 12589051207877381941UL) + ((uint64_t)op[5] * 11546390512710933485UL) + ((((uint64_t)op[6] * 16043893474339552361UL) + ((uint64_t)op[7] * 13380327731213645298UL) + ((uint64_t)op[8] * 3433701996778468474UL) + ((uint64_t)op[9] * 7218873199623368031UL) + ((uint64_t)op[10] * 2117895820833350699UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 2117895820833350699UL) + ((uint64_t)op[1] * 16813832461427639530UL) + ((uint64_t)op[2] * 15066400468617644908UL) + ((uint64_t)op[3] * 2118195139002276687UL) + ((uint64_t)op[4] * 9670333886905090383UL) + ((uint64_t)op[5] * 12589051207877381941UL) + ((uint64_t)op[6] * 11546390512710933485UL) + ((((uint64_t)op[7] * 16043893474339552361UL) + ((uint64_t)op[8] * 13380327731213645298UL) + ((uint64_t)op[9] * 3433701996778468474UL) + ((uint64_t)op[10] * 7218873199623368031UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 7218873199623368031UL) + ((uint64_t)op[1] * 2117895820833350699UL) + ((uint64_t)op[2] * 16813832461427639530UL) + ((uint64_t)op[3] * 15066400468617644908UL) + ((uint64_t)op[4] * 2118195139002276687UL) + ((uint64_t)op[5] * 9670333886905090383UL) + ((uint64_t)op[6] * 12589051207877381941UL) + ((uint64_t)op[7] * 11546390512710933485UL) + ((((uint64_t)op[8] * 16043893474339552361UL) + ((uint64_t)op[9] * 13380327731213645298UL) + ((uint64_t)op[10] * 3433701996778468474UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 3433701996778468474UL) + ((uint64_t)op[1] * 7218873199623368031UL) + ((uint64_t)op[2] * 2117895820833350699UL) + ((uint64_t)op[3] * 16813832461427639530UL) + ((uint64_t)op[4] * 15066400468617644908UL) + ((uint64_t)op[5] * 2118195139002276687UL) + ((uint64_t)op[6] * 9670333886905090383UL) + ((uint64_t)op[7] * 12589051207877381941UL) + ((uint64_t)op[8] * 11546390512710933485UL) + ((((uint64_t)op[9] * 16043893474339552361UL) + ((uint64_t)op[10] * 13380327731213645298UL)) * 18446744073709551614);
	tmp_q[9] = ((uint64_t)op[0] * 13380327731213645298UL) + ((uint64_t)op[1] * 3433701996778468474UL) + ((uint64_t)op[2] * 7218873199623368031UL) + ((uint64_t)op[3] * 2117895820833350699UL) + ((uint64_t)op[4] * 16813832461427639530UL) + ((uint64_t)op[5] * 15066400468617644908UL) + ((uint64_t)op[6] * 2118195139002276687UL) + ((uint64_t)op[7] * 9670333886905090383UL) + ((uint64_t)op[8] * 12589051207877381941UL) + ((uint64_t)op[9] * 11546390512710933485UL) + ((uint64_t)op[10] * 4805701198739998510UL);
	tmp_q[10] = ((uint64_t)op[0] * 16043893474339552361UL) + ((uint64_t)op[1] * 13380327731213645298UL) + ((uint64_t)op[2] * 3433701996778468474UL) + ((uint64_t)op[3] * 7218873199623368031UL) + ((uint64_t)op[4] * 2117895820833350699UL) + ((uint64_t)op[5] * 16813832461427639530UL) + ((uint64_t)op[6] * 15066400468617644908UL) + ((uint64_t)op[7] * 2118195139002276687UL) + ((uint64_t)op[8] * 9670333886905090383UL) + ((uint64_t)op[9] * 12589051207877381941UL) + ((uint64_t)op[10] * 11546390512710933485UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 87636852303333L) - ((((int128)tmp_q[1] * 11574255093607L) - ((int128)tmp_q[2] * 1048400507752L) - ((int128)tmp_q[3] * 105569470142256L) - ((int128)tmp_q[4] * 10303553938673L) - ((int128)tmp_q[5] * 51037885194973L) - ((int128)tmp_q[6] * 7188560526295L) - ((int128)tmp_q[7] * 62748404245801L) + ((int128)tmp_q[8] * 7746748493436L) + ((int128)tmp_q[9] * 67080104432018L) + ((int128)tmp_q[10] * 73497300892739L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 73497300892739L) - ((int128)tmp_q[1] * 87636852303333L) - ((((int128)tmp_q[2] * 11574255093607L) - ((int128)tmp_q[3] * 1048400507752L) - ((int128)tmp_q[4] * 105569470142256L) - ((int128)tmp_q[5] * 10303553938673L) - ((int128)tmp_q[6] * 51037885194973L) - ((int128)tmp_q[7] * 7188560526295L) - ((int128)tmp_q[8] * 62748404245801L) + ((int128)tmp_q[9] * 7746748493436L) + ((int128)tmp_q[10] * 67080104432018L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 67080104432018L) + ((int128)tmp_q[1] * 73497300892739L) - ((int128)tmp_q[2] * 87636852303333L) - ((((int128)tmp_q[3] * 11574255093607L) - ((int128)tmp_q[4] * 1048400507752L) - ((int128)tmp_q[5] * 105569470142256L) - ((int128)tmp_q[6] * 10303553938673L) - ((int128)tmp_q[7] * 51037885194973L) - ((int128)tmp_q[8] * 7188560526295L) - ((int128)tmp_q[9] * 62748404245801L) + ((int128)tmp_q[10] * 7746748493436L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 7746748493436L) + ((int128)tmp_q[1] * 67080104432018L) + ((int128)tmp_q[2] * 73497300892739L) - ((int128)tmp_q[3] * 87636852303333L) - ((((int128)tmp_q[4] * 11574255093607L) - ((int128)tmp_q[5] * 1048400507752L) - ((int128)tmp_q[6] * 105569470142256L) - ((int128)tmp_q[7] * 10303553938673L) - ((int128)tmp_q[8] * 51037885194973L) - ((int128)tmp_q[9] * 7188560526295L) - ((int128)tmp_q[10] * 62748404245801L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 62748404245801L) + ((int128)tmp_q[1] * 7746748493436L) + ((int128)tmp_q[2] * 67080104432018L) + ((int128)tmp_q[3] * 73497300892739L) - ((int128)tmp_q[4] * 87636852303333L) - ((((int128)tmp_q[5] * 11574255093607L) - ((int128)tmp_q[6] * 1048400507752L) - ((int128)tmp_q[7] * 105569470142256L) - ((int128)tmp_q[8] * 10303553938673L) - ((int128)tmp_q[9] * 51037885194973L) - ((int128)tmp_q[10] * 7188560526295L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 7188560526295L) - ((int128)tmp_q[1] * 62748404245801L) + ((int128)tmp_q[2] * 7746748493436L) + ((int128)tmp_q[3] * 67080104432018L) + ((int128)tmp_q[4] * 73497300892739L) - ((int128)tmp_q[5] * 87636852303333L) - ((((int128)tmp_q[6] * 11574255093607L) - ((int128)tmp_q[7] * 1048400507752L) - ((int128)tmp_q[8] * 105569470142256L) - ((int128)tmp_q[9] * 10303553938673L) - ((int128)tmp_q[10] * 51037885194973L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 51037885194973L) - ((int128)tmp_q[1] * 7188560526295L) - ((int128)tmp_q[2] * 62748404245801L) + ((int128)tmp_q[3] * 7746748493436L) + ((int128)tmp_q[4] * 67080104432018L) + ((int128)tmp_q[5] * 73497300892739L) - ((int128)tmp_q[6] * 87636852303333L) - ((((int128)tmp_q[7] * 11574255093607L) - ((int128)tmp_q[8] * 1048400507752L) - ((int128)tmp_q[9] * 105569470142256L) - ((int128)tmp_q[10] * 10303553938673L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 10303553938673L) - ((int128)tmp_q[1] * 51037885194973L) - ((int128)tmp_q[2] * 7188560526295L) - ((int128)tmp_q[3] * 62748404245801L) + ((int128)tmp_q[4] * 7746748493436L) + ((int128)tmp_q[5] * 67080104432018L) + ((int128)tmp_q[6] * 73497300892739L) - ((int128)tmp_q[7] * 87636852303333L) - ((((int128)tmp_q[8] * 11574255093607L) - ((int128)tmp_q[9] * 1048400507752L) - ((int128)tmp_q[10] * 105569470142256L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 105569470142256L) - ((int128)tmp_q[1] * 10303553938673L) - ((int128)tmp_q[2] * 51037885194973L) - ((int128)tmp_q[3] * 7188560526295L) - ((int128)tmp_q[4] * 62748404245801L) + ((int128)tmp_q[5] * 7746748493436L) + ((int128)tmp_q[6] * 67080104432018L) + ((int128)tmp_q[7] * 73497300892739L) - ((int128)tmp_q[8] * 87636852303333L) - ((((int128)tmp_q[9] * 11574255093607L) - ((int128)tmp_q[10] * 1048400507752L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 1048400507752L) - ((int128)tmp_q[1] * 105569470142256L) - ((int128)tmp_q[2] * 10303553938673L) - ((int128)tmp_q[3] * 51037885194973L) - ((int128)tmp_q[4] * 7188560526295L) - ((int128)tmp_q[5] * 62748404245801L) + ((int128)tmp_q[6] * 7746748493436L) + ((int128)tmp_q[7] * 67080104432018L) + ((int128)tmp_q[8] * 73497300892739L) - ((int128)tmp_q[9] * 87636852303333L) - ((int128)tmp_q[10] * 23148510187214L);
	tmp_zero[10] = ((int128)tmp_q[0] * 11574255093607L) - ((int128)tmp_q[1] * 1048400507752L) - ((int128)tmp_q[2] * 105569470142256L) - ((int128)tmp_q[3] * 10303553938673L) - ((int128)tmp_q[4] * 51037885194973L) - ((int128)tmp_q[5] * 7188560526295L) - ((int128)tmp_q[6] * 62748404245801L) + ((int128)tmp_q[7] * 7746748493436L) + ((int128)tmp_q[8] * 67080104432018L) + ((int128)tmp_q[9] * 73497300892739L) - ((int128)tmp_q[10] * 87636852303333L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
}

