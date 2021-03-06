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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4417629338939384745UL) + ((((uint64_t)op[1] * 5916825596024929524UL) + ((uint64_t)op[2] * 2057476723817462962UL) + ((uint64_t)op[3] * 13035376361473421173UL) + ((uint64_t)op[4] * 16576494055286377060UL) + ((uint64_t)op[5] * 2474003743431100546UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 2474003743431100546UL) + ((uint64_t)op[1] * 4417629338939384745UL) + ((((uint64_t)op[2] * 5916825596024929524UL) + ((uint64_t)op[3] * 2057476723817462962UL) + ((uint64_t)op[4] * 13035376361473421173UL) + ((uint64_t)op[5] * 16576494055286377060UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 16576494055286377060UL) + ((uint64_t)op[1] * 2474003743431100546UL) + ((uint64_t)op[2] * 4417629338939384745UL) + ((((uint64_t)op[3] * 5916825596024929524UL) + ((uint64_t)op[4] * 2057476723817462962UL) + ((uint64_t)op[5] * 13035376361473421173UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 13035376361473421173UL) + ((uint64_t)op[1] * 16576494055286377060UL) + ((uint64_t)op[2] * 2474003743431100546UL) + ((uint64_t)op[3] * 4417629338939384745UL) + ((((uint64_t)op[4] * 5916825596024929524UL) + ((uint64_t)op[5] * 2057476723817462962UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 2057476723817462962UL) + ((uint64_t)op[1] * 13035376361473421173UL) + ((uint64_t)op[2] * 16576494055286377060UL) + ((uint64_t)op[3] * 2474003743431100546UL) + ((uint64_t)op[4] * 4417629338939384745UL) + ((uint64_t)op[5] * 13226185763319385136UL);
	tmp_q[5] = ((uint64_t)op[0] * 5916825596024929524UL) + ((uint64_t)op[1] * 2057476723817462962UL) + ((uint64_t)op[2] * 13035376361473421173UL) + ((uint64_t)op[3] * 16576494055286377060UL) + ((uint64_t)op[4] * 2474003743431100546UL) + ((uint64_t)op[5] * 4417629338939384745UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 67640617339L) - ((((int128)tmp_q[1] * 4962311936L) + ((int128)tmp_q[2] * 89466452438L) - ((int128)tmp_q[3] * 120738820711L) + ((int128)tmp_q[4] * 24564750080L) + ((int128)tmp_q[5] * 12811667514L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 12811667514L) + ((int128)tmp_q[1] * 67640617339L) - ((((int128)tmp_q[2] * 4962311936L) + ((int128)tmp_q[3] * 89466452438L) - ((int128)tmp_q[4] * 120738820711L) + ((int128)tmp_q[5] * 24564750080L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 24564750080L) + ((int128)tmp_q[1] * 12811667514L) + ((int128)tmp_q[2] * 67640617339L) - ((((int128)tmp_q[3] * 4962311936L) + ((int128)tmp_q[4] * 89466452438L) - ((int128)tmp_q[5] * 120738820711L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 120738820711L) + ((int128)tmp_q[1] * 24564750080L) + ((int128)tmp_q[2] * 12811667514L) + ((int128)tmp_q[3] * 67640617339L) - ((((int128)tmp_q[4] * 4962311936L) + ((int128)tmp_q[5] * 89466452438L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 89466452438L) - ((int128)tmp_q[1] * 120738820711L) + ((int128)tmp_q[2] * 24564750080L) + ((int128)tmp_q[3] * 12811667514L) + ((int128)tmp_q[4] * 67640617339L) - ((int128)tmp_q[5] * 19849247744L);
	tmp_zero[5] = ((int128)tmp_q[0] * 4962311936L) + ((int128)tmp_q[1] * 89466452438L) - ((int128)tmp_q[2] * 120738820711L) + ((int128)tmp_q[3] * 24564750080L) + ((int128)tmp_q[4] * 12811667514L) + ((int128)tmp_q[5] * 67640617339L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

