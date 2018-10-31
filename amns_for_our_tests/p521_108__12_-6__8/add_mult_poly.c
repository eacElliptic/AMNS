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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 6);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 12);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 12);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 6);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16586577902293986413UL) + ((((uint64_t)op[1] * 8371076431049512751UL) + ((uint64_t)op[2] * 5644534091155847063UL) + ((uint64_t)op[3] * 13742810932566497955UL) + ((uint64_t)op[4] * 7509060529340663897UL) + ((uint64_t)op[5] * 6434915700725403400UL) + ((uint64_t)op[6] * 11375533929674314308UL) + ((uint64_t)op[7] * 6482590181570612186UL) + ((uint64_t)op[8] * 7186174428340229788UL) + ((uint64_t)op[9] * 13048132770825612550UL) + ((uint64_t)op[10] * 4113402058421614273UL) + ((uint64_t)op[11] * 13599984249787990699UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 13599984249787990699UL) + ((uint64_t)op[1] * 16586577902293986413UL) + ((((uint64_t)op[2] * 8371076431049512751UL) + ((uint64_t)op[3] * 5644534091155847063UL) + ((uint64_t)op[4] * 13742810932566497955UL) + ((uint64_t)op[5] * 7509060529340663897UL) + ((uint64_t)op[6] * 6434915700725403400UL) + ((uint64_t)op[7] * 11375533929674314308UL) + ((uint64_t)op[8] * 6482590181570612186UL) + ((uint64_t)op[9] * 7186174428340229788UL) + ((uint64_t)op[10] * 13048132770825612550UL) + ((uint64_t)op[11] * 4113402058421614273UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 4113402058421614273UL) + ((uint64_t)op[1] * 13599984249787990699UL) + ((uint64_t)op[2] * 16586577902293986413UL) + ((((uint64_t)op[3] * 8371076431049512751UL) + ((uint64_t)op[4] * 5644534091155847063UL) + ((uint64_t)op[5] * 13742810932566497955UL) + ((uint64_t)op[6] * 7509060529340663897UL) + ((uint64_t)op[7] * 6434915700725403400UL) + ((uint64_t)op[8] * 11375533929674314308UL) + ((uint64_t)op[9] * 6482590181570612186UL) + ((uint64_t)op[10] * 7186174428340229788UL) + ((uint64_t)op[11] * 13048132770825612550UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 13048132770825612550UL) + ((uint64_t)op[1] * 4113402058421614273UL) + ((uint64_t)op[2] * 13599984249787990699UL) + ((uint64_t)op[3] * 16586577902293986413UL) + ((((uint64_t)op[4] * 8371076431049512751UL) + ((uint64_t)op[5] * 5644534091155847063UL) + ((uint64_t)op[6] * 13742810932566497955UL) + ((uint64_t)op[7] * 7509060529340663897UL) + ((uint64_t)op[8] * 6434915700725403400UL) + ((uint64_t)op[9] * 11375533929674314308UL) + ((uint64_t)op[10] * 6482590181570612186UL) + ((uint64_t)op[11] * 7186174428340229788UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 7186174428340229788UL) + ((uint64_t)op[1] * 13048132770825612550UL) + ((uint64_t)op[2] * 4113402058421614273UL) + ((uint64_t)op[3] * 13599984249787990699UL) + ((uint64_t)op[4] * 16586577902293986413UL) + ((((uint64_t)op[5] * 8371076431049512751UL) + ((uint64_t)op[6] * 5644534091155847063UL) + ((uint64_t)op[7] * 13742810932566497955UL) + ((uint64_t)op[8] * 7509060529340663897UL) + ((uint64_t)op[9] * 6434915700725403400UL) + ((uint64_t)op[10] * 11375533929674314308UL) + ((uint64_t)op[11] * 6482590181570612186UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 6482590181570612186UL) + ((uint64_t)op[1] * 7186174428340229788UL) + ((uint64_t)op[2] * 13048132770825612550UL) + ((uint64_t)op[3] * 4113402058421614273UL) + ((uint64_t)op[4] * 13599984249787990699UL) + ((uint64_t)op[5] * 16586577902293986413UL) + ((((uint64_t)op[6] * 8371076431049512751UL) + ((uint64_t)op[7] * 5644534091155847063UL) + ((uint64_t)op[8] * 13742810932566497955UL) + ((uint64_t)op[9] * 7509060529340663897UL) + ((uint64_t)op[10] * 6434915700725403400UL) + ((uint64_t)op[11] * 11375533929674314308UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 11375533929674314308UL) + ((uint64_t)op[1] * 6482590181570612186UL) + ((uint64_t)op[2] * 7186174428340229788UL) + ((uint64_t)op[3] * 13048132770825612550UL) + ((uint64_t)op[4] * 4113402058421614273UL) + ((uint64_t)op[5] * 13599984249787990699UL) + ((uint64_t)op[6] * 16586577902293986413UL) + ((((uint64_t)op[7] * 8371076431049512751UL) + ((uint64_t)op[8] * 5644534091155847063UL) + ((uint64_t)op[9] * 13742810932566497955UL) + ((uint64_t)op[10] * 7509060529340663897UL) + ((uint64_t)op[11] * 6434915700725403400UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 6434915700725403400UL) + ((uint64_t)op[1] * 11375533929674314308UL) + ((uint64_t)op[2] * 6482590181570612186UL) + ((uint64_t)op[3] * 7186174428340229788UL) + ((uint64_t)op[4] * 13048132770825612550UL) + ((uint64_t)op[5] * 4113402058421614273UL) + ((uint64_t)op[6] * 13599984249787990699UL) + ((uint64_t)op[7] * 16586577902293986413UL) + ((((uint64_t)op[8] * 8371076431049512751UL) + ((uint64_t)op[9] * 5644534091155847063UL) + ((uint64_t)op[10] * 13742810932566497955UL) + ((uint64_t)op[11] * 7509060529340663897UL)) * 18446744073709551610);
	tmp_q[8] = ((uint64_t)op[0] * 7509060529340663897UL) + ((uint64_t)op[1] * 6434915700725403400UL) + ((uint64_t)op[2] * 11375533929674314308UL) + ((uint64_t)op[3] * 6482590181570612186UL) + ((uint64_t)op[4] * 7186174428340229788UL) + ((uint64_t)op[5] * 13048132770825612550UL) + ((uint64_t)op[6] * 4113402058421614273UL) + ((uint64_t)op[7] * 13599984249787990699UL) + ((uint64_t)op[8] * 16586577902293986413UL) + ((((uint64_t)op[9] * 8371076431049512751UL) + ((uint64_t)op[10] * 5644534091155847063UL) + ((uint64_t)op[11] * 13742810932566497955UL)) * 18446744073709551610);
	tmp_q[9] = ((uint64_t)op[0] * 13742810932566497955UL) + ((uint64_t)op[1] * 7509060529340663897UL) + ((uint64_t)op[2] * 6434915700725403400UL) + ((uint64_t)op[3] * 11375533929674314308UL) + ((uint64_t)op[4] * 6482590181570612186UL) + ((uint64_t)op[5] * 7186174428340229788UL) + ((uint64_t)op[6] * 13048132770825612550UL) + ((uint64_t)op[7] * 4113402058421614273UL) + ((uint64_t)op[8] * 13599984249787990699UL) + ((uint64_t)op[9] * 16586577902293986413UL) + ((((uint64_t)op[10] * 8371076431049512751UL) + ((uint64_t)op[11] * 5644534091155847063UL)) * 18446744073709551610);
	tmp_q[10] = ((uint64_t)op[0] * 5644534091155847063UL) + ((uint64_t)op[1] * 13742810932566497955UL) + ((uint64_t)op[2] * 7509060529340663897UL) + ((uint64_t)op[3] * 6434915700725403400UL) + ((uint64_t)op[4] * 11375533929674314308UL) + ((uint64_t)op[5] * 6482590181570612186UL) + ((uint64_t)op[6] * 7186174428340229788UL) + ((uint64_t)op[7] * 13048132770825612550UL) + ((uint64_t)op[8] * 4113402058421614273UL) + ((uint64_t)op[9] * 13599984249787990699UL) + ((uint64_t)op[10] * 16586577902293986413UL) + ((uint64_t)op[11] * 5113773634831578342UL);
	tmp_q[11] = ((uint64_t)op[0] * 8371076431049512751UL) + ((uint64_t)op[1] * 5644534091155847063UL) + ((uint64_t)op[2] * 13742810932566497955UL) + ((uint64_t)op[3] * 7509060529340663897UL) + ((uint64_t)op[4] * 6434915700725403400UL) + ((uint64_t)op[5] * 11375533929674314308UL) + ((uint64_t)op[6] * 6482590181570612186UL) + ((uint64_t)op[7] * 7186174428340229788UL) + ((uint64_t)op[8] * 13048132770825612550UL) + ((uint64_t)op[9] * 4113402058421614273UL) + ((uint64_t)op[10] * 13599984249787990699UL) + ((uint64_t)op[11] * 16586577902293986413UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4742939800031L) - ((((int128)tmp_q[1] * 4170175459214L) + ((int128)tmp_q[2] * 1979331652929L) + ((int128)tmp_q[3] * 2863729756160L) - ((int128)tmp_q[4] * 1747393599999L) - ((int128)tmp_q[5] * 325518559557L) - ((int128)tmp_q[6] * 1144383328749L) + ((int128)tmp_q[7] * 6024674647670L) + ((int128)tmp_q[8] * 1575548980185L) + ((int128)tmp_q[9] * 3055705258943L) - ((int128)tmp_q[10] * 2798222165386L) - ((int128)tmp_q[11] * 463785541333L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 463785541333L) + ((int128)tmp_q[1] * 4742939800031L) - ((((int128)tmp_q[2] * 4170175459214L) + ((int128)tmp_q[3] * 1979331652929L) + ((int128)tmp_q[4] * 2863729756160L) - ((int128)tmp_q[5] * 1747393599999L) - ((int128)tmp_q[6] * 325518559557L) - ((int128)tmp_q[7] * 1144383328749L) + ((int128)tmp_q[8] * 6024674647670L) + ((int128)tmp_q[9] * 1575548980185L) + ((int128)tmp_q[10] * 3055705258943L) - ((int128)tmp_q[11] * 2798222165386L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 2798222165386L) - ((int128)tmp_q[1] * 463785541333L) + ((int128)tmp_q[2] * 4742939800031L) - ((((int128)tmp_q[3] * 4170175459214L) + ((int128)tmp_q[4] * 1979331652929L) + ((int128)tmp_q[5] * 2863729756160L) - ((int128)tmp_q[6] * 1747393599999L) - ((int128)tmp_q[7] * 325518559557L) - ((int128)tmp_q[8] * 1144383328749L) + ((int128)tmp_q[9] * 6024674647670L) + ((int128)tmp_q[10] * 1575548980185L) + ((int128)tmp_q[11] * 3055705258943L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 3055705258943L) - ((int128)tmp_q[1] * 2798222165386L) - ((int128)tmp_q[2] * 463785541333L) + ((int128)tmp_q[3] * 4742939800031L) - ((((int128)tmp_q[4] * 4170175459214L) + ((int128)tmp_q[5] * 1979331652929L) + ((int128)tmp_q[6] * 2863729756160L) - ((int128)tmp_q[7] * 1747393599999L) - ((int128)tmp_q[8] * 325518559557L) - ((int128)tmp_q[9] * 1144383328749L) + ((int128)tmp_q[10] * 6024674647670L) + ((int128)tmp_q[11] * 1575548980185L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 1575548980185L) + ((int128)tmp_q[1] * 3055705258943L) - ((int128)tmp_q[2] * 2798222165386L) - ((int128)tmp_q[3] * 463785541333L) + ((int128)tmp_q[4] * 4742939800031L) - ((((int128)tmp_q[5] * 4170175459214L) + ((int128)tmp_q[6] * 1979331652929L) + ((int128)tmp_q[7] * 2863729756160L) - ((int128)tmp_q[8] * 1747393599999L) - ((int128)tmp_q[9] * 325518559557L) - ((int128)tmp_q[10] * 1144383328749L) + ((int128)tmp_q[11] * 6024674647670L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 6024674647670L) + ((int128)tmp_q[1] * 1575548980185L) + ((int128)tmp_q[2] * 3055705258943L) - ((int128)tmp_q[3] * 2798222165386L) - ((int128)tmp_q[4] * 463785541333L) + ((int128)tmp_q[5] * 4742939800031L) - ((((int128)tmp_q[6] * 4170175459214L) + ((int128)tmp_q[7] * 1979331652929L) + ((int128)tmp_q[8] * 2863729756160L) - ((int128)tmp_q[9] * 1747393599999L) - ((int128)tmp_q[10] * 325518559557L) - ((int128)tmp_q[11] * 1144383328749L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 1144383328749L) + ((int128)tmp_q[1] * 6024674647670L) + ((int128)tmp_q[2] * 1575548980185L) + ((int128)tmp_q[3] * 3055705258943L) - ((int128)tmp_q[4] * 2798222165386L) - ((int128)tmp_q[5] * 463785541333L) + ((int128)tmp_q[6] * 4742939800031L) - ((((int128)tmp_q[7] * 4170175459214L) + ((int128)tmp_q[8] * 1979331652929L) + ((int128)tmp_q[9] * 2863729756160L) - ((int128)tmp_q[10] * 1747393599999L) - ((int128)tmp_q[11] * 325518559557L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 325518559557L) - ((int128)tmp_q[1] * 1144383328749L) + ((int128)tmp_q[2] * 6024674647670L) + ((int128)tmp_q[3] * 1575548980185L) + ((int128)tmp_q[4] * 3055705258943L) - ((int128)tmp_q[5] * 2798222165386L) - ((int128)tmp_q[6] * 463785541333L) + ((int128)tmp_q[7] * 4742939800031L) - ((((int128)tmp_q[8] * 4170175459214L) + ((int128)tmp_q[9] * 1979331652929L) + ((int128)tmp_q[10] * 2863729756160L) - ((int128)tmp_q[11] * 1747393599999L)) * 6);
	tmp_zero[8] = -((int128)tmp_q[0] * 1747393599999L) - ((int128)tmp_q[1] * 325518559557L) - ((int128)tmp_q[2] * 1144383328749L) + ((int128)tmp_q[3] * 6024674647670L) + ((int128)tmp_q[4] * 1575548980185L) + ((int128)tmp_q[5] * 3055705258943L) - ((int128)tmp_q[6] * 2798222165386L) - ((int128)tmp_q[7] * 463785541333L) + ((int128)tmp_q[8] * 4742939800031L) - ((((int128)tmp_q[9] * 4170175459214L) + ((int128)tmp_q[10] * 1979331652929L) + ((int128)tmp_q[11] * 2863729756160L)) * 6);
	tmp_zero[9] = ((int128)tmp_q[0] * 2863729756160L) - ((int128)tmp_q[1] * 1747393599999L) - ((int128)tmp_q[2] * 325518559557L) - ((int128)tmp_q[3] * 1144383328749L) + ((int128)tmp_q[4] * 6024674647670L) + ((int128)tmp_q[5] * 1575548980185L) + ((int128)tmp_q[6] * 3055705258943L) - ((int128)tmp_q[7] * 2798222165386L) - ((int128)tmp_q[8] * 463785541333L) + ((int128)tmp_q[9] * 4742939800031L) - ((((int128)tmp_q[10] * 4170175459214L) + ((int128)tmp_q[11] * 1979331652929L)) * 6);
	tmp_zero[10] = ((int128)tmp_q[0] * 1979331652929L) + ((int128)tmp_q[1] * 2863729756160L) - ((int128)tmp_q[2] * 1747393599999L) - ((int128)tmp_q[3] * 325518559557L) - ((int128)tmp_q[4] * 1144383328749L) + ((int128)tmp_q[5] * 6024674647670L) + ((int128)tmp_q[6] * 1575548980185L) + ((int128)tmp_q[7] * 3055705258943L) - ((int128)tmp_q[8] * 2798222165386L) - ((int128)tmp_q[9] * 463785541333L) + ((int128)tmp_q[10] * 4742939800031L) - ((int128)tmp_q[11] * 25021052755284L);
	tmp_zero[11] = ((int128)tmp_q[0] * 4170175459214L) + ((int128)tmp_q[1] * 1979331652929L) + ((int128)tmp_q[2] * 2863729756160L) - ((int128)tmp_q[3] * 1747393599999L) - ((int128)tmp_q[4] * 325518559557L) - ((int128)tmp_q[5] * 1144383328749L) + ((int128)tmp_q[6] * 6024674647670L) + ((int128)tmp_q[7] * 1575548980185L) + ((int128)tmp_q[8] * 3055705258943L) - ((int128)tmp_q[9] * 2798222165386L) - ((int128)tmp_q[10] * 463785541333L) + ((int128)tmp_q[11] * 4742939800031L);

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
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

