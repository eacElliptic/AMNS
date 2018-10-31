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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10883422459652406185UL) + ((((uint64_t)op[1] * 5342609006528468788UL) + ((uint64_t)op[2] * 6739752439300693295UL) + ((uint64_t)op[3] * 17755411103475981349UL) + ((uint64_t)op[4] * 1196851098066101627UL) + ((uint64_t)op[5] * 5522947523484628917UL) + ((uint64_t)op[6] * 7399264221176564053UL) + ((uint64_t)op[7] * 2785037492349527960UL) + ((uint64_t)op[8] * 16839937630611475622UL) + ((uint64_t)op[9] * 15499340231624087516UL) + ((uint64_t)op[10] * 12522298333479763018UL) + ((uint64_t)op[11] * 12113416962807399125UL)) * 14);
	tmp_q[1] = ((uint64_t)op[0] * 12113416962807399125UL) + ((uint64_t)op[1] * 10883422459652406185UL) + ((((uint64_t)op[2] * 5342609006528468788UL) + ((uint64_t)op[3] * 6739752439300693295UL) + ((uint64_t)op[4] * 17755411103475981349UL) + ((uint64_t)op[5] * 1196851098066101627UL) + ((uint64_t)op[6] * 5522947523484628917UL) + ((uint64_t)op[7] * 7399264221176564053UL) + ((uint64_t)op[8] * 2785037492349527960UL) + ((uint64_t)op[9] * 16839937630611475622UL) + ((uint64_t)op[10] * 15499340231624087516UL) + ((uint64_t)op[11] * 12522298333479763018UL)) * 14);
	tmp_q[2] = ((uint64_t)op[0] * 12522298333479763018UL) + ((uint64_t)op[1] * 12113416962807399125UL) + ((uint64_t)op[2] * 10883422459652406185UL) + ((((uint64_t)op[3] * 5342609006528468788UL) + ((uint64_t)op[4] * 6739752439300693295UL) + ((uint64_t)op[5] * 17755411103475981349UL) + ((uint64_t)op[6] * 1196851098066101627UL) + ((uint64_t)op[7] * 5522947523484628917UL) + ((uint64_t)op[8] * 7399264221176564053UL) + ((uint64_t)op[9] * 2785037492349527960UL) + ((uint64_t)op[10] * 16839937630611475622UL) + ((uint64_t)op[11] * 15499340231624087516UL)) * 14);
	tmp_q[3] = ((uint64_t)op[0] * 15499340231624087516UL) + ((uint64_t)op[1] * 12522298333479763018UL) + ((uint64_t)op[2] * 12113416962807399125UL) + ((uint64_t)op[3] * 10883422459652406185UL) + ((((uint64_t)op[4] * 5342609006528468788UL) + ((uint64_t)op[5] * 6739752439300693295UL) + ((uint64_t)op[6] * 17755411103475981349UL) + ((uint64_t)op[7] * 1196851098066101627UL) + ((uint64_t)op[8] * 5522947523484628917UL) + ((uint64_t)op[9] * 7399264221176564053UL) + ((uint64_t)op[10] * 2785037492349527960UL) + ((uint64_t)op[11] * 16839937630611475622UL)) * 14);
	tmp_q[4] = ((uint64_t)op[0] * 16839937630611475622UL) + ((uint64_t)op[1] * 15499340231624087516UL) + ((uint64_t)op[2] * 12522298333479763018UL) + ((uint64_t)op[3] * 12113416962807399125UL) + ((uint64_t)op[4] * 10883422459652406185UL) + ((((uint64_t)op[5] * 5342609006528468788UL) + ((uint64_t)op[6] * 6739752439300693295UL) + ((uint64_t)op[7] * 17755411103475981349UL) + ((uint64_t)op[8] * 1196851098066101627UL) + ((uint64_t)op[9] * 5522947523484628917UL) + ((uint64_t)op[10] * 7399264221176564053UL) + ((uint64_t)op[11] * 2785037492349527960UL)) * 14);
	tmp_q[5] = ((uint64_t)op[0] * 2785037492349527960UL) + ((uint64_t)op[1] * 16839937630611475622UL) + ((uint64_t)op[2] * 15499340231624087516UL) + ((uint64_t)op[3] * 12522298333479763018UL) + ((uint64_t)op[4] * 12113416962807399125UL) + ((uint64_t)op[5] * 10883422459652406185UL) + ((((uint64_t)op[6] * 5342609006528468788UL) + ((uint64_t)op[7] * 6739752439300693295UL) + ((uint64_t)op[8] * 17755411103475981349UL) + ((uint64_t)op[9] * 1196851098066101627UL) + ((uint64_t)op[10] * 5522947523484628917UL) + ((uint64_t)op[11] * 7399264221176564053UL)) * 14);
	tmp_q[6] = ((uint64_t)op[0] * 7399264221176564053UL) + ((uint64_t)op[1] * 2785037492349527960UL) + ((uint64_t)op[2] * 16839937630611475622UL) + ((uint64_t)op[3] * 15499340231624087516UL) + ((uint64_t)op[4] * 12522298333479763018UL) + ((uint64_t)op[5] * 12113416962807399125UL) + ((uint64_t)op[6] * 10883422459652406185UL) + ((((uint64_t)op[7] * 5342609006528468788UL) + ((uint64_t)op[8] * 6739752439300693295UL) + ((uint64_t)op[9] * 17755411103475981349UL) + ((uint64_t)op[10] * 1196851098066101627UL) + ((uint64_t)op[11] * 5522947523484628917UL)) * 14);
	tmp_q[7] = ((uint64_t)op[0] * 5522947523484628917UL) + ((uint64_t)op[1] * 7399264221176564053UL) + ((uint64_t)op[2] * 2785037492349527960UL) + ((uint64_t)op[3] * 16839937630611475622UL) + ((uint64_t)op[4] * 15499340231624087516UL) + ((uint64_t)op[5] * 12522298333479763018UL) + ((uint64_t)op[6] * 12113416962807399125UL) + ((uint64_t)op[7] * 10883422459652406185UL) + ((((uint64_t)op[8] * 5342609006528468788UL) + ((uint64_t)op[9] * 6739752439300693295UL) + ((uint64_t)op[10] * 17755411103475981349UL) + ((uint64_t)op[11] * 1196851098066101627UL)) * 14);
	tmp_q[8] = ((uint64_t)op[0] * 1196851098066101627UL) + ((uint64_t)op[1] * 5522947523484628917UL) + ((uint64_t)op[2] * 7399264221176564053UL) + ((uint64_t)op[3] * 2785037492349527960UL) + ((uint64_t)op[4] * 16839937630611475622UL) + ((uint64_t)op[5] * 15499340231624087516UL) + ((uint64_t)op[6] * 12522298333479763018UL) + ((uint64_t)op[7] * 12113416962807399125UL) + ((uint64_t)op[8] * 10883422459652406185UL) + ((((uint64_t)op[9] * 5342609006528468788UL) + ((uint64_t)op[10] * 6739752439300693295UL) + ((uint64_t)op[11] * 17755411103475981349UL)) * 14);
	tmp_q[9] = ((uint64_t)op[0] * 17755411103475981349UL) + ((uint64_t)op[1] * 1196851098066101627UL) + ((uint64_t)op[2] * 5522947523484628917UL) + ((uint64_t)op[3] * 7399264221176564053UL) + ((uint64_t)op[4] * 2785037492349527960UL) + ((uint64_t)op[5] * 16839937630611475622UL) + ((uint64_t)op[6] * 15499340231624087516UL) + ((uint64_t)op[7] * 12522298333479763018UL) + ((uint64_t)op[8] * 12113416962807399125UL) + ((uint64_t)op[9] * 10883422459652406185UL) + ((((uint64_t)op[10] * 5342609006528468788UL) + ((uint64_t)op[11] * 6739752439300693295UL)) * 14);
	tmp_q[10] = ((uint64_t)op[0] * 6739752439300693295UL) + ((uint64_t)op[1] * 17755411103475981349UL) + ((uint64_t)op[2] * 1196851098066101627UL) + ((uint64_t)op[3] * 5522947523484628917UL) + ((uint64_t)op[4] * 7399264221176564053UL) + ((uint64_t)op[5] * 2785037492349527960UL) + ((uint64_t)op[6] * 16839937630611475622UL) + ((uint64_t)op[7] * 15499340231624087516UL) + ((uint64_t)op[8] * 12522298333479763018UL) + ((uint64_t)op[9] * 12113416962807399125UL) + ((uint64_t)op[10] * 10883422459652406185UL) + ((uint64_t)op[11] * 1009549796560356568UL);
	tmp_q[11] = ((uint64_t)op[0] * 5342609006528468788UL) + ((uint64_t)op[1] * 6739752439300693295UL) + ((uint64_t)op[2] * 17755411103475981349UL) + ((uint64_t)op[3] * 1196851098066101627UL) + ((uint64_t)op[4] * 5522947523484628917UL) + ((uint64_t)op[5] * 7399264221176564053UL) + ((uint64_t)op[6] * 2785037492349527960UL) + ((uint64_t)op[7] * 16839937630611475622UL) + ((uint64_t)op[8] * 15499340231624087516UL) + ((uint64_t)op[9] * 12522298333479763018UL) + ((uint64_t)op[10] * 12113416962807399125UL) + ((uint64_t)op[11] * 10883422459652406185UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5785270094587L) + ((((int128)tmp_q[1] * 6192656033489L) + ((int128)tmp_q[2] * 4229410366826L) - ((int128)tmp_q[3] * 3347363180073L) - ((int128)tmp_q[4] * 1146016764131L) + ((int128)tmp_q[5] * 4109476706016L) - ((int128)tmp_q[6] * 802969910418L) + ((int128)tmp_q[7] * 2920139703035L) + ((int128)tmp_q[8] * 4585184481551L) - ((int128)tmp_q[9] * 2441601366583L) - ((int128)tmp_q[10] * 956743807681L) - ((int128)tmp_q[11] * 343334231135L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 343334231135L) - ((int128)tmp_q[1] * 5785270094587L) + ((((int128)tmp_q[2] * 6192656033489L) + ((int128)tmp_q[3] * 4229410366826L) - ((int128)tmp_q[4] * 3347363180073L) - ((int128)tmp_q[5] * 1146016764131L) + ((int128)tmp_q[6] * 4109476706016L) - ((int128)tmp_q[7] * 802969910418L) + ((int128)tmp_q[8] * 2920139703035L) + ((int128)tmp_q[9] * 4585184481551L) - ((int128)tmp_q[10] * 2441601366583L) - ((int128)tmp_q[11] * 956743807681L)) * 14);
	tmp_zero[2] = -((int128)tmp_q[0] * 956743807681L) - ((int128)tmp_q[1] * 343334231135L) - ((int128)tmp_q[2] * 5785270094587L) + ((((int128)tmp_q[3] * 6192656033489L) + ((int128)tmp_q[4] * 4229410366826L) - ((int128)tmp_q[5] * 3347363180073L) - ((int128)tmp_q[6] * 1146016764131L) + ((int128)tmp_q[7] * 4109476706016L) - ((int128)tmp_q[8] * 802969910418L) + ((int128)tmp_q[9] * 2920139703035L) + ((int128)tmp_q[10] * 4585184481551L) - ((int128)tmp_q[11] * 2441601366583L)) * 14);
	tmp_zero[3] = -((int128)tmp_q[0] * 2441601366583L) - ((int128)tmp_q[1] * 956743807681L) - ((int128)tmp_q[2] * 343334231135L) - ((int128)tmp_q[3] * 5785270094587L) + ((((int128)tmp_q[4] * 6192656033489L) + ((int128)tmp_q[5] * 4229410366826L) - ((int128)tmp_q[6] * 3347363180073L) - ((int128)tmp_q[7] * 1146016764131L) + ((int128)tmp_q[8] * 4109476706016L) - ((int128)tmp_q[9] * 802969910418L) + ((int128)tmp_q[10] * 2920139703035L) + ((int128)tmp_q[11] * 4585184481551L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 4585184481551L) - ((int128)tmp_q[1] * 2441601366583L) - ((int128)tmp_q[2] * 956743807681L) - ((int128)tmp_q[3] * 343334231135L) - ((int128)tmp_q[4] * 5785270094587L) + ((((int128)tmp_q[5] * 6192656033489L) + ((int128)tmp_q[6] * 4229410366826L) - ((int128)tmp_q[7] * 3347363180073L) - ((int128)tmp_q[8] * 1146016764131L) + ((int128)tmp_q[9] * 4109476706016L) - ((int128)tmp_q[10] * 802969910418L) + ((int128)tmp_q[11] * 2920139703035L)) * 14);
	tmp_zero[5] = ((int128)tmp_q[0] * 2920139703035L) + ((int128)tmp_q[1] * 4585184481551L) - ((int128)tmp_q[2] * 2441601366583L) - ((int128)tmp_q[3] * 956743807681L) - ((int128)tmp_q[4] * 343334231135L) - ((int128)tmp_q[5] * 5785270094587L) + ((((int128)tmp_q[6] * 6192656033489L) + ((int128)tmp_q[7] * 4229410366826L) - ((int128)tmp_q[8] * 3347363180073L) - ((int128)tmp_q[9] * 1146016764131L) + ((int128)tmp_q[10] * 4109476706016L) - ((int128)tmp_q[11] * 802969910418L)) * 14);
	tmp_zero[6] = -((int128)tmp_q[0] * 802969910418L) + ((int128)tmp_q[1] * 2920139703035L) + ((int128)tmp_q[2] * 4585184481551L) - ((int128)tmp_q[3] * 2441601366583L) - ((int128)tmp_q[4] * 956743807681L) - ((int128)tmp_q[5] * 343334231135L) - ((int128)tmp_q[6] * 5785270094587L) + ((((int128)tmp_q[7] * 6192656033489L) + ((int128)tmp_q[8] * 4229410366826L) - ((int128)tmp_q[9] * 3347363180073L) - ((int128)tmp_q[10] * 1146016764131L) + ((int128)tmp_q[11] * 4109476706016L)) * 14);
	tmp_zero[7] = ((int128)tmp_q[0] * 4109476706016L) - ((int128)tmp_q[1] * 802969910418L) + ((int128)tmp_q[2] * 2920139703035L) + ((int128)tmp_q[3] * 4585184481551L) - ((int128)tmp_q[4] * 2441601366583L) - ((int128)tmp_q[5] * 956743807681L) - ((int128)tmp_q[6] * 343334231135L) - ((int128)tmp_q[7] * 5785270094587L) + ((((int128)tmp_q[8] * 6192656033489L) + ((int128)tmp_q[9] * 4229410366826L) - ((int128)tmp_q[10] * 3347363180073L) - ((int128)tmp_q[11] * 1146016764131L)) * 14);
	tmp_zero[8] = -((int128)tmp_q[0] * 1146016764131L) + ((int128)tmp_q[1] * 4109476706016L) - ((int128)tmp_q[2] * 802969910418L) + ((int128)tmp_q[3] * 2920139703035L) + ((int128)tmp_q[4] * 4585184481551L) - ((int128)tmp_q[5] * 2441601366583L) - ((int128)tmp_q[6] * 956743807681L) - ((int128)tmp_q[7] * 343334231135L) - ((int128)tmp_q[8] * 5785270094587L) + ((((int128)tmp_q[9] * 6192656033489L) + ((int128)tmp_q[10] * 4229410366826L) - ((int128)tmp_q[11] * 3347363180073L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 3347363180073L) - ((int128)tmp_q[1] * 1146016764131L) + ((int128)tmp_q[2] * 4109476706016L) - ((int128)tmp_q[3] * 802969910418L) + ((int128)tmp_q[4] * 2920139703035L) + ((int128)tmp_q[5] * 4585184481551L) - ((int128)tmp_q[6] * 2441601366583L) - ((int128)tmp_q[7] * 956743807681L) - ((int128)tmp_q[8] * 343334231135L) - ((int128)tmp_q[9] * 5785270094587L) + ((((int128)tmp_q[10] * 6192656033489L) + ((int128)tmp_q[11] * 4229410366826L)) * 14);
	tmp_zero[10] = ((int128)tmp_q[0] * 4229410366826L) - ((int128)tmp_q[1] * 3347363180073L) - ((int128)tmp_q[2] * 1146016764131L) + ((int128)tmp_q[3] * 4109476706016L) - ((int128)tmp_q[4] * 802969910418L) + ((int128)tmp_q[5] * 2920139703035L) + ((int128)tmp_q[6] * 4585184481551L) - ((int128)tmp_q[7] * 2441601366583L) - ((int128)tmp_q[8] * 956743807681L) - ((int128)tmp_q[9] * 343334231135L) - ((int128)tmp_q[10] * 5785270094587L) + ((int128)tmp_q[11] * 86697184468846L);
	tmp_zero[11] = ((int128)tmp_q[0] * 6192656033489L) + ((int128)tmp_q[1] * 4229410366826L) - ((int128)tmp_q[2] * 3347363180073L) - ((int128)tmp_q[3] * 1146016764131L) + ((int128)tmp_q[4] * 4109476706016L) - ((int128)tmp_q[5] * 802969910418L) + ((int128)tmp_q[6] * 2920139703035L) + ((int128)tmp_q[7] * 4585184481551L) - ((int128)tmp_q[8] * 2441601366583L) - ((int128)tmp_q[9] * 956743807681L) - ((int128)tmp_q[10] * 343334231135L) - ((int128)tmp_q[11] * 5785270094587L);

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

