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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16262177450320844943UL) + ((((uint64_t)op[1] * 2670187825735741494UL) + ((uint64_t)op[2] * 8472561835038459760UL) + ((uint64_t)op[3] * 11663133898145033421UL) + ((uint64_t)op[4] * 7051428226973299194UL) + ((uint64_t)op[5] * 14196757435020080270UL) + ((uint64_t)op[6] * 14588882102184077017UL) + ((uint64_t)op[7] * 4639115463508342799UL) + ((uint64_t)op[8] * 13800353308576437592UL) + ((uint64_t)op[9] * 18383604863254359164UL) + ((uint64_t)op[10] * 784066946546315601UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 784066946546315601UL) + ((uint64_t)op[1] * 16262177450320844943UL) + ((((uint64_t)op[2] * 2670187825735741494UL) + ((uint64_t)op[3] * 8472561835038459760UL) + ((uint64_t)op[4] * 11663133898145033421UL) + ((uint64_t)op[5] * 7051428226973299194UL) + ((uint64_t)op[6] * 14196757435020080270UL) + ((uint64_t)op[7] * 14588882102184077017UL) + ((uint64_t)op[8] * 4639115463508342799UL) + ((uint64_t)op[9] * 13800353308576437592UL) + ((uint64_t)op[10] * 18383604863254359164UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 18383604863254359164UL) + ((uint64_t)op[1] * 784066946546315601UL) + ((uint64_t)op[2] * 16262177450320844943UL) + ((((uint64_t)op[3] * 2670187825735741494UL) + ((uint64_t)op[4] * 8472561835038459760UL) + ((uint64_t)op[5] * 11663133898145033421UL) + ((uint64_t)op[6] * 7051428226973299194UL) + ((uint64_t)op[7] * 14196757435020080270UL) + ((uint64_t)op[8] * 14588882102184077017UL) + ((uint64_t)op[9] * 4639115463508342799UL) + ((uint64_t)op[10] * 13800353308576437592UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 13800353308576437592UL) + ((uint64_t)op[1] * 18383604863254359164UL) + ((uint64_t)op[2] * 784066946546315601UL) + ((uint64_t)op[3] * 16262177450320844943UL) + ((((uint64_t)op[4] * 2670187825735741494UL) + ((uint64_t)op[5] * 8472561835038459760UL) + ((uint64_t)op[6] * 11663133898145033421UL) + ((uint64_t)op[7] * 7051428226973299194UL) + ((uint64_t)op[8] * 14196757435020080270UL) + ((uint64_t)op[9] * 14588882102184077017UL) + ((uint64_t)op[10] * 4639115463508342799UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 4639115463508342799UL) + ((uint64_t)op[1] * 13800353308576437592UL) + ((uint64_t)op[2] * 18383604863254359164UL) + ((uint64_t)op[3] * 784066946546315601UL) + ((uint64_t)op[4] * 16262177450320844943UL) + ((((uint64_t)op[5] * 2670187825735741494UL) + ((uint64_t)op[6] * 8472561835038459760UL) + ((uint64_t)op[7] * 11663133898145033421UL) + ((uint64_t)op[8] * 7051428226973299194UL) + ((uint64_t)op[9] * 14196757435020080270UL) + ((uint64_t)op[10] * 14588882102184077017UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 14588882102184077017UL) + ((uint64_t)op[1] * 4639115463508342799UL) + ((uint64_t)op[2] * 13800353308576437592UL) + ((uint64_t)op[3] * 18383604863254359164UL) + ((uint64_t)op[4] * 784066946546315601UL) + ((uint64_t)op[5] * 16262177450320844943UL) + ((((uint64_t)op[6] * 2670187825735741494UL) + ((uint64_t)op[7] * 8472561835038459760UL) + ((uint64_t)op[8] * 11663133898145033421UL) + ((uint64_t)op[9] * 7051428226973299194UL) + ((uint64_t)op[10] * 14196757435020080270UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 14196757435020080270UL) + ((uint64_t)op[1] * 14588882102184077017UL) + ((uint64_t)op[2] * 4639115463508342799UL) + ((uint64_t)op[3] * 13800353308576437592UL) + ((uint64_t)op[4] * 18383604863254359164UL) + ((uint64_t)op[5] * 784066946546315601UL) + ((uint64_t)op[6] * 16262177450320844943UL) + ((((uint64_t)op[7] * 2670187825735741494UL) + ((uint64_t)op[8] * 8472561835038459760UL) + ((uint64_t)op[9] * 11663133898145033421UL) + ((uint64_t)op[10] * 7051428226973299194UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 7051428226973299194UL) + ((uint64_t)op[1] * 14196757435020080270UL) + ((uint64_t)op[2] * 14588882102184077017UL) + ((uint64_t)op[3] * 4639115463508342799UL) + ((uint64_t)op[4] * 13800353308576437592UL) + ((uint64_t)op[5] * 18383604863254359164UL) + ((uint64_t)op[6] * 784066946546315601UL) + ((uint64_t)op[7] * 16262177450320844943UL) + ((((uint64_t)op[8] * 2670187825735741494UL) + ((uint64_t)op[9] * 8472561835038459760UL) + ((uint64_t)op[10] * 11663133898145033421UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 11663133898145033421UL) + ((uint64_t)op[1] * 7051428226973299194UL) + ((uint64_t)op[2] * 14196757435020080270UL) + ((uint64_t)op[3] * 14588882102184077017UL) + ((uint64_t)op[4] * 4639115463508342799UL) + ((uint64_t)op[5] * 13800353308576437592UL) + ((uint64_t)op[6] * 18383604863254359164UL) + ((uint64_t)op[7] * 784066946546315601UL) + ((uint64_t)op[8] * 16262177450320844943UL) + ((((uint64_t)op[9] * 2670187825735741494UL) + ((uint64_t)op[10] * 8472561835038459760UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 8472561835038459760UL) + ((uint64_t)op[1] * 11663133898145033421UL) + ((uint64_t)op[2] * 7051428226973299194UL) + ((uint64_t)op[3] * 14196757435020080270UL) + ((uint64_t)op[4] * 14588882102184077017UL) + ((uint64_t)op[5] * 4639115463508342799UL) + ((uint64_t)op[6] * 13800353308576437592UL) + ((uint64_t)op[7] * 18383604863254359164UL) + ((uint64_t)op[8] * 784066946546315601UL) + ((uint64_t)op[9] * 16262177450320844943UL) + ((uint64_t)op[10] * 5340375651471482988UL);
	tmp_q[10] = ((uint64_t)op[0] * 2670187825735741494UL) + ((uint64_t)op[1] * 8472561835038459760UL) + ((uint64_t)op[2] * 11663133898145033421UL) + ((uint64_t)op[3] * 7051428226973299194UL) + ((uint64_t)op[4] * 14196757435020080270UL) + ((uint64_t)op[5] * 14588882102184077017UL) + ((uint64_t)op[6] * 4639115463508342799UL) + ((uint64_t)op[7] * 13800353308576437592UL) + ((uint64_t)op[8] * 18383604863254359164UL) + ((uint64_t)op[9] * 784066946546315601UL) + ((uint64_t)op[10] * 16262177450320844943UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4193342897775L) + ((((int128)tmp_q[1] * 46717218411936L) - ((int128)tmp_q[2] * 2852664485051L) - ((int128)tmp_q[3] * 64636343462574L) + ((int128)tmp_q[4] * 91860594163784L) - ((int128)tmp_q[5] * 17949733163902L) + ((int128)tmp_q[6] * 92786237768312L) + ((int128)tmp_q[7] * 4416472503452L) - ((int128)tmp_q[8] * 12491769983517L) - ((int128)tmp_q[9] * 76410259501951L) + ((int128)tmp_q[10] * 17540295409575L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 17540295409575L) + ((int128)tmp_q[1] * 4193342897775L) + ((((int128)tmp_q[2] * 46717218411936L) - ((int128)tmp_q[3] * 2852664485051L) - ((int128)tmp_q[4] * 64636343462574L) + ((int128)tmp_q[5] * 91860594163784L) - ((int128)tmp_q[6] * 17949733163902L) + ((int128)tmp_q[7] * 92786237768312L) + ((int128)tmp_q[8] * 4416472503452L) - ((int128)tmp_q[9] * 12491769983517L) - ((int128)tmp_q[10] * 76410259501951L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 76410259501951L) + ((int128)tmp_q[1] * 17540295409575L) + ((int128)tmp_q[2] * 4193342897775L) + ((((int128)tmp_q[3] * 46717218411936L) - ((int128)tmp_q[4] * 2852664485051L) - ((int128)tmp_q[5] * 64636343462574L) + ((int128)tmp_q[6] * 91860594163784L) - ((int128)tmp_q[7] * 17949733163902L) + ((int128)tmp_q[8] * 92786237768312L) + ((int128)tmp_q[9] * 4416472503452L) - ((int128)tmp_q[10] * 12491769983517L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 12491769983517L) - ((int128)tmp_q[1] * 76410259501951L) + ((int128)tmp_q[2] * 17540295409575L) + ((int128)tmp_q[3] * 4193342897775L) + ((((int128)tmp_q[4] * 46717218411936L) - ((int128)tmp_q[5] * 2852664485051L) - ((int128)tmp_q[6] * 64636343462574L) + ((int128)tmp_q[7] * 91860594163784L) - ((int128)tmp_q[8] * 17949733163902L) + ((int128)tmp_q[9] * 92786237768312L) + ((int128)tmp_q[10] * 4416472503452L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 4416472503452L) - ((int128)tmp_q[1] * 12491769983517L) - ((int128)tmp_q[2] * 76410259501951L) + ((int128)tmp_q[3] * 17540295409575L) + ((int128)tmp_q[4] * 4193342897775L) + ((((int128)tmp_q[5] * 46717218411936L) - ((int128)tmp_q[6] * 2852664485051L) - ((int128)tmp_q[7] * 64636343462574L) + ((int128)tmp_q[8] * 91860594163784L) - ((int128)tmp_q[9] * 17949733163902L) + ((int128)tmp_q[10] * 92786237768312L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 92786237768312L) + ((int128)tmp_q[1] * 4416472503452L) - ((int128)tmp_q[2] * 12491769983517L) - ((int128)tmp_q[3] * 76410259501951L) + ((int128)tmp_q[4] * 17540295409575L) + ((int128)tmp_q[5] * 4193342897775L) + ((((int128)tmp_q[6] * 46717218411936L) - ((int128)tmp_q[7] * 2852664485051L) - ((int128)tmp_q[8] * 64636343462574L) + ((int128)tmp_q[9] * 91860594163784L) - ((int128)tmp_q[10] * 17949733163902L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 17949733163902L) + ((int128)tmp_q[1] * 92786237768312L) + ((int128)tmp_q[2] * 4416472503452L) - ((int128)tmp_q[3] * 12491769983517L) - ((int128)tmp_q[4] * 76410259501951L) + ((int128)tmp_q[5] * 17540295409575L) + ((int128)tmp_q[6] * 4193342897775L) + ((((int128)tmp_q[7] * 46717218411936L) - ((int128)tmp_q[8] * 2852664485051L) - ((int128)tmp_q[9] * 64636343462574L) + ((int128)tmp_q[10] * 91860594163784L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 91860594163784L) - ((int128)tmp_q[1] * 17949733163902L) + ((int128)tmp_q[2] * 92786237768312L) + ((int128)tmp_q[3] * 4416472503452L) - ((int128)tmp_q[4] * 12491769983517L) - ((int128)tmp_q[5] * 76410259501951L) + ((int128)tmp_q[6] * 17540295409575L) + ((int128)tmp_q[7] * 4193342897775L) + ((((int128)tmp_q[8] * 46717218411936L) - ((int128)tmp_q[9] * 2852664485051L) - ((int128)tmp_q[10] * 64636343462574L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 64636343462574L) + ((int128)tmp_q[1] * 91860594163784L) - ((int128)tmp_q[2] * 17949733163902L) + ((int128)tmp_q[3] * 92786237768312L) + ((int128)tmp_q[4] * 4416472503452L) - ((int128)tmp_q[5] * 12491769983517L) - ((int128)tmp_q[6] * 76410259501951L) + ((int128)tmp_q[7] * 17540295409575L) + ((int128)tmp_q[8] * 4193342897775L) + ((((int128)tmp_q[9] * 46717218411936L) - ((int128)tmp_q[10] * 2852664485051L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 2852664485051L) - ((int128)tmp_q[1] * 64636343462574L) + ((int128)tmp_q[2] * 91860594163784L) - ((int128)tmp_q[3] * 17949733163902L) + ((int128)tmp_q[4] * 92786237768312L) + ((int128)tmp_q[5] * 4416472503452L) - ((int128)tmp_q[6] * 12491769983517L) - ((int128)tmp_q[7] * 76410259501951L) + ((int128)tmp_q[8] * 17540295409575L) + ((int128)tmp_q[9] * 4193342897775L) + ((int128)tmp_q[10] * 93434436823872L);
	tmp_zero[10] = ((int128)tmp_q[0] * 46717218411936L) - ((int128)tmp_q[1] * 2852664485051L) - ((int128)tmp_q[2] * 64636343462574L) + ((int128)tmp_q[3] * 91860594163784L) - ((int128)tmp_q[4] * 17949733163902L) + ((int128)tmp_q[5] * 92786237768312L) + ((int128)tmp_q[6] * 4416472503452L) - ((int128)tmp_q[7] * 12491769983517L) - ((int128)tmp_q[8] * 76410259501951L) + ((int128)tmp_q[9] * 17540295409575L) + ((int128)tmp_q[10] * 4193342897775L);

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

