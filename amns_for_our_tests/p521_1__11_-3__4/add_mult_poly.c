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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12290335185173767623UL) + ((((uint64_t)op[1] * 7109802905795275645UL) + ((uint64_t)op[2] * 9560756748099895475UL) + ((uint64_t)op[3] * 7320690403619197660UL) + ((uint64_t)op[4] * 4208506182406845564UL) + ((uint64_t)op[5] * 9232796789409382684UL) + ((uint64_t)op[6] * 6566418564624725768UL) + ((uint64_t)op[7] * 18023133085060502180UL) + ((uint64_t)op[8] * 11094896484681578238UL) + ((uint64_t)op[9] * 8451117314884735307UL) + ((uint64_t)op[10] * 5654787509874330299UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 5654787509874330299UL) + ((uint64_t)op[1] * 12290335185173767623UL) + ((((uint64_t)op[2] * 7109802905795275645UL) + ((uint64_t)op[3] * 9560756748099895475UL) + ((uint64_t)op[4] * 7320690403619197660UL) + ((uint64_t)op[5] * 4208506182406845564UL) + ((uint64_t)op[6] * 9232796789409382684UL) + ((uint64_t)op[7] * 6566418564624725768UL) + ((uint64_t)op[8] * 18023133085060502180UL) + ((uint64_t)op[9] * 11094896484681578238UL) + ((uint64_t)op[10] * 8451117314884735307UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 8451117314884735307UL) + ((uint64_t)op[1] * 5654787509874330299UL) + ((uint64_t)op[2] * 12290335185173767623UL) + ((((uint64_t)op[3] * 7109802905795275645UL) + ((uint64_t)op[4] * 9560756748099895475UL) + ((uint64_t)op[5] * 7320690403619197660UL) + ((uint64_t)op[6] * 4208506182406845564UL) + ((uint64_t)op[7] * 9232796789409382684UL) + ((uint64_t)op[8] * 6566418564624725768UL) + ((uint64_t)op[9] * 18023133085060502180UL) + ((uint64_t)op[10] * 11094896484681578238UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 11094896484681578238UL) + ((uint64_t)op[1] * 8451117314884735307UL) + ((uint64_t)op[2] * 5654787509874330299UL) + ((uint64_t)op[3] * 12290335185173767623UL) + ((((uint64_t)op[4] * 7109802905795275645UL) + ((uint64_t)op[5] * 9560756748099895475UL) + ((uint64_t)op[6] * 7320690403619197660UL) + ((uint64_t)op[7] * 4208506182406845564UL) + ((uint64_t)op[8] * 9232796789409382684UL) + ((uint64_t)op[9] * 6566418564624725768UL) + ((uint64_t)op[10] * 18023133085060502180UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 18023133085060502180UL) + ((uint64_t)op[1] * 11094896484681578238UL) + ((uint64_t)op[2] * 8451117314884735307UL) + ((uint64_t)op[3] * 5654787509874330299UL) + ((uint64_t)op[4] * 12290335185173767623UL) + ((((uint64_t)op[5] * 7109802905795275645UL) + ((uint64_t)op[6] * 9560756748099895475UL) + ((uint64_t)op[7] * 7320690403619197660UL) + ((uint64_t)op[8] * 4208506182406845564UL) + ((uint64_t)op[9] * 9232796789409382684UL) + ((uint64_t)op[10] * 6566418564624725768UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 6566418564624725768UL) + ((uint64_t)op[1] * 18023133085060502180UL) + ((uint64_t)op[2] * 11094896484681578238UL) + ((uint64_t)op[3] * 8451117314884735307UL) + ((uint64_t)op[4] * 5654787509874330299UL) + ((uint64_t)op[5] * 12290335185173767623UL) + ((((uint64_t)op[6] * 7109802905795275645UL) + ((uint64_t)op[7] * 9560756748099895475UL) + ((uint64_t)op[8] * 7320690403619197660UL) + ((uint64_t)op[9] * 4208506182406845564UL) + ((uint64_t)op[10] * 9232796789409382684UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 9232796789409382684UL) + ((uint64_t)op[1] * 6566418564624725768UL) + ((uint64_t)op[2] * 18023133085060502180UL) + ((uint64_t)op[3] * 11094896484681578238UL) + ((uint64_t)op[4] * 8451117314884735307UL) + ((uint64_t)op[5] * 5654787509874330299UL) + ((uint64_t)op[6] * 12290335185173767623UL) + ((((uint64_t)op[7] * 7109802905795275645UL) + ((uint64_t)op[8] * 9560756748099895475UL) + ((uint64_t)op[9] * 7320690403619197660UL) + ((uint64_t)op[10] * 4208506182406845564UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 4208506182406845564UL) + ((uint64_t)op[1] * 9232796789409382684UL) + ((uint64_t)op[2] * 6566418564624725768UL) + ((uint64_t)op[3] * 18023133085060502180UL) + ((uint64_t)op[4] * 11094896484681578238UL) + ((uint64_t)op[5] * 8451117314884735307UL) + ((uint64_t)op[6] * 5654787509874330299UL) + ((uint64_t)op[7] * 12290335185173767623UL) + ((((uint64_t)op[8] * 7109802905795275645UL) + ((uint64_t)op[9] * 9560756748099895475UL) + ((uint64_t)op[10] * 7320690403619197660UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 7320690403619197660UL) + ((uint64_t)op[1] * 4208506182406845564UL) + ((uint64_t)op[2] * 9232796789409382684UL) + ((uint64_t)op[3] * 6566418564624725768UL) + ((uint64_t)op[4] * 18023133085060502180UL) + ((uint64_t)op[5] * 11094896484681578238UL) + ((uint64_t)op[6] * 8451117314884735307UL) + ((uint64_t)op[7] * 5654787509874330299UL) + ((uint64_t)op[8] * 12290335185173767623UL) + ((((uint64_t)op[9] * 7109802905795275645UL) + ((uint64_t)op[10] * 9560756748099895475UL)) * 18446744073709551613);
	tmp_q[9] = ((uint64_t)op[0] * 9560756748099895475UL) + ((uint64_t)op[1] * 7320690403619197660UL) + ((uint64_t)op[2] * 4208506182406845564UL) + ((uint64_t)op[3] * 9232796789409382684UL) + ((uint64_t)op[4] * 6566418564624725768UL) + ((uint64_t)op[5] * 18023133085060502180UL) + ((uint64_t)op[6] * 11094896484681578238UL) + ((uint64_t)op[7] * 8451117314884735307UL) + ((uint64_t)op[8] * 5654787509874330299UL) + ((uint64_t)op[9] * 12290335185173767623UL) + ((uint64_t)op[10] * 15564079430033276297UL);
	tmp_q[10] = ((uint64_t)op[0] * 7109802905795275645UL) + ((uint64_t)op[1] * 9560756748099895475UL) + ((uint64_t)op[2] * 7320690403619197660UL) + ((uint64_t)op[3] * 4208506182406845564UL) + ((uint64_t)op[4] * 9232796789409382684UL) + ((uint64_t)op[5] * 6566418564624725768UL) + ((uint64_t)op[6] * 18023133085060502180UL) + ((uint64_t)op[7] * 11094896484681578238UL) + ((uint64_t)op[8] * 8451117314884735307UL) + ((uint64_t)op[9] * 5654787509874330299UL) + ((uint64_t)op[10] * 12290335185173767623UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 58085932624957L) - ((-((int128)tmp_q[1] * 56293702956467L) + ((int128)tmp_q[2] * 49670302403383L) + ((int128)tmp_q[3] * 69384316751074L) - ((int128)tmp_q[4] * 4350872300071L) + ((int128)tmp_q[5] * 47320695897401L) + ((int128)tmp_q[6] * 91984323060967L) - ((int128)tmp_q[7] * 61902444103497L) - ((int128)tmp_q[8] * 89670818449032L) - ((int128)tmp_q[9] * 38498057092931L) - ((int128)tmp_q[10] * 53938779284827L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 53938779284827L) + ((int128)tmp_q[1] * 58085932624957L) - ((-((int128)tmp_q[2] * 56293702956467L) + ((int128)tmp_q[3] * 49670302403383L) + ((int128)tmp_q[4] * 69384316751074L) - ((int128)tmp_q[5] * 4350872300071L) + ((int128)tmp_q[6] * 47320695897401L) + ((int128)tmp_q[7] * 91984323060967L) - ((int128)tmp_q[8] * 61902444103497L) - ((int128)tmp_q[9] * 89670818449032L) - ((int128)tmp_q[10] * 38498057092931L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 38498057092931L) - ((int128)tmp_q[1] * 53938779284827L) + ((int128)tmp_q[2] * 58085932624957L) - ((-((int128)tmp_q[3] * 56293702956467L) + ((int128)tmp_q[4] * 49670302403383L) + ((int128)tmp_q[5] * 69384316751074L) - ((int128)tmp_q[6] * 4350872300071L) + ((int128)tmp_q[7] * 47320695897401L) + ((int128)tmp_q[8] * 91984323060967L) - ((int128)tmp_q[9] * 61902444103497L) - ((int128)tmp_q[10] * 89670818449032L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 89670818449032L) - ((int128)tmp_q[1] * 38498057092931L) - ((int128)tmp_q[2] * 53938779284827L) + ((int128)tmp_q[3] * 58085932624957L) - ((-((int128)tmp_q[4] * 56293702956467L) + ((int128)tmp_q[5] * 49670302403383L) + ((int128)tmp_q[6] * 69384316751074L) - ((int128)tmp_q[7] * 4350872300071L) + ((int128)tmp_q[8] * 47320695897401L) + ((int128)tmp_q[9] * 91984323060967L) - ((int128)tmp_q[10] * 61902444103497L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 61902444103497L) - ((int128)tmp_q[1] * 89670818449032L) - ((int128)tmp_q[2] * 38498057092931L) - ((int128)tmp_q[3] * 53938779284827L) + ((int128)tmp_q[4] * 58085932624957L) - ((-((int128)tmp_q[5] * 56293702956467L) + ((int128)tmp_q[6] * 49670302403383L) + ((int128)tmp_q[7] * 69384316751074L) - ((int128)tmp_q[8] * 4350872300071L) + ((int128)tmp_q[9] * 47320695897401L) + ((int128)tmp_q[10] * 91984323060967L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 91984323060967L) - ((int128)tmp_q[1] * 61902444103497L) - ((int128)tmp_q[2] * 89670818449032L) - ((int128)tmp_q[3] * 38498057092931L) - ((int128)tmp_q[4] * 53938779284827L) + ((int128)tmp_q[5] * 58085932624957L) - ((-((int128)tmp_q[6] * 56293702956467L) + ((int128)tmp_q[7] * 49670302403383L) + ((int128)tmp_q[8] * 69384316751074L) - ((int128)tmp_q[9] * 4350872300071L) + ((int128)tmp_q[10] * 47320695897401L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 47320695897401L) + ((int128)tmp_q[1] * 91984323060967L) - ((int128)tmp_q[2] * 61902444103497L) - ((int128)tmp_q[3] * 89670818449032L) - ((int128)tmp_q[4] * 38498057092931L) - ((int128)tmp_q[5] * 53938779284827L) + ((int128)tmp_q[6] * 58085932624957L) - ((-((int128)tmp_q[7] * 56293702956467L) + ((int128)tmp_q[8] * 49670302403383L) + ((int128)tmp_q[9] * 69384316751074L) - ((int128)tmp_q[10] * 4350872300071L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 4350872300071L) + ((int128)tmp_q[1] * 47320695897401L) + ((int128)tmp_q[2] * 91984323060967L) - ((int128)tmp_q[3] * 61902444103497L) - ((int128)tmp_q[4] * 89670818449032L) - ((int128)tmp_q[5] * 38498057092931L) - ((int128)tmp_q[6] * 53938779284827L) + ((int128)tmp_q[7] * 58085932624957L) - ((-((int128)tmp_q[8] * 56293702956467L) + ((int128)tmp_q[9] * 49670302403383L) + ((int128)tmp_q[10] * 69384316751074L)) * 3);
	tmp_zero[8] = ((int128)tmp_q[0] * 69384316751074L) - ((int128)tmp_q[1] * 4350872300071L) + ((int128)tmp_q[2] * 47320695897401L) + ((int128)tmp_q[3] * 91984323060967L) - ((int128)tmp_q[4] * 61902444103497L) - ((int128)tmp_q[5] * 89670818449032L) - ((int128)tmp_q[6] * 38498057092931L) - ((int128)tmp_q[7] * 53938779284827L) + ((int128)tmp_q[8] * 58085932624957L) - ((-((int128)tmp_q[9] * 56293702956467L) + ((int128)tmp_q[10] * 49670302403383L)) * 3);
	tmp_zero[9] = ((int128)tmp_q[0] * 49670302403383L) + ((int128)tmp_q[1] * 69384316751074L) - ((int128)tmp_q[2] * 4350872300071L) + ((int128)tmp_q[3] * 47320695897401L) + ((int128)tmp_q[4] * 91984323060967L) - ((int128)tmp_q[5] * 61902444103497L) - ((int128)tmp_q[6] * 89670818449032L) - ((int128)tmp_q[7] * 38498057092931L) - ((int128)tmp_q[8] * 53938779284827L) + ((int128)tmp_q[9] * 58085932624957L) + ((int128)tmp_q[10] * 168881108869401L);
	tmp_zero[10] = -((int128)tmp_q[0] * 56293702956467L) + ((int128)tmp_q[1] * 49670302403383L) + ((int128)tmp_q[2] * 69384316751074L) - ((int128)tmp_q[3] * 4350872300071L) + ((int128)tmp_q[4] * 47320695897401L) + ((int128)tmp_q[5] * 91984323060967L) - ((int128)tmp_q[6] * 61902444103497L) - ((int128)tmp_q[7] * 89670818449032L) - ((int128)tmp_q[8] * 38498057092931L) - ((int128)tmp_q[9] * 53938779284827L) + ((int128)tmp_q[10] * 58085932624957L);

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

