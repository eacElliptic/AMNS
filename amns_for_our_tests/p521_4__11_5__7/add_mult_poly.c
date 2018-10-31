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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9421290090258827841UL) + ((((uint64_t)op[1] * 5986585843038591886UL) + ((uint64_t)op[2] * 14213659503220740940UL) + ((uint64_t)op[3] * 16982430686005027885UL) + ((uint64_t)op[4] * 14970932391691794141UL) + ((uint64_t)op[5] * 14033948716194374925UL) + ((uint64_t)op[6] * 10677038138303389335UL) + ((uint64_t)op[7] * 4298178193398895357UL) + ((uint64_t)op[8] * 459103241890573720UL) + ((uint64_t)op[9] * 420674507042013611UL) + ((uint64_t)op[10] * 7729600067392707748UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 7729600067392707748UL) + ((uint64_t)op[1] * 9421290090258827841UL) + ((((uint64_t)op[2] * 5986585843038591886UL) + ((uint64_t)op[3] * 14213659503220740940UL) + ((uint64_t)op[4] * 16982430686005027885UL) + ((uint64_t)op[5] * 14970932391691794141UL) + ((uint64_t)op[6] * 14033948716194374925UL) + ((uint64_t)op[7] * 10677038138303389335UL) + ((uint64_t)op[8] * 4298178193398895357UL) + ((uint64_t)op[9] * 459103241890573720UL) + ((uint64_t)op[10] * 420674507042013611UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 420674507042013611UL) + ((uint64_t)op[1] * 7729600067392707748UL) + ((uint64_t)op[2] * 9421290090258827841UL) + ((((uint64_t)op[3] * 5986585843038591886UL) + ((uint64_t)op[4] * 14213659503220740940UL) + ((uint64_t)op[5] * 16982430686005027885UL) + ((uint64_t)op[6] * 14970932391691794141UL) + ((uint64_t)op[7] * 14033948716194374925UL) + ((uint64_t)op[8] * 10677038138303389335UL) + ((uint64_t)op[9] * 4298178193398895357UL) + ((uint64_t)op[10] * 459103241890573720UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 459103241890573720UL) + ((uint64_t)op[1] * 420674507042013611UL) + ((uint64_t)op[2] * 7729600067392707748UL) + ((uint64_t)op[3] * 9421290090258827841UL) + ((((uint64_t)op[4] * 5986585843038591886UL) + ((uint64_t)op[5] * 14213659503220740940UL) + ((uint64_t)op[6] * 16982430686005027885UL) + ((uint64_t)op[7] * 14970932391691794141UL) + ((uint64_t)op[8] * 14033948716194374925UL) + ((uint64_t)op[9] * 10677038138303389335UL) + ((uint64_t)op[10] * 4298178193398895357UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 4298178193398895357UL) + ((uint64_t)op[1] * 459103241890573720UL) + ((uint64_t)op[2] * 420674507042013611UL) + ((uint64_t)op[3] * 7729600067392707748UL) + ((uint64_t)op[4] * 9421290090258827841UL) + ((((uint64_t)op[5] * 5986585843038591886UL) + ((uint64_t)op[6] * 14213659503220740940UL) + ((uint64_t)op[7] * 16982430686005027885UL) + ((uint64_t)op[8] * 14970932391691794141UL) + ((uint64_t)op[9] * 14033948716194374925UL) + ((uint64_t)op[10] * 10677038138303389335UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 10677038138303389335UL) + ((uint64_t)op[1] * 4298178193398895357UL) + ((uint64_t)op[2] * 459103241890573720UL) + ((uint64_t)op[3] * 420674507042013611UL) + ((uint64_t)op[4] * 7729600067392707748UL) + ((uint64_t)op[5] * 9421290090258827841UL) + ((((uint64_t)op[6] * 5986585843038591886UL) + ((uint64_t)op[7] * 14213659503220740940UL) + ((uint64_t)op[8] * 16982430686005027885UL) + ((uint64_t)op[9] * 14970932391691794141UL) + ((uint64_t)op[10] * 14033948716194374925UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 14033948716194374925UL) + ((uint64_t)op[1] * 10677038138303389335UL) + ((uint64_t)op[2] * 4298178193398895357UL) + ((uint64_t)op[3] * 459103241890573720UL) + ((uint64_t)op[4] * 420674507042013611UL) + ((uint64_t)op[5] * 7729600067392707748UL) + ((uint64_t)op[6] * 9421290090258827841UL) + ((((uint64_t)op[7] * 5986585843038591886UL) + ((uint64_t)op[8] * 14213659503220740940UL) + ((uint64_t)op[9] * 16982430686005027885UL) + ((uint64_t)op[10] * 14970932391691794141UL)) * 5);
	tmp_q[7] = ((uint64_t)op[0] * 14970932391691794141UL) + ((uint64_t)op[1] * 14033948716194374925UL) + ((uint64_t)op[2] * 10677038138303389335UL) + ((uint64_t)op[3] * 4298178193398895357UL) + ((uint64_t)op[4] * 459103241890573720UL) + ((uint64_t)op[5] * 420674507042013611UL) + ((uint64_t)op[6] * 7729600067392707748UL) + ((uint64_t)op[7] * 9421290090258827841UL) + ((((uint64_t)op[8] * 5986585843038591886UL) + ((uint64_t)op[9] * 14213659503220740940UL) + ((uint64_t)op[10] * 16982430686005027885UL)) * 5);
	tmp_q[8] = ((uint64_t)op[0] * 16982430686005027885UL) + ((uint64_t)op[1] * 14970932391691794141UL) + ((uint64_t)op[2] * 14033948716194374925UL) + ((uint64_t)op[3] * 10677038138303389335UL) + ((uint64_t)op[4] * 4298178193398895357UL) + ((uint64_t)op[5] * 459103241890573720UL) + ((uint64_t)op[6] * 420674507042013611UL) + ((uint64_t)op[7] * 7729600067392707748UL) + ((uint64_t)op[8] * 9421290090258827841UL) + ((((uint64_t)op[9] * 5986585843038591886UL) + ((uint64_t)op[10] * 14213659503220740940UL)) * 5);
	tmp_q[9] = ((uint64_t)op[0] * 14213659503220740940UL) + ((uint64_t)op[1] * 16982430686005027885UL) + ((uint64_t)op[2] * 14970932391691794141UL) + ((uint64_t)op[3] * 14033948716194374925UL) + ((uint64_t)op[4] * 10677038138303389335UL) + ((uint64_t)op[5] * 4298178193398895357UL) + ((uint64_t)op[6] * 459103241890573720UL) + ((uint64_t)op[7] * 420674507042013611UL) + ((uint64_t)op[8] * 7729600067392707748UL) + ((uint64_t)op[9] * 9421290090258827841UL) + ((uint64_t)op[10] * 11486185141483407814UL);
	tmp_q[10] = ((uint64_t)op[0] * 5986585843038591886UL) + ((uint64_t)op[1] * 14213659503220740940UL) + ((uint64_t)op[2] * 16982430686005027885UL) + ((uint64_t)op[3] * 14970932391691794141UL) + ((uint64_t)op[4] * 14033948716194374925UL) + ((uint64_t)op[5] * 10677038138303389335UL) + ((uint64_t)op[6] * 4298178193398895357UL) + ((uint64_t)op[7] * 459103241890573720UL) + ((uint64_t)op[8] * 420674507042013611UL) + ((uint64_t)op[9] * 7729600067392707748UL) + ((uint64_t)op[10] * 9421290090258827841UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 52504023431909L) + ((-((int128)tmp_q[1] * 53180276692910L) - ((int128)tmp_q[2] * 74911817787848L) + ((int128)tmp_q[3] * 62277769917104L) + ((int128)tmp_q[4] * 26175613713206L) + ((int128)tmp_q[5] * 2215605435586L) + ((int128)tmp_q[6] * 48930471118107L) + ((int128)tmp_q[7] * 16969694708707L) + ((int128)tmp_q[8] * 11025890776264L) + ((int128)tmp_q[9] * 43347723845914L) - ((int128)tmp_q[10] * 49324602339088L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 49324602339088L) + ((int128)tmp_q[1] * 52504023431909L) + ((-((int128)tmp_q[2] * 53180276692910L) - ((int128)tmp_q[3] * 74911817787848L) + ((int128)tmp_q[4] * 62277769917104L) + ((int128)tmp_q[5] * 26175613713206L) + ((int128)tmp_q[6] * 2215605435586L) + ((int128)tmp_q[7] * 48930471118107L) + ((int128)tmp_q[8] * 16969694708707L) + ((int128)tmp_q[9] * 11025890776264L) + ((int128)tmp_q[10] * 43347723845914L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 43347723845914L) - ((int128)tmp_q[1] * 49324602339088L) + ((int128)tmp_q[2] * 52504023431909L) + ((-((int128)tmp_q[3] * 53180276692910L) - ((int128)tmp_q[4] * 74911817787848L) + ((int128)tmp_q[5] * 62277769917104L) + ((int128)tmp_q[6] * 26175613713206L) + ((int128)tmp_q[7] * 2215605435586L) + ((int128)tmp_q[8] * 48930471118107L) + ((int128)tmp_q[9] * 16969694708707L) + ((int128)tmp_q[10] * 11025890776264L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 11025890776264L) + ((int128)tmp_q[1] * 43347723845914L) - ((int128)tmp_q[2] * 49324602339088L) + ((int128)tmp_q[3] * 52504023431909L) + ((-((int128)tmp_q[4] * 53180276692910L) - ((int128)tmp_q[5] * 74911817787848L) + ((int128)tmp_q[6] * 62277769917104L) + ((int128)tmp_q[7] * 26175613713206L) + ((int128)tmp_q[8] * 2215605435586L) + ((int128)tmp_q[9] * 48930471118107L) + ((int128)tmp_q[10] * 16969694708707L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 16969694708707L) + ((int128)tmp_q[1] * 11025890776264L) + ((int128)tmp_q[2] * 43347723845914L) - ((int128)tmp_q[3] * 49324602339088L) + ((int128)tmp_q[4] * 52504023431909L) + ((-((int128)tmp_q[5] * 53180276692910L) - ((int128)tmp_q[6] * 74911817787848L) + ((int128)tmp_q[7] * 62277769917104L) + ((int128)tmp_q[8] * 26175613713206L) + ((int128)tmp_q[9] * 2215605435586L) + ((int128)tmp_q[10] * 48930471118107L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 48930471118107L) + ((int128)tmp_q[1] * 16969694708707L) + ((int128)tmp_q[2] * 11025890776264L) + ((int128)tmp_q[3] * 43347723845914L) - ((int128)tmp_q[4] * 49324602339088L) + ((int128)tmp_q[5] * 52504023431909L) + ((-((int128)tmp_q[6] * 53180276692910L) - ((int128)tmp_q[7] * 74911817787848L) + ((int128)tmp_q[8] * 62277769917104L) + ((int128)tmp_q[9] * 26175613713206L) + ((int128)tmp_q[10] * 2215605435586L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 2215605435586L) + ((int128)tmp_q[1] * 48930471118107L) + ((int128)tmp_q[2] * 16969694708707L) + ((int128)tmp_q[3] * 11025890776264L) + ((int128)tmp_q[4] * 43347723845914L) - ((int128)tmp_q[5] * 49324602339088L) + ((int128)tmp_q[6] * 52504023431909L) + ((-((int128)tmp_q[7] * 53180276692910L) - ((int128)tmp_q[8] * 74911817787848L) + ((int128)tmp_q[9] * 62277769917104L) + ((int128)tmp_q[10] * 26175613713206L)) * 5);
	tmp_zero[7] = ((int128)tmp_q[0] * 26175613713206L) + ((int128)tmp_q[1] * 2215605435586L) + ((int128)tmp_q[2] * 48930471118107L) + ((int128)tmp_q[3] * 16969694708707L) + ((int128)tmp_q[4] * 11025890776264L) + ((int128)tmp_q[5] * 43347723845914L) - ((int128)tmp_q[6] * 49324602339088L) + ((int128)tmp_q[7] * 52504023431909L) + ((-((int128)tmp_q[8] * 53180276692910L) - ((int128)tmp_q[9] * 74911817787848L) + ((int128)tmp_q[10] * 62277769917104L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 62277769917104L) + ((int128)tmp_q[1] * 26175613713206L) + ((int128)tmp_q[2] * 2215605435586L) + ((int128)tmp_q[3] * 48930471118107L) + ((int128)tmp_q[4] * 16969694708707L) + ((int128)tmp_q[5] * 11025890776264L) + ((int128)tmp_q[6] * 43347723845914L) - ((int128)tmp_q[7] * 49324602339088L) + ((int128)tmp_q[8] * 52504023431909L) + ((-((int128)tmp_q[9] * 53180276692910L) - ((int128)tmp_q[10] * 74911817787848L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 74911817787848L) + ((int128)tmp_q[1] * 62277769917104L) + ((int128)tmp_q[2] * 26175613713206L) + ((int128)tmp_q[3] * 2215605435586L) + ((int128)tmp_q[4] * 48930471118107L) + ((int128)tmp_q[5] * 16969694708707L) + ((int128)tmp_q[6] * 11025890776264L) + ((int128)tmp_q[7] * 43347723845914L) - ((int128)tmp_q[8] * 49324602339088L) + ((int128)tmp_q[9] * 52504023431909L) - ((int128)tmp_q[10] * 265901383464550L);
	tmp_zero[10] = -((int128)tmp_q[0] * 53180276692910L) - ((int128)tmp_q[1] * 74911817787848L) + ((int128)tmp_q[2] * 62277769917104L) + ((int128)tmp_q[3] * 26175613713206L) + ((int128)tmp_q[4] * 2215605435586L) + ((int128)tmp_q[5] * 48930471118107L) + ((int128)tmp_q[6] * 16969694708707L) + ((int128)tmp_q[7] * 11025890776264L) + ((int128)tmp_q[8] * 43347723845914L) - ((int128)tmp_q[9] * 49324602339088L) + ((int128)tmp_q[10] * 52504023431909L);

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

