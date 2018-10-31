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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10014518919023341277UL) + ((((uint64_t)op[1] * 17560870390019442219UL) + ((uint64_t)op[2] * 11448234668118483157UL) + ((uint64_t)op[3] * 16824160668370844516UL) + ((uint64_t)op[4] * 830384849022258780UL) + ((uint64_t)op[5] * 17357284534051305828UL) + ((uint64_t)op[6] * 1999320084379395252UL) + ((uint64_t)op[7] * 648762251525377864UL) + ((uint64_t)op[8] * 17388451311684903007UL) + ((uint64_t)op[9] * 14757455458042504508UL) + ((uint64_t)op[10] * 18037587355836302125UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 18037587355836302125UL) + ((uint64_t)op[1] * 10014518919023341277UL) + ((((uint64_t)op[2] * 17560870390019442219UL) + ((uint64_t)op[3] * 11448234668118483157UL) + ((uint64_t)op[4] * 16824160668370844516UL) + ((uint64_t)op[5] * 830384849022258780UL) + ((uint64_t)op[6] * 17357284534051305828UL) + ((uint64_t)op[7] * 1999320084379395252UL) + ((uint64_t)op[8] * 648762251525377864UL) + ((uint64_t)op[9] * 17388451311684903007UL) + ((uint64_t)op[10] * 14757455458042504508UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 14757455458042504508UL) + ((uint64_t)op[1] * 18037587355836302125UL) + ((uint64_t)op[2] * 10014518919023341277UL) + ((((uint64_t)op[3] * 17560870390019442219UL) + ((uint64_t)op[4] * 11448234668118483157UL) + ((uint64_t)op[5] * 16824160668370844516UL) + ((uint64_t)op[6] * 830384849022258780UL) + ((uint64_t)op[7] * 17357284534051305828UL) + ((uint64_t)op[8] * 1999320084379395252UL) + ((uint64_t)op[9] * 648762251525377864UL) + ((uint64_t)op[10] * 17388451311684903007UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 17388451311684903007UL) + ((uint64_t)op[1] * 14757455458042504508UL) + ((uint64_t)op[2] * 18037587355836302125UL) + ((uint64_t)op[3] * 10014518919023341277UL) + ((((uint64_t)op[4] * 17560870390019442219UL) + ((uint64_t)op[5] * 11448234668118483157UL) + ((uint64_t)op[6] * 16824160668370844516UL) + ((uint64_t)op[7] * 830384849022258780UL) + ((uint64_t)op[8] * 17357284534051305828UL) + ((uint64_t)op[9] * 1999320084379395252UL) + ((uint64_t)op[10] * 648762251525377864UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 648762251525377864UL) + ((uint64_t)op[1] * 17388451311684903007UL) + ((uint64_t)op[2] * 14757455458042504508UL) + ((uint64_t)op[3] * 18037587355836302125UL) + ((uint64_t)op[4] * 10014518919023341277UL) + ((((uint64_t)op[5] * 17560870390019442219UL) + ((uint64_t)op[6] * 11448234668118483157UL) + ((uint64_t)op[7] * 16824160668370844516UL) + ((uint64_t)op[8] * 830384849022258780UL) + ((uint64_t)op[9] * 17357284534051305828UL) + ((uint64_t)op[10] * 1999320084379395252UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 1999320084379395252UL) + ((uint64_t)op[1] * 648762251525377864UL) + ((uint64_t)op[2] * 17388451311684903007UL) + ((uint64_t)op[3] * 14757455458042504508UL) + ((uint64_t)op[4] * 18037587355836302125UL) + ((uint64_t)op[5] * 10014518919023341277UL) + ((((uint64_t)op[6] * 17560870390019442219UL) + ((uint64_t)op[7] * 11448234668118483157UL) + ((uint64_t)op[8] * 16824160668370844516UL) + ((uint64_t)op[9] * 830384849022258780UL) + ((uint64_t)op[10] * 17357284534051305828UL)) * 18446744073709551612);
	tmp_q[6] = ((uint64_t)op[0] * 17357284534051305828UL) + ((uint64_t)op[1] * 1999320084379395252UL) + ((uint64_t)op[2] * 648762251525377864UL) + ((uint64_t)op[3] * 17388451311684903007UL) + ((uint64_t)op[4] * 14757455458042504508UL) + ((uint64_t)op[5] * 18037587355836302125UL) + ((uint64_t)op[6] * 10014518919023341277UL) + ((((uint64_t)op[7] * 17560870390019442219UL) + ((uint64_t)op[8] * 11448234668118483157UL) + ((uint64_t)op[9] * 16824160668370844516UL) + ((uint64_t)op[10] * 830384849022258780UL)) * 18446744073709551612);
	tmp_q[7] = ((uint64_t)op[0] * 830384849022258780UL) + ((uint64_t)op[1] * 17357284534051305828UL) + ((uint64_t)op[2] * 1999320084379395252UL) + ((uint64_t)op[3] * 648762251525377864UL) + ((uint64_t)op[4] * 17388451311684903007UL) + ((uint64_t)op[5] * 14757455458042504508UL) + ((uint64_t)op[6] * 18037587355836302125UL) + ((uint64_t)op[7] * 10014518919023341277UL) + ((((uint64_t)op[8] * 17560870390019442219UL) + ((uint64_t)op[9] * 11448234668118483157UL) + ((uint64_t)op[10] * 16824160668370844516UL)) * 18446744073709551612);
	tmp_q[8] = ((uint64_t)op[0] * 16824160668370844516UL) + ((uint64_t)op[1] * 830384849022258780UL) + ((uint64_t)op[2] * 17357284534051305828UL) + ((uint64_t)op[3] * 1999320084379395252UL) + ((uint64_t)op[4] * 648762251525377864UL) + ((uint64_t)op[5] * 17388451311684903007UL) + ((uint64_t)op[6] * 14757455458042504508UL) + ((uint64_t)op[7] * 18037587355836302125UL) + ((uint64_t)op[8] * 10014518919023341277UL) + ((((uint64_t)op[9] * 17560870390019442219UL) + ((uint64_t)op[10] * 11448234668118483157UL)) * 18446744073709551612);
	tmp_q[9] = ((uint64_t)op[0] * 11448234668118483157UL) + ((uint64_t)op[1] * 16824160668370844516UL) + ((uint64_t)op[2] * 830384849022258780UL) + ((uint64_t)op[3] * 17357284534051305828UL) + ((uint64_t)op[4] * 1999320084379395252UL) + ((uint64_t)op[5] * 648762251525377864UL) + ((uint64_t)op[6] * 17388451311684903007UL) + ((uint64_t)op[7] * 14757455458042504508UL) + ((uint64_t)op[8] * 18037587355836302125UL) + ((uint64_t)op[9] * 10014518919023341277UL) + ((uint64_t)op[10] * 3543494734760437588UL);
	tmp_q[10] = ((uint64_t)op[0] * 17560870390019442219UL) + ((uint64_t)op[1] * 11448234668118483157UL) + ((uint64_t)op[2] * 16824160668370844516UL) + ((uint64_t)op[3] * 830384849022258780UL) + ((uint64_t)op[4] * 17357284534051305828UL) + ((uint64_t)op[5] * 1999320084379395252UL) + ((uint64_t)op[6] * 648762251525377864UL) + ((uint64_t)op[7] * 17388451311684903007UL) + ((uint64_t)op[8] * 14757455458042504508UL) + ((uint64_t)op[9] * 18037587355836302125UL) + ((uint64_t)op[10] * 10014518919023341277UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 58172758618453L) - ((-((int128)tmp_q[1] * 50519821707227L) + ((int128)tmp_q[2] * 77921646312768L) + ((int128)tmp_q[3] * 34977450810031L) + ((int128)tmp_q[4] * 17533215089699L) + ((int128)tmp_q[5] * 50888306895870L) + ((int128)tmp_q[6] * 69508849792362L) + ((int128)tmp_q[7] * 10895827869093L) - ((int128)tmp_q[8] * 60078887577812L) + ((int128)tmp_q[9] * 22978784946131L) - ((int128)tmp_q[10] * 62308065051591L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 62308065051591L) - ((int128)tmp_q[1] * 58172758618453L) - ((-((int128)tmp_q[2] * 50519821707227L) + ((int128)tmp_q[3] * 77921646312768L) + ((int128)tmp_q[4] * 34977450810031L) + ((int128)tmp_q[5] * 17533215089699L) + ((int128)tmp_q[6] * 50888306895870L) + ((int128)tmp_q[7] * 69508849792362L) + ((int128)tmp_q[8] * 10895827869093L) - ((int128)tmp_q[9] * 60078887577812L) + ((int128)tmp_q[10] * 22978784946131L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 22978784946131L) - ((int128)tmp_q[1] * 62308065051591L) - ((int128)tmp_q[2] * 58172758618453L) - ((-((int128)tmp_q[3] * 50519821707227L) + ((int128)tmp_q[4] * 77921646312768L) + ((int128)tmp_q[5] * 34977450810031L) + ((int128)tmp_q[6] * 17533215089699L) + ((int128)tmp_q[7] * 50888306895870L) + ((int128)tmp_q[8] * 69508849792362L) + ((int128)tmp_q[9] * 10895827869093L) - ((int128)tmp_q[10] * 60078887577812L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 60078887577812L) + ((int128)tmp_q[1] * 22978784946131L) - ((int128)tmp_q[2] * 62308065051591L) - ((int128)tmp_q[3] * 58172758618453L) - ((-((int128)tmp_q[4] * 50519821707227L) + ((int128)tmp_q[5] * 77921646312768L) + ((int128)tmp_q[6] * 34977450810031L) + ((int128)tmp_q[7] * 17533215089699L) + ((int128)tmp_q[8] * 50888306895870L) + ((int128)tmp_q[9] * 69508849792362L) + ((int128)tmp_q[10] * 10895827869093L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 10895827869093L) - ((int128)tmp_q[1] * 60078887577812L) + ((int128)tmp_q[2] * 22978784946131L) - ((int128)tmp_q[3] * 62308065051591L) - ((int128)tmp_q[4] * 58172758618453L) - ((-((int128)tmp_q[5] * 50519821707227L) + ((int128)tmp_q[6] * 77921646312768L) + ((int128)tmp_q[7] * 34977450810031L) + ((int128)tmp_q[8] * 17533215089699L) + ((int128)tmp_q[9] * 50888306895870L) + ((int128)tmp_q[10] * 69508849792362L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 69508849792362L) + ((int128)tmp_q[1] * 10895827869093L) - ((int128)tmp_q[2] * 60078887577812L) + ((int128)tmp_q[3] * 22978784946131L) - ((int128)tmp_q[4] * 62308065051591L) - ((int128)tmp_q[5] * 58172758618453L) - ((-((int128)tmp_q[6] * 50519821707227L) + ((int128)tmp_q[7] * 77921646312768L) + ((int128)tmp_q[8] * 34977450810031L) + ((int128)tmp_q[9] * 17533215089699L) + ((int128)tmp_q[10] * 50888306895870L)) * 4);
	tmp_zero[6] = ((int128)tmp_q[0] * 50888306895870L) + ((int128)tmp_q[1] * 69508849792362L) + ((int128)tmp_q[2] * 10895827869093L) - ((int128)tmp_q[3] * 60078887577812L) + ((int128)tmp_q[4] * 22978784946131L) - ((int128)tmp_q[5] * 62308065051591L) - ((int128)tmp_q[6] * 58172758618453L) - ((-((int128)tmp_q[7] * 50519821707227L) + ((int128)tmp_q[8] * 77921646312768L) + ((int128)tmp_q[9] * 34977450810031L) + ((int128)tmp_q[10] * 17533215089699L)) * 4);
	tmp_zero[7] = ((int128)tmp_q[0] * 17533215089699L) + ((int128)tmp_q[1] * 50888306895870L) + ((int128)tmp_q[2] * 69508849792362L) + ((int128)tmp_q[3] * 10895827869093L) - ((int128)tmp_q[4] * 60078887577812L) + ((int128)tmp_q[5] * 22978784946131L) - ((int128)tmp_q[6] * 62308065051591L) - ((int128)tmp_q[7] * 58172758618453L) - ((-((int128)tmp_q[8] * 50519821707227L) + ((int128)tmp_q[9] * 77921646312768L) + ((int128)tmp_q[10] * 34977450810031L)) * 4);
	tmp_zero[8] = ((int128)tmp_q[0] * 34977450810031L) + ((int128)tmp_q[1] * 17533215089699L) + ((int128)tmp_q[2] * 50888306895870L) + ((int128)tmp_q[3] * 69508849792362L) + ((int128)tmp_q[4] * 10895827869093L) - ((int128)tmp_q[5] * 60078887577812L) + ((int128)tmp_q[6] * 22978784946131L) - ((int128)tmp_q[7] * 62308065051591L) - ((int128)tmp_q[8] * 58172758618453L) - ((-((int128)tmp_q[9] * 50519821707227L) + ((int128)tmp_q[10] * 77921646312768L)) * 4);
	tmp_zero[9] = ((int128)tmp_q[0] * 77921646312768L) + ((int128)tmp_q[1] * 34977450810031L) + ((int128)tmp_q[2] * 17533215089699L) + ((int128)tmp_q[3] * 50888306895870L) + ((int128)tmp_q[4] * 69508849792362L) + ((int128)tmp_q[5] * 10895827869093L) - ((int128)tmp_q[6] * 60078887577812L) + ((int128)tmp_q[7] * 22978784946131L) - ((int128)tmp_q[8] * 62308065051591L) - ((int128)tmp_q[9] * 58172758618453L) + ((int128)tmp_q[10] * 202079286828908L);
	tmp_zero[10] = -((int128)tmp_q[0] * 50519821707227L) + ((int128)tmp_q[1] * 77921646312768L) + ((int128)tmp_q[2] * 34977450810031L) + ((int128)tmp_q[3] * 17533215089699L) + ((int128)tmp_q[4] * 50888306895870L) + ((int128)tmp_q[5] * 69508849792362L) + ((int128)tmp_q[6] * 10895827869093L) - ((int128)tmp_q[7] * 60078887577812L) + ((int128)tmp_q[8] * 22978784946131L) - ((int128)tmp_q[9] * 62308065051591L) - ((int128)tmp_q[10] * 58172758618453L);

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

