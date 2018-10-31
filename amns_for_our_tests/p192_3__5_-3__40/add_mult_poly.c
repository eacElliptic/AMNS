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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10652602987574287389UL) + ((((uint64_t)op[1] * 16113709045207005224UL) + ((uint64_t)op[2] * 11061857503245314597UL) + ((uint64_t)op[3] * 9842603582955748185UL) + ((uint64_t)op[4] * 11837415693332255598UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 11837415693332255598UL) + ((uint64_t)op[1] * 10652602987574287389UL) + ((((uint64_t)op[2] * 16113709045207005224UL) + ((uint64_t)op[3] * 11061857503245314597UL) + ((uint64_t)op[4] * 9842603582955748185UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 9842603582955748185UL) + ((uint64_t)op[1] * 11837415693332255598UL) + ((uint64_t)op[2] * 10652602987574287389UL) + ((((uint64_t)op[3] * 16113709045207005224UL) + ((uint64_t)op[4] * 11061857503245314597UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 11061857503245314597UL) + ((uint64_t)op[1] * 9842603582955748185UL) + ((uint64_t)op[2] * 11837415693332255598UL) + ((uint64_t)op[3] * 10652602987574287389UL) + ((uint64_t)op[4] * 6999105085507639176UL);
	tmp_q[4] = ((uint64_t)op[0] * 16113709045207005224UL) + ((uint64_t)op[1] * 11061857503245314597UL) + ((uint64_t)op[2] * 9842603582955748185UL) + ((uint64_t)op[3] * 11837415693332255598UL) + ((uint64_t)op[4] * 10652602987574287389UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 373727665381L) - ((-((int128)tmp_q[1] * 34838666351L) + ((int128)tmp_q[2] * 105894189676L) - ((int128)tmp_q[3] * 44013271466L) + ((int128)tmp_q[4] * 169546971593L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 169546971593L) - ((int128)tmp_q[1] * 373727665381L) - ((-((int128)tmp_q[2] * 34838666351L) + ((int128)tmp_q[3] * 105894189676L) - ((int128)tmp_q[4] * 44013271466L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 44013271466L) + ((int128)tmp_q[1] * 169546971593L) - ((int128)tmp_q[2] * 373727665381L) - ((-((int128)tmp_q[3] * 34838666351L) + ((int128)tmp_q[4] * 105894189676L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 105894189676L) - ((int128)tmp_q[1] * 44013271466L) + ((int128)tmp_q[2] * 169546971593L) - ((int128)tmp_q[3] * 373727665381L) + ((int128)tmp_q[4] * 104515999053L);
	tmp_zero[4] = -((int128)tmp_q[0] * 34838666351L) + ((int128)tmp_q[1] * 105894189676L) - ((int128)tmp_q[2] * 44013271466L) + ((int128)tmp_q[3] * 169546971593L) - ((int128)tmp_q[4] * 373727665381L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

