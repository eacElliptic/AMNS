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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17342627491456163767UL) + ((((uint64_t)op[1] * 4650668328407513906UL) + ((uint64_t)op[2] * 3135107812845445692UL) + ((uint64_t)op[3] * 2734926780066899882UL) + ((uint64_t)op[4] * 1801354727202932945UL) + ((uint64_t)op[5] * 9837381880291603952UL) + ((uint64_t)op[6] * 15504498958281627786UL) + ((uint64_t)op[7] * 12458371419583601566UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 12458371419583601566UL) + ((uint64_t)op[1] * 17342627491456163767UL) + ((((uint64_t)op[2] * 4650668328407513906UL) + ((uint64_t)op[3] * 3135107812845445692UL) + ((uint64_t)op[4] * 2734926780066899882UL) + ((uint64_t)op[5] * 1801354727202932945UL) + ((uint64_t)op[6] * 9837381880291603952UL) + ((uint64_t)op[7] * 15504498958281627786UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 15504498958281627786UL) + ((uint64_t)op[1] * 12458371419583601566UL) + ((uint64_t)op[2] * 17342627491456163767UL) + ((((uint64_t)op[3] * 4650668328407513906UL) + ((uint64_t)op[4] * 3135107812845445692UL) + ((uint64_t)op[5] * 2734926780066899882UL) + ((uint64_t)op[6] * 1801354727202932945UL) + ((uint64_t)op[7] * 9837381880291603952UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 9837381880291603952UL) + ((uint64_t)op[1] * 15504498958281627786UL) + ((uint64_t)op[2] * 12458371419583601566UL) + ((uint64_t)op[3] * 17342627491456163767UL) + ((((uint64_t)op[4] * 4650668328407513906UL) + ((uint64_t)op[5] * 3135107812845445692UL) + ((uint64_t)op[6] * 2734926780066899882UL) + ((uint64_t)op[7] * 1801354727202932945UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 1801354727202932945UL) + ((uint64_t)op[1] * 9837381880291603952UL) + ((uint64_t)op[2] * 15504498958281627786UL) + ((uint64_t)op[3] * 12458371419583601566UL) + ((uint64_t)op[4] * 17342627491456163767UL) + ((((uint64_t)op[5] * 4650668328407513906UL) + ((uint64_t)op[6] * 3135107812845445692UL) + ((uint64_t)op[7] * 2734926780066899882UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 2734926780066899882UL) + ((uint64_t)op[1] * 1801354727202932945UL) + ((uint64_t)op[2] * 9837381880291603952UL) + ((uint64_t)op[3] * 15504498958281627786UL) + ((uint64_t)op[4] * 12458371419583601566UL) + ((uint64_t)op[5] * 17342627491456163767UL) + ((((uint64_t)op[6] * 4650668328407513906UL) + ((uint64_t)op[7] * 3135107812845445692UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 3135107812845445692UL) + ((uint64_t)op[1] * 2734926780066899882UL) + ((uint64_t)op[2] * 1801354727202932945UL) + ((uint64_t)op[3] * 9837381880291603952UL) + ((uint64_t)op[4] * 15504498958281627786UL) + ((uint64_t)op[5] * 12458371419583601566UL) + ((uint64_t)op[6] * 17342627491456163767UL) + ((uint64_t)op[7] * 9145407416894523804UL);
	tmp_q[7] = ((uint64_t)op[0] * 4650668328407513906UL) + ((uint64_t)op[1] * 3135107812845445692UL) + ((uint64_t)op[2] * 2734926780066899882UL) + ((uint64_t)op[3] * 1801354727202932945UL) + ((uint64_t)op[4] * 9837381880291603952UL) + ((uint64_t)op[5] * 15504498958281627786UL) + ((uint64_t)op[6] * 12458371419583601566UL) + ((uint64_t)op[7] * 17342627491456163767UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 8673078415733L) - ((-((int128)tmp_q[1] * 24435120248554L) - ((int128)tmp_q[2] * 29600845977212L) - ((int128)tmp_q[3] * 104110181510726L) + ((int128)tmp_q[4] * 56650789276039L) + ((int128)tmp_q[5] * 126144696304648L) - ((int128)tmp_q[6] * 99850844940046L) + ((int128)tmp_q[7] * 171348083429370L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 171348083429370L) - ((int128)tmp_q[1] * 8673078415733L) - ((-((int128)tmp_q[2] * 24435120248554L) - ((int128)tmp_q[3] * 29600845977212L) - ((int128)tmp_q[4] * 104110181510726L) + ((int128)tmp_q[5] * 56650789276039L) + ((int128)tmp_q[6] * 126144696304648L) - ((int128)tmp_q[7] * 99850844940046L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 99850844940046L) + ((int128)tmp_q[1] * 171348083429370L) - ((int128)tmp_q[2] * 8673078415733L) - ((-((int128)tmp_q[3] * 24435120248554L) - ((int128)tmp_q[4] * 29600845977212L) - ((int128)tmp_q[5] * 104110181510726L) + ((int128)tmp_q[6] * 56650789276039L) + ((int128)tmp_q[7] * 126144696304648L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 126144696304648L) - ((int128)tmp_q[1] * 99850844940046L) + ((int128)tmp_q[2] * 171348083429370L) - ((int128)tmp_q[3] * 8673078415733L) - ((-((int128)tmp_q[4] * 24435120248554L) - ((int128)tmp_q[5] * 29600845977212L) - ((int128)tmp_q[6] * 104110181510726L) + ((int128)tmp_q[7] * 56650789276039L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 56650789276039L) + ((int128)tmp_q[1] * 126144696304648L) - ((int128)tmp_q[2] * 99850844940046L) + ((int128)tmp_q[3] * 171348083429370L) - ((int128)tmp_q[4] * 8673078415733L) - ((-((int128)tmp_q[5] * 24435120248554L) - ((int128)tmp_q[6] * 29600845977212L) - ((int128)tmp_q[7] * 104110181510726L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 104110181510726L) + ((int128)tmp_q[1] * 56650789276039L) + ((int128)tmp_q[2] * 126144696304648L) - ((int128)tmp_q[3] * 99850844940046L) + ((int128)tmp_q[4] * 171348083429370L) - ((int128)tmp_q[5] * 8673078415733L) - ((-((int128)tmp_q[6] * 24435120248554L) - ((int128)tmp_q[7] * 29600845977212L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 29600845977212L) - ((int128)tmp_q[1] * 104110181510726L) + ((int128)tmp_q[2] * 56650789276039L) + ((int128)tmp_q[3] * 126144696304648L) - ((int128)tmp_q[4] * 99850844940046L) + ((int128)tmp_q[5] * 171348083429370L) - ((int128)tmp_q[6] * 8673078415733L) + ((int128)tmp_q[7] * 48870240497108L);
	tmp_zero[7] = -((int128)tmp_q[0] * 24435120248554L) - ((int128)tmp_q[1] * 29600845977212L) - ((int128)tmp_q[2] * 104110181510726L) + ((int128)tmp_q[3] * 56650789276039L) + ((int128)tmp_q[4] * 126144696304648L) - ((int128)tmp_q[5] * 99850844940046L) + ((int128)tmp_q[6] * 171348083429370L) - ((int128)tmp_q[7] * 8673078415733L);

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

