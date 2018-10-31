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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16073276782297802907UL) + ((((uint64_t)op[1] * 6825513530090253624UL) + ((uint64_t)op[2] * 15922561173780807793UL) + ((uint64_t)op[3] * 3256724135995114520UL) + ((uint64_t)op[4] * 13998813185930542666UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 13998813185930542666UL) + ((uint64_t)op[1] * 16073276782297802907UL) + ((((uint64_t)op[2] * 6825513530090253624UL) + ((uint64_t)op[3] * 15922561173780807793UL) + ((uint64_t)op[4] * 3256724135995114520UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 3256724135995114520UL) + ((uint64_t)op[1] * 13998813185930542666UL) + ((uint64_t)op[2] * 16073276782297802907UL) + ((((uint64_t)op[3] * 6825513530090253624UL) + ((uint64_t)op[4] * 15922561173780807793UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 15922561173780807793UL) + ((uint64_t)op[1] * 3256724135995114520UL) + ((uint64_t)op[2] * 13998813185930542666UL) + ((uint64_t)op[3] * 16073276782297802907UL) + ((uint64_t)op[4] * 736123980406625856UL);
	tmp_q[4] = ((uint64_t)op[0] * 6825513530090253624UL) + ((uint64_t)op[1] * 15922561173780807793UL) + ((uint64_t)op[2] * 3256724135995114520UL) + ((uint64_t)op[3] * 13998813185930542666UL) + ((uint64_t)op[4] * 16073276782297802907UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 260012889325683L) - ((-((int128)tmp_q[1] * 482048312838652L) + ((int128)tmp_q[2] * 199248123841697L) - ((int128)tmp_q[3] * 51305683964196L) + ((int128)tmp_q[4] * 197340263445682L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 197340263445682L) - ((int128)tmp_q[1] * 260012889325683L) - ((-((int128)tmp_q[2] * 482048312838652L) + ((int128)tmp_q[3] * 199248123841697L) - ((int128)tmp_q[4] * 51305683964196L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 51305683964196L) + ((int128)tmp_q[1] * 197340263445682L) - ((int128)tmp_q[2] * 260012889325683L) - ((-((int128)tmp_q[3] * 482048312838652L) + ((int128)tmp_q[4] * 199248123841697L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 199248123841697L) - ((int128)tmp_q[1] * 51305683964196L) + ((int128)tmp_q[2] * 197340263445682L) - ((int128)tmp_q[3] * 260012889325683L) + ((int128)tmp_q[4] * 3856386502709216L);
	tmp_zero[4] = -((int128)tmp_q[0] * 482048312838652L) + ((int128)tmp_q[1] * 199248123841697L) - ((int128)tmp_q[2] * 51305683964196L) + ((int128)tmp_q[3] * 197340263445682L) - ((int128)tmp_q[4] * 260012889325683L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

