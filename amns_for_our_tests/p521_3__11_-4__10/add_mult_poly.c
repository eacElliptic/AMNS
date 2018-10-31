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
	tmp_q[0] = ((uint64_t)op[0] * 362033763411498393UL) + ((((uint64_t)op[1] * 17270937074950858103UL) + ((uint64_t)op[2] * 3051111590016292667UL) + ((uint64_t)op[3] * 940624508921270999UL) + ((uint64_t)op[4] * 10869836600158573566UL) + ((uint64_t)op[5] * 1682341512372201259UL) + ((uint64_t)op[6] * 16315846956822539324UL) + ((uint64_t)op[7] * 7366561357038173502UL) + ((uint64_t)op[8] * 207659372662257641UL) + ((uint64_t)op[9] * 11831601022938119727UL) + ((uint64_t)op[10] * 5838084059091458956UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 5838084059091458956UL) + ((uint64_t)op[1] * 362033763411498393UL) + ((((uint64_t)op[2] * 17270937074950858103UL) + ((uint64_t)op[3] * 3051111590016292667UL) + ((uint64_t)op[4] * 940624508921270999UL) + ((uint64_t)op[5] * 10869836600158573566UL) + ((uint64_t)op[6] * 1682341512372201259UL) + ((uint64_t)op[7] * 16315846956822539324UL) + ((uint64_t)op[8] * 7366561357038173502UL) + ((uint64_t)op[9] * 207659372662257641UL) + ((uint64_t)op[10] * 11831601022938119727UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 11831601022938119727UL) + ((uint64_t)op[1] * 5838084059091458956UL) + ((uint64_t)op[2] * 362033763411498393UL) + ((((uint64_t)op[3] * 17270937074950858103UL) + ((uint64_t)op[4] * 3051111590016292667UL) + ((uint64_t)op[5] * 940624508921270999UL) + ((uint64_t)op[6] * 10869836600158573566UL) + ((uint64_t)op[7] * 1682341512372201259UL) + ((uint64_t)op[8] * 16315846956822539324UL) + ((uint64_t)op[9] * 7366561357038173502UL) + ((uint64_t)op[10] * 207659372662257641UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 207659372662257641UL) + ((uint64_t)op[1] * 11831601022938119727UL) + ((uint64_t)op[2] * 5838084059091458956UL) + ((uint64_t)op[3] * 362033763411498393UL) + ((((uint64_t)op[4] * 17270937074950858103UL) + ((uint64_t)op[5] * 3051111590016292667UL) + ((uint64_t)op[6] * 940624508921270999UL) + ((uint64_t)op[7] * 10869836600158573566UL) + ((uint64_t)op[8] * 1682341512372201259UL) + ((uint64_t)op[9] * 16315846956822539324UL) + ((uint64_t)op[10] * 7366561357038173502UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 7366561357038173502UL) + ((uint64_t)op[1] * 207659372662257641UL) + ((uint64_t)op[2] * 11831601022938119727UL) + ((uint64_t)op[3] * 5838084059091458956UL) + ((uint64_t)op[4] * 362033763411498393UL) + ((((uint64_t)op[5] * 17270937074950858103UL) + ((uint64_t)op[6] * 3051111590016292667UL) + ((uint64_t)op[7] * 940624508921270999UL) + ((uint64_t)op[8] * 10869836600158573566UL) + ((uint64_t)op[9] * 1682341512372201259UL) + ((uint64_t)op[10] * 16315846956822539324UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 16315846956822539324UL) + ((uint64_t)op[1] * 7366561357038173502UL) + ((uint64_t)op[2] * 207659372662257641UL) + ((uint64_t)op[3] * 11831601022938119727UL) + ((uint64_t)op[4] * 5838084059091458956UL) + ((uint64_t)op[5] * 362033763411498393UL) + ((((uint64_t)op[6] * 17270937074950858103UL) + ((uint64_t)op[7] * 3051111590016292667UL) + ((uint64_t)op[8] * 940624508921270999UL) + ((uint64_t)op[9] * 10869836600158573566UL) + ((uint64_t)op[10] * 1682341512372201259UL)) * 18446744073709551612);
	tmp_q[6] = ((uint64_t)op[0] * 1682341512372201259UL) + ((uint64_t)op[1] * 16315846956822539324UL) + ((uint64_t)op[2] * 7366561357038173502UL) + ((uint64_t)op[3] * 207659372662257641UL) + ((uint64_t)op[4] * 11831601022938119727UL) + ((uint64_t)op[5] * 5838084059091458956UL) + ((uint64_t)op[6] * 362033763411498393UL) + ((((uint64_t)op[7] * 17270937074950858103UL) + ((uint64_t)op[8] * 3051111590016292667UL) + ((uint64_t)op[9] * 940624508921270999UL) + ((uint64_t)op[10] * 10869836600158573566UL)) * 18446744073709551612);
	tmp_q[7] = ((uint64_t)op[0] * 10869836600158573566UL) + ((uint64_t)op[1] * 1682341512372201259UL) + ((uint64_t)op[2] * 16315846956822539324UL) + ((uint64_t)op[3] * 7366561357038173502UL) + ((uint64_t)op[4] * 207659372662257641UL) + ((uint64_t)op[5] * 11831601022938119727UL) + ((uint64_t)op[6] * 5838084059091458956UL) + ((uint64_t)op[7] * 362033763411498393UL) + ((((uint64_t)op[8] * 17270937074950858103UL) + ((uint64_t)op[9] * 3051111590016292667UL) + ((uint64_t)op[10] * 940624508921270999UL)) * 18446744073709551612);
	tmp_q[8] = ((uint64_t)op[0] * 940624508921270999UL) + ((uint64_t)op[1] * 10869836600158573566UL) + ((uint64_t)op[2] * 1682341512372201259UL) + ((uint64_t)op[3] * 16315846956822539324UL) + ((uint64_t)op[4] * 7366561357038173502UL) + ((uint64_t)op[5] * 207659372662257641UL) + ((uint64_t)op[6] * 11831601022938119727UL) + ((uint64_t)op[7] * 5838084059091458956UL) + ((uint64_t)op[8] * 362033763411498393UL) + ((((uint64_t)op[9] * 17270937074950858103UL) + ((uint64_t)op[10] * 3051111590016292667UL)) * 18446744073709551612);
	tmp_q[9] = ((uint64_t)op[0] * 3051111590016292667UL) + ((uint64_t)op[1] * 940624508921270999UL) + ((uint64_t)op[2] * 10869836600158573566UL) + ((uint64_t)op[3] * 1682341512372201259UL) + ((uint64_t)op[4] * 16315846956822539324UL) + ((uint64_t)op[5] * 7366561357038173502UL) + ((uint64_t)op[6] * 207659372662257641UL) + ((uint64_t)op[7] * 11831601022938119727UL) + ((uint64_t)op[8] * 5838084059091458956UL) + ((uint64_t)op[9] * 362033763411498393UL) + ((uint64_t)op[10] * 4703227995034774052UL);
	tmp_q[10] = ((uint64_t)op[0] * 17270937074950858103UL) + ((uint64_t)op[1] * 3051111590016292667UL) + ((uint64_t)op[2] * 940624508921270999UL) + ((uint64_t)op[3] * 10869836600158573566UL) + ((uint64_t)op[4] * 1682341512372201259UL) + ((uint64_t)op[5] * 16315846956822539324UL) + ((uint64_t)op[6] * 7366561357038173502UL) + ((uint64_t)op[7] * 207659372662257641UL) + ((uint64_t)op[8] * 11831601022938119727UL) + ((uint64_t)op[9] * 5838084059091458956UL) + ((uint64_t)op[10] * 362033763411498393UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 40943142513533L) - ((-((int128)tmp_q[1] * 7148361980215L) + ((int128)tmp_q[2] * 90529968971046L) - ((int128)tmp_q[3] * 31846499259517L) - ((int128)tmp_q[4] * 107941585658099L) - ((int128)tmp_q[5] * 36375890759163L) - ((int128)tmp_q[6] * 32954281998326L) - ((int128)tmp_q[7] * 47857097858667L) - ((int128)tmp_q[8] * 75318355649475L) - ((int128)tmp_q[9] * 18800035801149L) + ((int128)tmp_q[10] * 8286696503120L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 8286696503120L) - ((int128)tmp_q[1] * 40943142513533L) - ((-((int128)tmp_q[2] * 7148361980215L) + ((int128)tmp_q[3] * 90529968971046L) - ((int128)tmp_q[4] * 31846499259517L) - ((int128)tmp_q[5] * 107941585658099L) - ((int128)tmp_q[6] * 36375890759163L) - ((int128)tmp_q[7] * 32954281998326L) - ((int128)tmp_q[8] * 47857097858667L) - ((int128)tmp_q[9] * 75318355649475L) - ((int128)tmp_q[10] * 18800035801149L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 18800035801149L) + ((int128)tmp_q[1] * 8286696503120L) - ((int128)tmp_q[2] * 40943142513533L) - ((-((int128)tmp_q[3] * 7148361980215L) + ((int128)tmp_q[4] * 90529968971046L) - ((int128)tmp_q[5] * 31846499259517L) - ((int128)tmp_q[6] * 107941585658099L) - ((int128)tmp_q[7] * 36375890759163L) - ((int128)tmp_q[8] * 32954281998326L) - ((int128)tmp_q[9] * 47857097858667L) - ((int128)tmp_q[10] * 75318355649475L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 75318355649475L) - ((int128)tmp_q[1] * 18800035801149L) + ((int128)tmp_q[2] * 8286696503120L) - ((int128)tmp_q[3] * 40943142513533L) - ((-((int128)tmp_q[4] * 7148361980215L) + ((int128)tmp_q[5] * 90529968971046L) - ((int128)tmp_q[6] * 31846499259517L) - ((int128)tmp_q[7] * 107941585658099L) - ((int128)tmp_q[8] * 36375890759163L) - ((int128)tmp_q[9] * 32954281998326L) - ((int128)tmp_q[10] * 47857097858667L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 47857097858667L) - ((int128)tmp_q[1] * 75318355649475L) - ((int128)tmp_q[2] * 18800035801149L) + ((int128)tmp_q[3] * 8286696503120L) - ((int128)tmp_q[4] * 40943142513533L) - ((-((int128)tmp_q[5] * 7148361980215L) + ((int128)tmp_q[6] * 90529968971046L) - ((int128)tmp_q[7] * 31846499259517L) - ((int128)tmp_q[8] * 107941585658099L) - ((int128)tmp_q[9] * 36375890759163L) - ((int128)tmp_q[10] * 32954281998326L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 32954281998326L) - ((int128)tmp_q[1] * 47857097858667L) - ((int128)tmp_q[2] * 75318355649475L) - ((int128)tmp_q[3] * 18800035801149L) + ((int128)tmp_q[4] * 8286696503120L) - ((int128)tmp_q[5] * 40943142513533L) - ((-((int128)tmp_q[6] * 7148361980215L) + ((int128)tmp_q[7] * 90529968971046L) - ((int128)tmp_q[8] * 31846499259517L) - ((int128)tmp_q[9] * 107941585658099L) - ((int128)tmp_q[10] * 36375890759163L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 36375890759163L) - ((int128)tmp_q[1] * 32954281998326L) - ((int128)tmp_q[2] * 47857097858667L) - ((int128)tmp_q[3] * 75318355649475L) - ((int128)tmp_q[4] * 18800035801149L) + ((int128)tmp_q[5] * 8286696503120L) - ((int128)tmp_q[6] * 40943142513533L) - ((-((int128)tmp_q[7] * 7148361980215L) + ((int128)tmp_q[8] * 90529968971046L) - ((int128)tmp_q[9] * 31846499259517L) - ((int128)tmp_q[10] * 107941585658099L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 107941585658099L) - ((int128)tmp_q[1] * 36375890759163L) - ((int128)tmp_q[2] * 32954281998326L) - ((int128)tmp_q[3] * 47857097858667L) - ((int128)tmp_q[4] * 75318355649475L) - ((int128)tmp_q[5] * 18800035801149L) + ((int128)tmp_q[6] * 8286696503120L) - ((int128)tmp_q[7] * 40943142513533L) - ((-((int128)tmp_q[8] * 7148361980215L) + ((int128)tmp_q[9] * 90529968971046L) - ((int128)tmp_q[10] * 31846499259517L)) * 4);
	tmp_zero[8] = -((int128)tmp_q[0] * 31846499259517L) - ((int128)tmp_q[1] * 107941585658099L) - ((int128)tmp_q[2] * 36375890759163L) - ((int128)tmp_q[3] * 32954281998326L) - ((int128)tmp_q[4] * 47857097858667L) - ((int128)tmp_q[5] * 75318355649475L) - ((int128)tmp_q[6] * 18800035801149L) + ((int128)tmp_q[7] * 8286696503120L) - ((int128)tmp_q[8] * 40943142513533L) - ((-((int128)tmp_q[9] * 7148361980215L) + ((int128)tmp_q[10] * 90529968971046L)) * 4);
	tmp_zero[9] = ((int128)tmp_q[0] * 90529968971046L) - ((int128)tmp_q[1] * 31846499259517L) - ((int128)tmp_q[2] * 107941585658099L) - ((int128)tmp_q[3] * 36375890759163L) - ((int128)tmp_q[4] * 32954281998326L) - ((int128)tmp_q[5] * 47857097858667L) - ((int128)tmp_q[6] * 75318355649475L) - ((int128)tmp_q[7] * 18800035801149L) + ((int128)tmp_q[8] * 8286696503120L) - ((int128)tmp_q[9] * 40943142513533L) + ((int128)tmp_q[10] * 28593447920860L);
	tmp_zero[10] = -((int128)tmp_q[0] * 7148361980215L) + ((int128)tmp_q[1] * 90529968971046L) - ((int128)tmp_q[2] * 31846499259517L) - ((int128)tmp_q[3] * 107941585658099L) - ((int128)tmp_q[4] * 36375890759163L) - ((int128)tmp_q[5] * 32954281998326L) - ((int128)tmp_q[6] * 47857097858667L) - ((int128)tmp_q[7] * 75318355649475L) - ((int128)tmp_q[8] * 18800035801149L) + ((int128)tmp_q[9] * 8286696503120L) - ((int128)tmp_q[10] * 40943142513533L);

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

