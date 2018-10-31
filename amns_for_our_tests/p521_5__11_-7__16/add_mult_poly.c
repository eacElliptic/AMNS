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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9398408090827270617UL) + ((((uint64_t)op[1] * 7314918444598304289UL) + ((uint64_t)op[2] * 15010174984106496879UL) + ((uint64_t)op[3] * 8695639198555464415UL) + ((uint64_t)op[4] * 415742185652332287UL) + ((uint64_t)op[5] * 14971330595751891135UL) + ((uint64_t)op[6] * 1199911456779291470UL) + ((uint64_t)op[7] * 16079450188461417004UL) + ((uint64_t)op[8] * 13104052030274545521UL) + ((uint64_t)op[9] * 1873511197073093991UL) + ((uint64_t)op[10] * 6771460835154525505UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 6771460835154525505UL) + ((uint64_t)op[1] * 9398408090827270617UL) + ((((uint64_t)op[2] * 7314918444598304289UL) + ((uint64_t)op[3] * 15010174984106496879UL) + ((uint64_t)op[4] * 8695639198555464415UL) + ((uint64_t)op[5] * 415742185652332287UL) + ((uint64_t)op[6] * 14971330595751891135UL) + ((uint64_t)op[7] * 1199911456779291470UL) + ((uint64_t)op[8] * 16079450188461417004UL) + ((uint64_t)op[9] * 13104052030274545521UL) + ((uint64_t)op[10] * 1873511197073093991UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 1873511197073093991UL) + ((uint64_t)op[1] * 6771460835154525505UL) + ((uint64_t)op[2] * 9398408090827270617UL) + ((((uint64_t)op[3] * 7314918444598304289UL) + ((uint64_t)op[4] * 15010174984106496879UL) + ((uint64_t)op[5] * 8695639198555464415UL) + ((uint64_t)op[6] * 415742185652332287UL) + ((uint64_t)op[7] * 14971330595751891135UL) + ((uint64_t)op[8] * 1199911456779291470UL) + ((uint64_t)op[9] * 16079450188461417004UL) + ((uint64_t)op[10] * 13104052030274545521UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 13104052030274545521UL) + ((uint64_t)op[1] * 1873511197073093991UL) + ((uint64_t)op[2] * 6771460835154525505UL) + ((uint64_t)op[3] * 9398408090827270617UL) + ((((uint64_t)op[4] * 7314918444598304289UL) + ((uint64_t)op[5] * 15010174984106496879UL) + ((uint64_t)op[6] * 8695639198555464415UL) + ((uint64_t)op[7] * 415742185652332287UL) + ((uint64_t)op[8] * 14971330595751891135UL) + ((uint64_t)op[9] * 1199911456779291470UL) + ((uint64_t)op[10] * 16079450188461417004UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 16079450188461417004UL) + ((uint64_t)op[1] * 13104052030274545521UL) + ((uint64_t)op[2] * 1873511197073093991UL) + ((uint64_t)op[3] * 6771460835154525505UL) + ((uint64_t)op[4] * 9398408090827270617UL) + ((((uint64_t)op[5] * 7314918444598304289UL) + ((uint64_t)op[6] * 15010174984106496879UL) + ((uint64_t)op[7] * 8695639198555464415UL) + ((uint64_t)op[8] * 415742185652332287UL) + ((uint64_t)op[9] * 14971330595751891135UL) + ((uint64_t)op[10] * 1199911456779291470UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 1199911456779291470UL) + ((uint64_t)op[1] * 16079450188461417004UL) + ((uint64_t)op[2] * 13104052030274545521UL) + ((uint64_t)op[3] * 1873511197073093991UL) + ((uint64_t)op[4] * 6771460835154525505UL) + ((uint64_t)op[5] * 9398408090827270617UL) + ((((uint64_t)op[6] * 7314918444598304289UL) + ((uint64_t)op[7] * 15010174984106496879UL) + ((uint64_t)op[8] * 8695639198555464415UL) + ((uint64_t)op[9] * 415742185652332287UL) + ((uint64_t)op[10] * 14971330595751891135UL)) * 18446744073709551609);
	tmp_q[6] = ((uint64_t)op[0] * 14971330595751891135UL) + ((uint64_t)op[1] * 1199911456779291470UL) + ((uint64_t)op[2] * 16079450188461417004UL) + ((uint64_t)op[3] * 13104052030274545521UL) + ((uint64_t)op[4] * 1873511197073093991UL) + ((uint64_t)op[5] * 6771460835154525505UL) + ((uint64_t)op[6] * 9398408090827270617UL) + ((((uint64_t)op[7] * 7314918444598304289UL) + ((uint64_t)op[8] * 15010174984106496879UL) + ((uint64_t)op[9] * 8695639198555464415UL) + ((uint64_t)op[10] * 415742185652332287UL)) * 18446744073709551609);
	tmp_q[7] = ((uint64_t)op[0] * 415742185652332287UL) + ((uint64_t)op[1] * 14971330595751891135UL) + ((uint64_t)op[2] * 1199911456779291470UL) + ((uint64_t)op[3] * 16079450188461417004UL) + ((uint64_t)op[4] * 13104052030274545521UL) + ((uint64_t)op[5] * 1873511197073093991UL) + ((uint64_t)op[6] * 6771460835154525505UL) + ((uint64_t)op[7] * 9398408090827270617UL) + ((((uint64_t)op[8] * 7314918444598304289UL) + ((uint64_t)op[9] * 15010174984106496879UL) + ((uint64_t)op[10] * 8695639198555464415UL)) * 18446744073709551609);
	tmp_q[8] = ((uint64_t)op[0] * 8695639198555464415UL) + ((uint64_t)op[1] * 415742185652332287UL) + ((uint64_t)op[2] * 14971330595751891135UL) + ((uint64_t)op[3] * 1199911456779291470UL) + ((uint64_t)op[4] * 16079450188461417004UL) + ((uint64_t)op[5] * 13104052030274545521UL) + ((uint64_t)op[6] * 1873511197073093991UL) + ((uint64_t)op[7] * 6771460835154525505UL) + ((uint64_t)op[8] * 9398408090827270617UL) + ((((uint64_t)op[9] * 7314918444598304289UL) + ((uint64_t)op[10] * 15010174984106496879UL)) * 18446744073709551609);
	tmp_q[9] = ((uint64_t)op[0] * 15010174984106496879UL) + ((uint64_t)op[1] * 8695639198555464415UL) + ((uint64_t)op[2] * 415742185652332287UL) + ((uint64_t)op[3] * 14971330595751891135UL) + ((uint64_t)op[4] * 1199911456779291470UL) + ((uint64_t)op[5] * 16079450188461417004UL) + ((uint64_t)op[6] * 13104052030274545521UL) + ((uint64_t)op[7] * 1873511197073093991UL) + ((uint64_t)op[8] * 6771460835154525505UL) + ((uint64_t)op[9] * 9398408090827270617UL) + ((uint64_t)op[10] * 4135803108940524825UL);
	tmp_q[10] = ((uint64_t)op[0] * 7314918444598304289UL) + ((uint64_t)op[1] * 15010174984106496879UL) + ((uint64_t)op[2] * 8695639198555464415UL) + ((uint64_t)op[3] * 415742185652332287UL) + ((uint64_t)op[4] * 14971330595751891135UL) + ((uint64_t)op[5] * 1199911456779291470UL) + ((uint64_t)op[6] * 16079450188461417004UL) + ((uint64_t)op[7] * 13104052030274545521UL) + ((uint64_t)op[8] * 1873511197073093991UL) + ((uint64_t)op[9] * 6771460835154525505UL) + ((uint64_t)op[10] * 9398408090827270617UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 15718637948236L) - ((-((int128)tmp_q[1] * 22722399257249L) - ((int128)tmp_q[2] * 82295024771034L) - ((int128)tmp_q[3] * 36464659919831L) - ((int128)tmp_q[4] * 10833866126300L) - ((int128)tmp_q[5] * 42287080185528L) - ((int128)tmp_q[6] * 34147093906685L) + ((int128)tmp_q[7] * 2165000142572L) - ((int128)tmp_q[8] * 23977602794247L) + ((int128)tmp_q[9] * 7494427281824L) + ((int128)tmp_q[10] * 86434688149905L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 86434688149905L) + ((int128)tmp_q[1] * 15718637948236L) - ((-((int128)tmp_q[2] * 22722399257249L) - ((int128)tmp_q[3] * 82295024771034L) - ((int128)tmp_q[4] * 36464659919831L) - ((int128)tmp_q[5] * 10833866126300L) - ((int128)tmp_q[6] * 42287080185528L) - ((int128)tmp_q[7] * 34147093906685L) + ((int128)tmp_q[8] * 2165000142572L) - ((int128)tmp_q[9] * 23977602794247L) + ((int128)tmp_q[10] * 7494427281824L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 7494427281824L) + ((int128)tmp_q[1] * 86434688149905L) + ((int128)tmp_q[2] * 15718637948236L) - ((-((int128)tmp_q[3] * 22722399257249L) - ((int128)tmp_q[4] * 82295024771034L) - ((int128)tmp_q[5] * 36464659919831L) - ((int128)tmp_q[6] * 10833866126300L) - ((int128)tmp_q[7] * 42287080185528L) - ((int128)tmp_q[8] * 34147093906685L) + ((int128)tmp_q[9] * 2165000142572L) - ((int128)tmp_q[10] * 23977602794247L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 23977602794247L) + ((int128)tmp_q[1] * 7494427281824L) + ((int128)tmp_q[2] * 86434688149905L) + ((int128)tmp_q[3] * 15718637948236L) - ((-((int128)tmp_q[4] * 22722399257249L) - ((int128)tmp_q[5] * 82295024771034L) - ((int128)tmp_q[6] * 36464659919831L) - ((int128)tmp_q[7] * 10833866126300L) - ((int128)tmp_q[8] * 42287080185528L) - ((int128)tmp_q[9] * 34147093906685L) + ((int128)tmp_q[10] * 2165000142572L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 2165000142572L) - ((int128)tmp_q[1] * 23977602794247L) + ((int128)tmp_q[2] * 7494427281824L) + ((int128)tmp_q[3] * 86434688149905L) + ((int128)tmp_q[4] * 15718637948236L) - ((-((int128)tmp_q[5] * 22722399257249L) - ((int128)tmp_q[6] * 82295024771034L) - ((int128)tmp_q[7] * 36464659919831L) - ((int128)tmp_q[8] * 10833866126300L) - ((int128)tmp_q[9] * 42287080185528L) - ((int128)tmp_q[10] * 34147093906685L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 34147093906685L) + ((int128)tmp_q[1] * 2165000142572L) - ((int128)tmp_q[2] * 23977602794247L) + ((int128)tmp_q[3] * 7494427281824L) + ((int128)tmp_q[4] * 86434688149905L) + ((int128)tmp_q[5] * 15718637948236L) - ((-((int128)tmp_q[6] * 22722399257249L) - ((int128)tmp_q[7] * 82295024771034L) - ((int128)tmp_q[8] * 36464659919831L) - ((int128)tmp_q[9] * 10833866126300L) - ((int128)tmp_q[10] * 42287080185528L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 42287080185528L) - ((int128)tmp_q[1] * 34147093906685L) + ((int128)tmp_q[2] * 2165000142572L) - ((int128)tmp_q[3] * 23977602794247L) + ((int128)tmp_q[4] * 7494427281824L) + ((int128)tmp_q[5] * 86434688149905L) + ((int128)tmp_q[6] * 15718637948236L) - ((-((int128)tmp_q[7] * 22722399257249L) - ((int128)tmp_q[8] * 82295024771034L) - ((int128)tmp_q[9] * 36464659919831L) - ((int128)tmp_q[10] * 10833866126300L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 10833866126300L) - ((int128)tmp_q[1] * 42287080185528L) - ((int128)tmp_q[2] * 34147093906685L) + ((int128)tmp_q[3] * 2165000142572L) - ((int128)tmp_q[4] * 23977602794247L) + ((int128)tmp_q[5] * 7494427281824L) + ((int128)tmp_q[6] * 86434688149905L) + ((int128)tmp_q[7] * 15718637948236L) - ((-((int128)tmp_q[8] * 22722399257249L) - ((int128)tmp_q[9] * 82295024771034L) - ((int128)tmp_q[10] * 36464659919831L)) * 7);
	tmp_zero[8] = -((int128)tmp_q[0] * 36464659919831L) - ((int128)tmp_q[1] * 10833866126300L) - ((int128)tmp_q[2] * 42287080185528L) - ((int128)tmp_q[3] * 34147093906685L) + ((int128)tmp_q[4] * 2165000142572L) - ((int128)tmp_q[5] * 23977602794247L) + ((int128)tmp_q[6] * 7494427281824L) + ((int128)tmp_q[7] * 86434688149905L) + ((int128)tmp_q[8] * 15718637948236L) - ((-((int128)tmp_q[9] * 22722399257249L) - ((int128)tmp_q[10] * 82295024771034L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 82295024771034L) - ((int128)tmp_q[1] * 36464659919831L) - ((int128)tmp_q[2] * 10833866126300L) - ((int128)tmp_q[3] * 42287080185528L) - ((int128)tmp_q[4] * 34147093906685L) + ((int128)tmp_q[5] * 2165000142572L) - ((int128)tmp_q[6] * 23977602794247L) + ((int128)tmp_q[7] * 7494427281824L) + ((int128)tmp_q[8] * 86434688149905L) + ((int128)tmp_q[9] * 15718637948236L) + ((int128)tmp_q[10] * 159056794800743L);
	tmp_zero[10] = -((int128)tmp_q[0] * 22722399257249L) - ((int128)tmp_q[1] * 82295024771034L) - ((int128)tmp_q[2] * 36464659919831L) - ((int128)tmp_q[3] * 10833866126300L) - ((int128)tmp_q[4] * 42287080185528L) - ((int128)tmp_q[5] * 34147093906685L) + ((int128)tmp_q[6] * 2165000142572L) - ((int128)tmp_q[7] * 23977602794247L) + ((int128)tmp_q[8] * 7494427281824L) + ((int128)tmp_q[9] * 86434688149905L) + ((int128)tmp_q[10] * 15718637948236L);

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

