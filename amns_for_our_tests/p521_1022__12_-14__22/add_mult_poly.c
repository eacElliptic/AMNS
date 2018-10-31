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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7669496064376994257UL) + ((((uint64_t)op[1] * 4042920637653388146UL) + ((uint64_t)op[2] * 11786028307832907115UL) + ((uint64_t)op[3] * 4598063758621063706UL) + ((uint64_t)op[4] * 8970259282742467521UL) + ((uint64_t)op[5] * 7083946420626452705UL) + ((uint64_t)op[6] * 15984750690151706086UL) + ((uint64_t)op[7] * 10003303695410170883UL) + ((uint64_t)op[8] * 1753511218199071413UL) + ((uint64_t)op[9] * 10279625485865553400UL) + ((uint64_t)op[10] * 10392849457081946955UL) + ((uint64_t)op[11] * 9400431246639801692UL)) * 18446744073709551602);
	tmp_q[1] = ((uint64_t)op[0] * 9400431246639801692UL) + ((uint64_t)op[1] * 7669496064376994257UL) + ((((uint64_t)op[2] * 4042920637653388146UL) + ((uint64_t)op[3] * 11786028307832907115UL) + ((uint64_t)op[4] * 4598063758621063706UL) + ((uint64_t)op[5] * 8970259282742467521UL) + ((uint64_t)op[6] * 7083946420626452705UL) + ((uint64_t)op[7] * 15984750690151706086UL) + ((uint64_t)op[8] * 10003303695410170883UL) + ((uint64_t)op[9] * 1753511218199071413UL) + ((uint64_t)op[10] * 10279625485865553400UL) + ((uint64_t)op[11] * 10392849457081946955UL)) * 18446744073709551602);
	tmp_q[2] = ((uint64_t)op[0] * 10392849457081946955UL) + ((uint64_t)op[1] * 9400431246639801692UL) + ((uint64_t)op[2] * 7669496064376994257UL) + ((((uint64_t)op[3] * 4042920637653388146UL) + ((uint64_t)op[4] * 11786028307832907115UL) + ((uint64_t)op[5] * 4598063758621063706UL) + ((uint64_t)op[6] * 8970259282742467521UL) + ((uint64_t)op[7] * 7083946420626452705UL) + ((uint64_t)op[8] * 15984750690151706086UL) + ((uint64_t)op[9] * 10003303695410170883UL) + ((uint64_t)op[10] * 1753511218199071413UL) + ((uint64_t)op[11] * 10279625485865553400UL)) * 18446744073709551602);
	tmp_q[3] = ((uint64_t)op[0] * 10279625485865553400UL) + ((uint64_t)op[1] * 10392849457081946955UL) + ((uint64_t)op[2] * 9400431246639801692UL) + ((uint64_t)op[3] * 7669496064376994257UL) + ((((uint64_t)op[4] * 4042920637653388146UL) + ((uint64_t)op[5] * 11786028307832907115UL) + ((uint64_t)op[6] * 4598063758621063706UL) + ((uint64_t)op[7] * 8970259282742467521UL) + ((uint64_t)op[8] * 7083946420626452705UL) + ((uint64_t)op[9] * 15984750690151706086UL) + ((uint64_t)op[10] * 10003303695410170883UL) + ((uint64_t)op[11] * 1753511218199071413UL)) * 18446744073709551602);
	tmp_q[4] = ((uint64_t)op[0] * 1753511218199071413UL) + ((uint64_t)op[1] * 10279625485865553400UL) + ((uint64_t)op[2] * 10392849457081946955UL) + ((uint64_t)op[3] * 9400431246639801692UL) + ((uint64_t)op[4] * 7669496064376994257UL) + ((((uint64_t)op[5] * 4042920637653388146UL) + ((uint64_t)op[6] * 11786028307832907115UL) + ((uint64_t)op[7] * 4598063758621063706UL) + ((uint64_t)op[8] * 8970259282742467521UL) + ((uint64_t)op[9] * 7083946420626452705UL) + ((uint64_t)op[10] * 15984750690151706086UL) + ((uint64_t)op[11] * 10003303695410170883UL)) * 18446744073709551602);
	tmp_q[5] = ((uint64_t)op[0] * 10003303695410170883UL) + ((uint64_t)op[1] * 1753511218199071413UL) + ((uint64_t)op[2] * 10279625485865553400UL) + ((uint64_t)op[3] * 10392849457081946955UL) + ((uint64_t)op[4] * 9400431246639801692UL) + ((uint64_t)op[5] * 7669496064376994257UL) + ((((uint64_t)op[6] * 4042920637653388146UL) + ((uint64_t)op[7] * 11786028307832907115UL) + ((uint64_t)op[8] * 4598063758621063706UL) + ((uint64_t)op[9] * 8970259282742467521UL) + ((uint64_t)op[10] * 7083946420626452705UL) + ((uint64_t)op[11] * 15984750690151706086UL)) * 18446744073709551602);
	tmp_q[6] = ((uint64_t)op[0] * 15984750690151706086UL) + ((uint64_t)op[1] * 10003303695410170883UL) + ((uint64_t)op[2] * 1753511218199071413UL) + ((uint64_t)op[3] * 10279625485865553400UL) + ((uint64_t)op[4] * 10392849457081946955UL) + ((uint64_t)op[5] * 9400431246639801692UL) + ((uint64_t)op[6] * 7669496064376994257UL) + ((((uint64_t)op[7] * 4042920637653388146UL) + ((uint64_t)op[8] * 11786028307832907115UL) + ((uint64_t)op[9] * 4598063758621063706UL) + ((uint64_t)op[10] * 8970259282742467521UL) + ((uint64_t)op[11] * 7083946420626452705UL)) * 18446744073709551602);
	tmp_q[7] = ((uint64_t)op[0] * 7083946420626452705UL) + ((uint64_t)op[1] * 15984750690151706086UL) + ((uint64_t)op[2] * 10003303695410170883UL) + ((uint64_t)op[3] * 1753511218199071413UL) + ((uint64_t)op[4] * 10279625485865553400UL) + ((uint64_t)op[5] * 10392849457081946955UL) + ((uint64_t)op[6] * 9400431246639801692UL) + ((uint64_t)op[7] * 7669496064376994257UL) + ((((uint64_t)op[8] * 4042920637653388146UL) + ((uint64_t)op[9] * 11786028307832907115UL) + ((uint64_t)op[10] * 4598063758621063706UL) + ((uint64_t)op[11] * 8970259282742467521UL)) * 18446744073709551602);
	tmp_q[8] = ((uint64_t)op[0] * 8970259282742467521UL) + ((uint64_t)op[1] * 7083946420626452705UL) + ((uint64_t)op[2] * 15984750690151706086UL) + ((uint64_t)op[3] * 10003303695410170883UL) + ((uint64_t)op[4] * 1753511218199071413UL) + ((uint64_t)op[5] * 10279625485865553400UL) + ((uint64_t)op[6] * 10392849457081946955UL) + ((uint64_t)op[7] * 9400431246639801692UL) + ((uint64_t)op[8] * 7669496064376994257UL) + ((((uint64_t)op[9] * 4042920637653388146UL) + ((uint64_t)op[10] * 11786028307832907115UL) + ((uint64_t)op[11] * 4598063758621063706UL)) * 18446744073709551602);
	tmp_q[9] = ((uint64_t)op[0] * 4598063758621063706UL) + ((uint64_t)op[1] * 8970259282742467521UL) + ((uint64_t)op[2] * 7083946420626452705UL) + ((uint64_t)op[3] * 15984750690151706086UL) + ((uint64_t)op[4] * 10003303695410170883UL) + ((uint64_t)op[5] * 1753511218199071413UL) + ((uint64_t)op[6] * 10279625485865553400UL) + ((uint64_t)op[7] * 10392849457081946955UL) + ((uint64_t)op[8] * 9400431246639801692UL) + ((uint64_t)op[9] * 7669496064376994257UL) + ((((uint64_t)op[10] * 4042920637653388146UL) + ((uint64_t)op[11] * 11786028307832907115UL)) * 18446744073709551602);
	tmp_q[10] = ((uint64_t)op[0] * 11786028307832907115UL) + ((uint64_t)op[1] * 4598063758621063706UL) + ((uint64_t)op[2] * 8970259282742467521UL) + ((uint64_t)op[3] * 7083946420626452705UL) + ((uint64_t)op[4] * 15984750690151706086UL) + ((uint64_t)op[5] * 10003303695410170883UL) + ((uint64_t)op[6] * 1753511218199071413UL) + ((uint64_t)op[7] * 10279625485865553400UL) + ((uint64_t)op[8] * 10392849457081946955UL) + ((uint64_t)op[9] * 9400431246639801692UL) + ((uint64_t)op[10] * 7669496064376994257UL) + ((uint64_t)op[11] * 17186087367690772420UL);
	tmp_q[11] = ((uint64_t)op[0] * 4042920637653388146UL) + ((uint64_t)op[1] * 11786028307832907115UL) + ((uint64_t)op[2] * 4598063758621063706UL) + ((uint64_t)op[3] * 8970259282742467521UL) + ((uint64_t)op[4] * 7083946420626452705UL) + ((uint64_t)op[5] * 15984750690151706086UL) + ((uint64_t)op[6] * 10003303695410170883UL) + ((uint64_t)op[7] * 1753511218199071413UL) + ((uint64_t)op[8] * 10279625485865553400UL) + ((uint64_t)op[9] * 10392849457081946955UL) + ((uint64_t)op[10] * 9400431246639801692UL) + ((uint64_t)op[11] * 7669496064376994257UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 260181053373L) - ((-((int128)tmp_q[1] * 3428502208085L) + ((int128)tmp_q[2] * 3035699839798L) - ((int128)tmp_q[3] * 888377295563L) - ((int128)tmp_q[4] * 1243912747944L) + ((int128)tmp_q[5] * 5737314091699L) + ((int128)tmp_q[6] * 5756691767429L) + ((int128)tmp_q[7] * 6883920369715L) + ((int128)tmp_q[8] * 7622497925556L) + ((int128)tmp_q[9] * 4266134365062L) - ((int128)tmp_q[10] * 242803750101L) - ((int128)tmp_q[11] * 6556381185396L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 6556381185396L) + ((int128)tmp_q[1] * 260181053373L) - ((-((int128)tmp_q[2] * 3428502208085L) + ((int128)tmp_q[3] * 3035699839798L) - ((int128)tmp_q[4] * 888377295563L) - ((int128)tmp_q[5] * 1243912747944L) + ((int128)tmp_q[6] * 5737314091699L) + ((int128)tmp_q[7] * 5756691767429L) + ((int128)tmp_q[8] * 6883920369715L) + ((int128)tmp_q[9] * 7622497925556L) + ((int128)tmp_q[10] * 4266134365062L) - ((int128)tmp_q[11] * 242803750101L)) * 14);
	tmp_zero[2] = -((int128)tmp_q[0] * 242803750101L) - ((int128)tmp_q[1] * 6556381185396L) + ((int128)tmp_q[2] * 260181053373L) - ((-((int128)tmp_q[3] * 3428502208085L) + ((int128)tmp_q[4] * 3035699839798L) - ((int128)tmp_q[5] * 888377295563L) - ((int128)tmp_q[6] * 1243912747944L) + ((int128)tmp_q[7] * 5737314091699L) + ((int128)tmp_q[8] * 5756691767429L) + ((int128)tmp_q[9] * 6883920369715L) + ((int128)tmp_q[10] * 7622497925556L) + ((int128)tmp_q[11] * 4266134365062L)) * 14);
	tmp_zero[3] = ((int128)tmp_q[0] * 4266134365062L) - ((int128)tmp_q[1] * 242803750101L) - ((int128)tmp_q[2] * 6556381185396L) + ((int128)tmp_q[3] * 260181053373L) - ((-((int128)tmp_q[4] * 3428502208085L) + ((int128)tmp_q[5] * 3035699839798L) - ((int128)tmp_q[6] * 888377295563L) - ((int128)tmp_q[7] * 1243912747944L) + ((int128)tmp_q[8] * 5737314091699L) + ((int128)tmp_q[9] * 5756691767429L) + ((int128)tmp_q[10] * 6883920369715L) + ((int128)tmp_q[11] * 7622497925556L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 7622497925556L) + ((int128)tmp_q[1] * 4266134365062L) - ((int128)tmp_q[2] * 242803750101L) - ((int128)tmp_q[3] * 6556381185396L) + ((int128)tmp_q[4] * 260181053373L) - ((-((int128)tmp_q[5] * 3428502208085L) + ((int128)tmp_q[6] * 3035699839798L) - ((int128)tmp_q[7] * 888377295563L) - ((int128)tmp_q[8] * 1243912747944L) + ((int128)tmp_q[9] * 5737314091699L) + ((int128)tmp_q[10] * 5756691767429L) + ((int128)tmp_q[11] * 6883920369715L)) * 14);
	tmp_zero[5] = ((int128)tmp_q[0] * 6883920369715L) + ((int128)tmp_q[1] * 7622497925556L) + ((int128)tmp_q[2] * 4266134365062L) - ((int128)tmp_q[3] * 242803750101L) - ((int128)tmp_q[4] * 6556381185396L) + ((int128)tmp_q[5] * 260181053373L) - ((-((int128)tmp_q[6] * 3428502208085L) + ((int128)tmp_q[7] * 3035699839798L) - ((int128)tmp_q[8] * 888377295563L) - ((int128)tmp_q[9] * 1243912747944L) + ((int128)tmp_q[10] * 5737314091699L) + ((int128)tmp_q[11] * 5756691767429L)) * 14);
	tmp_zero[6] = ((int128)tmp_q[0] * 5756691767429L) + ((int128)tmp_q[1] * 6883920369715L) + ((int128)tmp_q[2] * 7622497925556L) + ((int128)tmp_q[3] * 4266134365062L) - ((int128)tmp_q[4] * 242803750101L) - ((int128)tmp_q[5] * 6556381185396L) + ((int128)tmp_q[6] * 260181053373L) - ((-((int128)tmp_q[7] * 3428502208085L) + ((int128)tmp_q[8] * 3035699839798L) - ((int128)tmp_q[9] * 888377295563L) - ((int128)tmp_q[10] * 1243912747944L) + ((int128)tmp_q[11] * 5737314091699L)) * 14);
	tmp_zero[7] = ((int128)tmp_q[0] * 5737314091699L) + ((int128)tmp_q[1] * 5756691767429L) + ((int128)tmp_q[2] * 6883920369715L) + ((int128)tmp_q[3] * 7622497925556L) + ((int128)tmp_q[4] * 4266134365062L) - ((int128)tmp_q[5] * 242803750101L) - ((int128)tmp_q[6] * 6556381185396L) + ((int128)tmp_q[7] * 260181053373L) - ((-((int128)tmp_q[8] * 3428502208085L) + ((int128)tmp_q[9] * 3035699839798L) - ((int128)tmp_q[10] * 888377295563L) - ((int128)tmp_q[11] * 1243912747944L)) * 14);
	tmp_zero[8] = -((int128)tmp_q[0] * 1243912747944L) + ((int128)tmp_q[1] * 5737314091699L) + ((int128)tmp_q[2] * 5756691767429L) + ((int128)tmp_q[3] * 6883920369715L) + ((int128)tmp_q[4] * 7622497925556L) + ((int128)tmp_q[5] * 4266134365062L) - ((int128)tmp_q[6] * 242803750101L) - ((int128)tmp_q[7] * 6556381185396L) + ((int128)tmp_q[8] * 260181053373L) - ((-((int128)tmp_q[9] * 3428502208085L) + ((int128)tmp_q[10] * 3035699839798L) - ((int128)tmp_q[11] * 888377295563L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 888377295563L) - ((int128)tmp_q[1] * 1243912747944L) + ((int128)tmp_q[2] * 5737314091699L) + ((int128)tmp_q[3] * 5756691767429L) + ((int128)tmp_q[4] * 6883920369715L) + ((int128)tmp_q[5] * 7622497925556L) + ((int128)tmp_q[6] * 4266134365062L) - ((int128)tmp_q[7] * 242803750101L) - ((int128)tmp_q[8] * 6556381185396L) + ((int128)tmp_q[9] * 260181053373L) - ((-((int128)tmp_q[10] * 3428502208085L) + ((int128)tmp_q[11] * 3035699839798L)) * 14);
	tmp_zero[10] = ((int128)tmp_q[0] * 3035699839798L) - ((int128)tmp_q[1] * 888377295563L) - ((int128)tmp_q[2] * 1243912747944L) + ((int128)tmp_q[3] * 5737314091699L) + ((int128)tmp_q[4] * 5756691767429L) + ((int128)tmp_q[5] * 6883920369715L) + ((int128)tmp_q[6] * 7622497925556L) + ((int128)tmp_q[7] * 4266134365062L) - ((int128)tmp_q[8] * 242803750101L) - ((int128)tmp_q[9] * 6556381185396L) + ((int128)tmp_q[10] * 260181053373L) + ((int128)tmp_q[11] * 47999030913190L);
	tmp_zero[11] = -((int128)tmp_q[0] * 3428502208085L) + ((int128)tmp_q[1] * 3035699839798L) - ((int128)tmp_q[2] * 888377295563L) - ((int128)tmp_q[3] * 1243912747944L) + ((int128)tmp_q[4] * 5737314091699L) + ((int128)tmp_q[5] * 5756691767429L) + ((int128)tmp_q[6] * 6883920369715L) + ((int128)tmp_q[7] * 7622497925556L) + ((int128)tmp_q[8] * 4266134365062L) - ((int128)tmp_q[9] * 242803750101L) - ((int128)tmp_q[10] * 6556381185396L) + ((int128)tmp_q[11] * 260181053373L);

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

