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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) << 4);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10063660970701225915UL) + ((((uint64_t)op[1] * 13979567360937708182UL) + ((uint64_t)op[2] * 9292056772890500392UL) + ((uint64_t)op[3] * 12859200970704531568UL) + ((uint64_t)op[4] * 16764371915592266321UL) + ((uint64_t)op[5] * 8350102646421872521UL) + ((uint64_t)op[6] * 13473507727689569943UL) + ((uint64_t)op[7] * 4777879096736719616UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 4777879096736719616UL) + ((uint64_t)op[1] * 10063660970701225915UL) + ((((uint64_t)op[2] * 13979567360937708182UL) + ((uint64_t)op[3] * 9292056772890500392UL) + ((uint64_t)op[4] * 12859200970704531568UL) + ((uint64_t)op[5] * 16764371915592266321UL) + ((uint64_t)op[6] * 8350102646421872521UL) + ((uint64_t)op[7] * 13473507727689569943UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 13473507727689569943UL) + ((uint64_t)op[1] * 4777879096736719616UL) + ((uint64_t)op[2] * 10063660970701225915UL) + ((((uint64_t)op[3] * 13979567360937708182UL) + ((uint64_t)op[4] * 9292056772890500392UL) + ((uint64_t)op[5] * 12859200970704531568UL) + ((uint64_t)op[6] * 16764371915592266321UL) + ((uint64_t)op[7] * 8350102646421872521UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 8350102646421872521UL) + ((uint64_t)op[1] * 13473507727689569943UL) + ((uint64_t)op[2] * 4777879096736719616UL) + ((uint64_t)op[3] * 10063660970701225915UL) + ((((uint64_t)op[4] * 13979567360937708182UL) + ((uint64_t)op[5] * 9292056772890500392UL) + ((uint64_t)op[6] * 12859200970704531568UL) + ((uint64_t)op[7] * 16764371915592266321UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 16764371915592266321UL) + ((uint64_t)op[1] * 8350102646421872521UL) + ((uint64_t)op[2] * 13473507727689569943UL) + ((uint64_t)op[3] * 4777879096736719616UL) + ((uint64_t)op[4] * 10063660970701225915UL) + ((((uint64_t)op[5] * 13979567360937708182UL) + ((uint64_t)op[6] * 9292056772890500392UL) + ((uint64_t)op[7] * 12859200970704531568UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 12859200970704531568UL) + ((uint64_t)op[1] * 16764371915592266321UL) + ((uint64_t)op[2] * 8350102646421872521UL) + ((uint64_t)op[3] * 13473507727689569943UL) + ((uint64_t)op[4] * 4777879096736719616UL) + ((uint64_t)op[5] * 10063660970701225915UL) + ((((uint64_t)op[6] * 13979567360937708182UL) + ((uint64_t)op[7] * 9292056772890500392UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 9292056772890500392UL) + ((uint64_t)op[1] * 12859200970704531568UL) + ((uint64_t)op[2] * 16764371915592266321UL) + ((uint64_t)op[3] * 8350102646421872521UL) + ((uint64_t)op[4] * 13473507727689569943UL) + ((uint64_t)op[5] * 4777879096736719616UL) + ((uint64_t)op[6] * 10063660970701225915UL) + ((uint64_t)op[7] * 17290669628465195856UL);
	tmp_q[7] = ((uint64_t)op[0] * 13979567360937708182UL) + ((uint64_t)op[1] * 9292056772890500392UL) + ((uint64_t)op[2] * 12859200970704531568UL) + ((uint64_t)op[3] * 16764371915592266321UL) + ((uint64_t)op[4] * 8350102646421872521UL) + ((uint64_t)op[5] * 13473507727689569943UL) + ((uint64_t)op[6] * 4777879096736719616UL) + ((uint64_t)op[7] * 10063660970701225915UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 95073189743747L) - ((-((int128)tmp_q[1] * 81082126359301L) - ((int128)tmp_q[2] * 11437739585598L) + ((int128)tmp_q[3] * 160963949991286L) - ((int128)tmp_q[4] * 3061087868674L) - ((int128)tmp_q[5] * 107708556508607L) + ((int128)tmp_q[6] * 78923077378567L) - ((int128)tmp_q[7] * 130945907213912L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 130945907213912L) - ((int128)tmp_q[1] * 95073189743747L) - ((-((int128)tmp_q[2] * 81082126359301L) - ((int128)tmp_q[3] * 11437739585598L) + ((int128)tmp_q[4] * 160963949991286L) - ((int128)tmp_q[5] * 3061087868674L) - ((int128)tmp_q[6] * 107708556508607L) + ((int128)tmp_q[7] * 78923077378567L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 78923077378567L) - ((int128)tmp_q[1] * 130945907213912L) - ((int128)tmp_q[2] * 95073189743747L) - ((-((int128)tmp_q[3] * 81082126359301L) - ((int128)tmp_q[4] * 11437739585598L) + ((int128)tmp_q[5] * 160963949991286L) - ((int128)tmp_q[6] * 3061087868674L) - ((int128)tmp_q[7] * 107708556508607L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 107708556508607L) + ((int128)tmp_q[1] * 78923077378567L) - ((int128)tmp_q[2] * 130945907213912L) - ((int128)tmp_q[3] * 95073189743747L) - ((-((int128)tmp_q[4] * 81082126359301L) - ((int128)tmp_q[5] * 11437739585598L) + ((int128)tmp_q[6] * 160963949991286L) - ((int128)tmp_q[7] * 3061087868674L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 3061087868674L) - ((int128)tmp_q[1] * 107708556508607L) + ((int128)tmp_q[2] * 78923077378567L) - ((int128)tmp_q[3] * 130945907213912L) - ((int128)tmp_q[4] * 95073189743747L) - ((-((int128)tmp_q[5] * 81082126359301L) - ((int128)tmp_q[6] * 11437739585598L) + ((int128)tmp_q[7] * 160963949991286L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 160963949991286L) - ((int128)tmp_q[1] * 3061087868674L) - ((int128)tmp_q[2] * 107708556508607L) + ((int128)tmp_q[3] * 78923077378567L) - ((int128)tmp_q[4] * 130945907213912L) - ((int128)tmp_q[5] * 95073189743747L) - ((-((int128)tmp_q[6] * 81082126359301L) - ((int128)tmp_q[7] * 11437739585598L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 11437739585598L) + ((int128)tmp_q[1] * 160963949991286L) - ((int128)tmp_q[2] * 3061087868674L) - ((int128)tmp_q[3] * 107708556508607L) + ((int128)tmp_q[4] * 78923077378567L) - ((int128)tmp_q[5] * 130945907213912L) - ((int128)tmp_q[6] * 95073189743747L) + ((int128)tmp_q[7] * 648657010874408L);
	tmp_zero[7] = -((int128)tmp_q[0] * 81082126359301L) - ((int128)tmp_q[1] * 11437739585598L) + ((int128)tmp_q[2] * 160963949991286L) - ((int128)tmp_q[3] * 3061087868674L) - ((int128)tmp_q[4] * 107708556508607L) + ((int128)tmp_q[5] * 78923077378567L) - ((int128)tmp_q[6] * 130945907213912L) - ((int128)tmp_q[7] * 95073189743747L);

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

