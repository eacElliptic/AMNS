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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2778164821557608479UL) + ((((uint64_t)op[1] * 13039887820147563247UL) + ((uint64_t)op[2] * 16011737580140636103UL) + ((uint64_t)op[3] * 7082985834773122653UL) + ((uint64_t)op[4] * 2590950582624587966UL) + ((uint64_t)op[5] * 9573456165025359943UL) + ((uint64_t)op[6] * 11610372857797856826UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 11610372857797856826UL) + ((uint64_t)op[1] * 2778164821557608479UL) + ((((uint64_t)op[2] * 13039887820147563247UL) + ((uint64_t)op[3] * 16011737580140636103UL) + ((uint64_t)op[4] * 7082985834773122653UL) + ((uint64_t)op[5] * 2590950582624587966UL) + ((uint64_t)op[6] * 9573456165025359943UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 9573456165025359943UL) + ((uint64_t)op[1] * 11610372857797856826UL) + ((uint64_t)op[2] * 2778164821557608479UL) + ((((uint64_t)op[3] * 13039887820147563247UL) + ((uint64_t)op[4] * 16011737580140636103UL) + ((uint64_t)op[5] * 7082985834773122653UL) + ((uint64_t)op[6] * 2590950582624587966UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 2590950582624587966UL) + ((uint64_t)op[1] * 9573456165025359943UL) + ((uint64_t)op[2] * 11610372857797856826UL) + ((uint64_t)op[3] * 2778164821557608479UL) + ((((uint64_t)op[4] * 13039887820147563247UL) + ((uint64_t)op[5] * 16011737580140636103UL) + ((uint64_t)op[6] * 7082985834773122653UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 7082985834773122653UL) + ((uint64_t)op[1] * 2590950582624587966UL) + ((uint64_t)op[2] * 9573456165025359943UL) + ((uint64_t)op[3] * 11610372857797856826UL) + ((uint64_t)op[4] * 2778164821557608479UL) + ((((uint64_t)op[5] * 13039887820147563247UL) + ((uint64_t)op[6] * 16011737580140636103UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 16011737580140636103UL) + ((uint64_t)op[1] * 7082985834773122653UL) + ((uint64_t)op[2] * 2590950582624587966UL) + ((uint64_t)op[3] * 9573456165025359943UL) + ((uint64_t)op[4] * 11610372857797856826UL) + ((uint64_t)op[5] * 2778164821557608479UL) + ((uint64_t)op[6] * 3180680940538401860UL);
	tmp_q[6] = ((uint64_t)op[0] * 13039887820147563247UL) + ((uint64_t)op[1] * 16011737580140636103UL) + ((uint64_t)op[2] * 7082985834773122653UL) + ((uint64_t)op[3] * 2590950582624587966UL) + ((uint64_t)op[4] * 9573456165025359943UL) + ((uint64_t)op[5] * 11610372857797856826UL) + ((uint64_t)op[6] * 2778164821557608479UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 31792774401L) - ((-((int128)tmp_q[1] * 6181916200L) - ((int128)tmp_q[2] * 32132807983L) + ((int128)tmp_q[3] * 61679716962L) - ((int128)tmp_q[4] * 55837479326L) - ((int128)tmp_q[5] * 442226545L) + ((int128)tmp_q[6] * 37463788222L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 37463788222L) + ((int128)tmp_q[1] * 31792774401L) - ((-((int128)tmp_q[2] * 6181916200L) - ((int128)tmp_q[3] * 32132807983L) + ((int128)tmp_q[4] * 61679716962L) - ((int128)tmp_q[5] * 55837479326L) - ((int128)tmp_q[6] * 442226545L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 442226545L) + ((int128)tmp_q[1] * 37463788222L) + ((int128)tmp_q[2] * 31792774401L) - ((-((int128)tmp_q[3] * 6181916200L) - ((int128)tmp_q[4] * 32132807983L) + ((int128)tmp_q[5] * 61679716962L) - ((int128)tmp_q[6] * 55837479326L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 55837479326L) - ((int128)tmp_q[1] * 442226545L) + ((int128)tmp_q[2] * 37463788222L) + ((int128)tmp_q[3] * 31792774401L) - ((-((int128)tmp_q[4] * 6181916200L) - ((int128)tmp_q[5] * 32132807983L) + ((int128)tmp_q[6] * 61679716962L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 61679716962L) - ((int128)tmp_q[1] * 55837479326L) - ((int128)tmp_q[2] * 442226545L) + ((int128)tmp_q[3] * 37463788222L) + ((int128)tmp_q[4] * 31792774401L) - ((-((int128)tmp_q[5] * 6181916200L) - ((int128)tmp_q[6] * 32132807983L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 32132807983L) + ((int128)tmp_q[1] * 61679716962L) - ((int128)tmp_q[2] * 55837479326L) - ((int128)tmp_q[3] * 442226545L) + ((int128)tmp_q[4] * 37463788222L) + ((int128)tmp_q[5] * 31792774401L) + ((int128)tmp_q[6] * 24727664800L);
	tmp_zero[6] = -((int128)tmp_q[0] * 6181916200L) - ((int128)tmp_q[1] * 32132807983L) + ((int128)tmp_q[2] * 61679716962L) - ((int128)tmp_q[3] * 55837479326L) - ((int128)tmp_q[4] * 442226545L) + ((int128)tmp_q[5] * 37463788222L) + ((int128)tmp_q[6] * 31792774401L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

