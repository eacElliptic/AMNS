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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17668037991674163917UL) + ((((uint64_t)op[1] * 15748325767709809803UL) + ((uint64_t)op[2] * 11416321829371204223UL) + ((uint64_t)op[3] * 4555312442794048614UL) + ((uint64_t)op[4] * 7662390325980935779UL) + ((uint64_t)op[5] * 4035293862314035958UL) + ((uint64_t)op[6] * 16219903615691092167UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 16219903615691092167UL) + ((uint64_t)op[1] * 17668037991674163917UL) + ((((uint64_t)op[2] * 15748325767709809803UL) + ((uint64_t)op[3] * 11416321829371204223UL) + ((uint64_t)op[4] * 4555312442794048614UL) + ((uint64_t)op[5] * 7662390325980935779UL) + ((uint64_t)op[6] * 4035293862314035958UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4035293862314035958UL) + ((uint64_t)op[1] * 16219903615691092167UL) + ((uint64_t)op[2] * 17668037991674163917UL) + ((((uint64_t)op[3] * 15748325767709809803UL) + ((uint64_t)op[4] * 11416321829371204223UL) + ((uint64_t)op[5] * 4555312442794048614UL) + ((uint64_t)op[6] * 7662390325980935779UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 7662390325980935779UL) + ((uint64_t)op[1] * 4035293862314035958UL) + ((uint64_t)op[2] * 16219903615691092167UL) + ((uint64_t)op[3] * 17668037991674163917UL) + ((((uint64_t)op[4] * 15748325767709809803UL) + ((uint64_t)op[5] * 11416321829371204223UL) + ((uint64_t)op[6] * 4555312442794048614UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 4555312442794048614UL) + ((uint64_t)op[1] * 7662390325980935779UL) + ((uint64_t)op[2] * 4035293862314035958UL) + ((uint64_t)op[3] * 16219903615691092167UL) + ((uint64_t)op[4] * 17668037991674163917UL) + ((((uint64_t)op[5] * 15748325767709809803UL) + ((uint64_t)op[6] * 11416321829371204223UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 11416321829371204223UL) + ((uint64_t)op[1] * 4555312442794048614UL) + ((uint64_t)op[2] * 7662390325980935779UL) + ((uint64_t)op[3] * 4035293862314035958UL) + ((uint64_t)op[4] * 16219903615691092167UL) + ((uint64_t)op[5] * 17668037991674163917UL) + ((uint64_t)op[6] * 5396836611999483626UL);
	tmp_q[6] = ((uint64_t)op[0] * 15748325767709809803UL) + ((uint64_t)op[1] * 11416321829371204223UL) + ((uint64_t)op[2] * 4555312442794048614UL) + ((uint64_t)op[3] * 7662390325980935779UL) + ((uint64_t)op[4] * 4035293862314035958UL) + ((uint64_t)op[5] * 16219903615691092167UL) + ((uint64_t)op[6] * 17668037991674163917UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 26896998749L) - ((-((int128)tmp_q[1] * 33899479163L) - ((int128)tmp_q[2] * 22250971507L) - ((int128)tmp_q[3] * 13309629139L) - ((int128)tmp_q[4] * 19012136598L) - ((int128)tmp_q[5] * 66203003891L) - ((int128)tmp_q[6] * 20685354105L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 20685354105L) - ((int128)tmp_q[1] * 26896998749L) - ((-((int128)tmp_q[2] * 33899479163L) - ((int128)tmp_q[3] * 22250971507L) - ((int128)tmp_q[4] * 13309629139L) - ((int128)tmp_q[5] * 19012136598L) - ((int128)tmp_q[6] * 66203003891L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 66203003891L) - ((int128)tmp_q[1] * 20685354105L) - ((int128)tmp_q[2] * 26896998749L) - ((-((int128)tmp_q[3] * 33899479163L) - ((int128)tmp_q[4] * 22250971507L) - ((int128)tmp_q[5] * 13309629139L) - ((int128)tmp_q[6] * 19012136598L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 19012136598L) - ((int128)tmp_q[1] * 66203003891L) - ((int128)tmp_q[2] * 20685354105L) - ((int128)tmp_q[3] * 26896998749L) - ((-((int128)tmp_q[4] * 33899479163L) - ((int128)tmp_q[5] * 22250971507L) - ((int128)tmp_q[6] * 13309629139L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 13309629139L) - ((int128)tmp_q[1] * 19012136598L) - ((int128)tmp_q[2] * 66203003891L) - ((int128)tmp_q[3] * 20685354105L) - ((int128)tmp_q[4] * 26896998749L) - ((-((int128)tmp_q[5] * 33899479163L) - ((int128)tmp_q[6] * 22250971507L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 22250971507L) - ((int128)tmp_q[1] * 13309629139L) - ((int128)tmp_q[2] * 19012136598L) - ((int128)tmp_q[3] * 66203003891L) - ((int128)tmp_q[4] * 20685354105L) - ((int128)tmp_q[5] * 26896998749L) + ((int128)tmp_q[6] * 67798958326L);
	tmp_zero[6] = -((int128)tmp_q[0] * 33899479163L) - ((int128)tmp_q[1] * 22250971507L) - ((int128)tmp_q[2] * 13309629139L) - ((int128)tmp_q[3] * 19012136598L) - ((int128)tmp_q[4] * 66203003891L) - ((int128)tmp_q[5] * 20685354105L) - ((int128)tmp_q[6] * 26896998749L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

