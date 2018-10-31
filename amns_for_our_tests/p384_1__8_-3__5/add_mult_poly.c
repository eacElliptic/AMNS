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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7026126602821087541UL) + ((((uint64_t)op[1] * 4697535196755499310UL) + ((uint64_t)op[2] * 8238100885460117944UL) + ((uint64_t)op[3] * 11088148055710922426UL) + ((uint64_t)op[4] * 703632857184084734UL) + ((uint64_t)op[5] * 12013657049686004605UL) + ((uint64_t)op[6] * 9335170955979142009UL) + ((uint64_t)op[7] * 7147212304322312342UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 7147212304322312342UL) + ((uint64_t)op[1] * 7026126602821087541UL) + ((((uint64_t)op[2] * 4697535196755499310UL) + ((uint64_t)op[3] * 8238100885460117944UL) + ((uint64_t)op[4] * 11088148055710922426UL) + ((uint64_t)op[5] * 703632857184084734UL) + ((uint64_t)op[6] * 12013657049686004605UL) + ((uint64_t)op[7] * 9335170955979142009UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 9335170955979142009UL) + ((uint64_t)op[1] * 7147212304322312342UL) + ((uint64_t)op[2] * 7026126602821087541UL) + ((((uint64_t)op[3] * 4697535196755499310UL) + ((uint64_t)op[4] * 8238100885460117944UL) + ((uint64_t)op[5] * 11088148055710922426UL) + ((uint64_t)op[6] * 703632857184084734UL) + ((uint64_t)op[7] * 12013657049686004605UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 12013657049686004605UL) + ((uint64_t)op[1] * 9335170955979142009UL) + ((uint64_t)op[2] * 7147212304322312342UL) + ((uint64_t)op[3] * 7026126602821087541UL) + ((((uint64_t)op[4] * 4697535196755499310UL) + ((uint64_t)op[5] * 8238100885460117944UL) + ((uint64_t)op[6] * 11088148055710922426UL) + ((uint64_t)op[7] * 703632857184084734UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 703632857184084734UL) + ((uint64_t)op[1] * 12013657049686004605UL) + ((uint64_t)op[2] * 9335170955979142009UL) + ((uint64_t)op[3] * 7147212304322312342UL) + ((uint64_t)op[4] * 7026126602821087541UL) + ((((uint64_t)op[5] * 4697535196755499310UL) + ((uint64_t)op[6] * 8238100885460117944UL) + ((uint64_t)op[7] * 11088148055710922426UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 11088148055710922426UL) + ((uint64_t)op[1] * 703632857184084734UL) + ((uint64_t)op[2] * 12013657049686004605UL) + ((uint64_t)op[3] * 9335170955979142009UL) + ((uint64_t)op[4] * 7147212304322312342UL) + ((uint64_t)op[5] * 7026126602821087541UL) + ((((uint64_t)op[6] * 4697535196755499310UL) + ((uint64_t)op[7] * 8238100885460117944UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 8238100885460117944UL) + ((uint64_t)op[1] * 11088148055710922426UL) + ((uint64_t)op[2] * 703632857184084734UL) + ((uint64_t)op[3] * 12013657049686004605UL) + ((uint64_t)op[4] * 9335170955979142009UL) + ((uint64_t)op[5] * 7147212304322312342UL) + ((uint64_t)op[6] * 7026126602821087541UL) + ((uint64_t)op[7] * 4354138483443053686UL);
	tmp_q[7] = ((uint64_t)op[0] * 4697535196755499310UL) + ((uint64_t)op[1] * 8238100885460117944UL) + ((uint64_t)op[2] * 11088148055710922426UL) + ((uint64_t)op[3] * 703632857184084734UL) + ((uint64_t)op[4] * 12013657049686004605UL) + ((uint64_t)op[5] * 9335170955979142009UL) + ((uint64_t)op[6] * 7147212304322312342UL) + ((uint64_t)op[7] * 7026126602821087541UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 37619074256863L) - ((-((int128)tmp_q[1] * 117198226354665L) + ((int128)tmp_q[2] * 84753094924003L) + ((int128)tmp_q[3] * 40728271448029L) - ((int128)tmp_q[4] * 109515759211526L) + ((int128)tmp_q[5] * 61029722294103L) + ((int128)tmp_q[6] * 23905816200512L) - ((int128)tmp_q[7] * 117424359647530L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 117424359647530L) + ((int128)tmp_q[1] * 37619074256863L) - ((-((int128)tmp_q[2] * 117198226354665L) + ((int128)tmp_q[3] * 84753094924003L) + ((int128)tmp_q[4] * 40728271448029L) - ((int128)tmp_q[5] * 109515759211526L) + ((int128)tmp_q[6] * 61029722294103L) + ((int128)tmp_q[7] * 23905816200512L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 23905816200512L) - ((int128)tmp_q[1] * 117424359647530L) + ((int128)tmp_q[2] * 37619074256863L) - ((-((int128)tmp_q[3] * 117198226354665L) + ((int128)tmp_q[4] * 84753094924003L) + ((int128)tmp_q[5] * 40728271448029L) - ((int128)tmp_q[6] * 109515759211526L) + ((int128)tmp_q[7] * 61029722294103L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 61029722294103L) + ((int128)tmp_q[1] * 23905816200512L) - ((int128)tmp_q[2] * 117424359647530L) + ((int128)tmp_q[3] * 37619074256863L) - ((-((int128)tmp_q[4] * 117198226354665L) + ((int128)tmp_q[5] * 84753094924003L) + ((int128)tmp_q[6] * 40728271448029L) - ((int128)tmp_q[7] * 109515759211526L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 109515759211526L) + ((int128)tmp_q[1] * 61029722294103L) + ((int128)tmp_q[2] * 23905816200512L) - ((int128)tmp_q[3] * 117424359647530L) + ((int128)tmp_q[4] * 37619074256863L) - ((-((int128)tmp_q[5] * 117198226354665L) + ((int128)tmp_q[6] * 84753094924003L) + ((int128)tmp_q[7] * 40728271448029L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 40728271448029L) - ((int128)tmp_q[1] * 109515759211526L) + ((int128)tmp_q[2] * 61029722294103L) + ((int128)tmp_q[3] * 23905816200512L) - ((int128)tmp_q[4] * 117424359647530L) + ((int128)tmp_q[5] * 37619074256863L) - ((-((int128)tmp_q[6] * 117198226354665L) + ((int128)tmp_q[7] * 84753094924003L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 84753094924003L) + ((int128)tmp_q[1] * 40728271448029L) - ((int128)tmp_q[2] * 109515759211526L) + ((int128)tmp_q[3] * 61029722294103L) + ((int128)tmp_q[4] * 23905816200512L) - ((int128)tmp_q[5] * 117424359647530L) + ((int128)tmp_q[6] * 37619074256863L) + ((int128)tmp_q[7] * 351594679063995L);
	tmp_zero[7] = -((int128)tmp_q[0] * 117198226354665L) + ((int128)tmp_q[1] * 84753094924003L) + ((int128)tmp_q[2] * 40728271448029L) - ((int128)tmp_q[3] * 109515759211526L) + ((int128)tmp_q[4] * 61029722294103L) + ((int128)tmp_q[5] * 23905816200512L) - ((int128)tmp_q[6] * 117424359647530L) + ((int128)tmp_q[7] * 37619074256863L);

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

