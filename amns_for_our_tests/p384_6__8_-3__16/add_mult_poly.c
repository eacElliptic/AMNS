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
	tmp_q[0] = ((uint64_t)op[0] * 18227213348349103701UL) + ((((uint64_t)op[1] * 12215885458869462689UL) + ((uint64_t)op[2] * 11492138037089485889UL) + ((uint64_t)op[3] * 1463231905857229975UL) + ((uint64_t)op[4] * 10552420420846247146UL) + ((uint64_t)op[5] * 8980866077307913606UL) + ((uint64_t)op[6] * 9520839157442733196UL) + ((uint64_t)op[7] * 16697259823162269853UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16697259823162269853UL) + ((uint64_t)op[1] * 18227213348349103701UL) + ((((uint64_t)op[2] * 12215885458869462689UL) + ((uint64_t)op[3] * 11492138037089485889UL) + ((uint64_t)op[4] * 1463231905857229975UL) + ((uint64_t)op[5] * 10552420420846247146UL) + ((uint64_t)op[6] * 8980866077307913606UL) + ((uint64_t)op[7] * 9520839157442733196UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 9520839157442733196UL) + ((uint64_t)op[1] * 16697259823162269853UL) + ((uint64_t)op[2] * 18227213348349103701UL) + ((((uint64_t)op[3] * 12215885458869462689UL) + ((uint64_t)op[4] * 11492138037089485889UL) + ((uint64_t)op[5] * 1463231905857229975UL) + ((uint64_t)op[6] * 10552420420846247146UL) + ((uint64_t)op[7] * 8980866077307913606UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 8980866077307913606UL) + ((uint64_t)op[1] * 9520839157442733196UL) + ((uint64_t)op[2] * 16697259823162269853UL) + ((uint64_t)op[3] * 18227213348349103701UL) + ((((uint64_t)op[4] * 12215885458869462689UL) + ((uint64_t)op[5] * 11492138037089485889UL) + ((uint64_t)op[6] * 1463231905857229975UL) + ((uint64_t)op[7] * 10552420420846247146UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 10552420420846247146UL) + ((uint64_t)op[1] * 8980866077307913606UL) + ((uint64_t)op[2] * 9520839157442733196UL) + ((uint64_t)op[3] * 16697259823162269853UL) + ((uint64_t)op[4] * 18227213348349103701UL) + ((((uint64_t)op[5] * 12215885458869462689UL) + ((uint64_t)op[6] * 11492138037089485889UL) + ((uint64_t)op[7] * 1463231905857229975UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 1463231905857229975UL) + ((uint64_t)op[1] * 10552420420846247146UL) + ((uint64_t)op[2] * 8980866077307913606UL) + ((uint64_t)op[3] * 9520839157442733196UL) + ((uint64_t)op[4] * 16697259823162269853UL) + ((uint64_t)op[5] * 18227213348349103701UL) + ((((uint64_t)op[6] * 12215885458869462689UL) + ((uint64_t)op[7] * 11492138037089485889UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 11492138037089485889UL) + ((uint64_t)op[1] * 1463231905857229975UL) + ((uint64_t)op[2] * 10552420420846247146UL) + ((uint64_t)op[3] * 8980866077307913606UL) + ((uint64_t)op[4] * 9520839157442733196UL) + ((uint64_t)op[5] * 16697259823162269853UL) + ((uint64_t)op[6] * 18227213348349103701UL) + ((uint64_t)op[7] * 245831770810715165UL);
	tmp_q[7] = ((uint64_t)op[0] * 12215885458869462689UL) + ((uint64_t)op[1] * 11492138037089485889UL) + ((uint64_t)op[2] * 1463231905857229975UL) + ((uint64_t)op[3] * 10552420420846247146UL) + ((uint64_t)op[4] * 8980866077307913606UL) + ((uint64_t)op[5] * 9520839157442733196UL) + ((uint64_t)op[6] * 16697259823162269853UL) + ((uint64_t)op[7] * 18227213348349103701UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 66508233584156L) - ((((int128)tmp_q[1] * 59712631156562L) + ((int128)tmp_q[2] * 16333892610417L) - ((int128)tmp_q[3] * 125366480739526L) + ((int128)tmp_q[4] * 153651705679153L) - ((int128)tmp_q[5] * 127971554440120L) - ((int128)tmp_q[6] * 79389370587842L) + ((int128)tmp_q[7] * 46457761480413L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 46457761480413L) + ((int128)tmp_q[1] * 66508233584156L) - ((((int128)tmp_q[2] * 59712631156562L) + ((int128)tmp_q[3] * 16333892610417L) - ((int128)tmp_q[4] * 125366480739526L) + ((int128)tmp_q[5] * 153651705679153L) - ((int128)tmp_q[6] * 127971554440120L) - ((int128)tmp_q[7] * 79389370587842L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 79389370587842L) + ((int128)tmp_q[1] * 46457761480413L) + ((int128)tmp_q[2] * 66508233584156L) - ((((int128)tmp_q[3] * 59712631156562L) + ((int128)tmp_q[4] * 16333892610417L) - ((int128)tmp_q[5] * 125366480739526L) + ((int128)tmp_q[6] * 153651705679153L) - ((int128)tmp_q[7] * 127971554440120L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 127971554440120L) - ((int128)tmp_q[1] * 79389370587842L) + ((int128)tmp_q[2] * 46457761480413L) + ((int128)tmp_q[3] * 66508233584156L) - ((((int128)tmp_q[4] * 59712631156562L) + ((int128)tmp_q[5] * 16333892610417L) - ((int128)tmp_q[6] * 125366480739526L) + ((int128)tmp_q[7] * 153651705679153L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 153651705679153L) - ((int128)tmp_q[1] * 127971554440120L) - ((int128)tmp_q[2] * 79389370587842L) + ((int128)tmp_q[3] * 46457761480413L) + ((int128)tmp_q[4] * 66508233584156L) - ((((int128)tmp_q[5] * 59712631156562L) + ((int128)tmp_q[6] * 16333892610417L) - ((int128)tmp_q[7] * 125366480739526L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 125366480739526L) + ((int128)tmp_q[1] * 153651705679153L) - ((int128)tmp_q[2] * 127971554440120L) - ((int128)tmp_q[3] * 79389370587842L) + ((int128)tmp_q[4] * 46457761480413L) + ((int128)tmp_q[5] * 66508233584156L) - ((((int128)tmp_q[6] * 59712631156562L) + ((int128)tmp_q[7] * 16333892610417L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 16333892610417L) - ((int128)tmp_q[1] * 125366480739526L) + ((int128)tmp_q[2] * 153651705679153L) - ((int128)tmp_q[3] * 127971554440120L) - ((int128)tmp_q[4] * 79389370587842L) + ((int128)tmp_q[5] * 46457761480413L) + ((int128)tmp_q[6] * 66508233584156L) - ((int128)tmp_q[7] * 179137893469686L);
	tmp_zero[7] = ((int128)tmp_q[0] * 59712631156562L) + ((int128)tmp_q[1] * 16333892610417L) - ((int128)tmp_q[2] * 125366480739526L) + ((int128)tmp_q[3] * 153651705679153L) - ((int128)tmp_q[4] * 127971554440120L) - ((int128)tmp_q[5] * 79389370587842L) + ((int128)tmp_q[6] * 46457761480413L) + ((int128)tmp_q[7] * 66508233584156L);

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

