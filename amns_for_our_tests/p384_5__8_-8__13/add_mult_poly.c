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
	tmp_q[0] = ((uint64_t)op[0] * 3337648998760853193UL) + ((((uint64_t)op[1] * 7227717734689569545UL) + ((uint64_t)op[2] * 6877501086195889675UL) + ((uint64_t)op[3] * 13966677085885815708UL) + ((uint64_t)op[4] * 16170750529072111897UL) + ((uint64_t)op[5] * 10008440245622709940UL) + ((uint64_t)op[6] * 12033635644162101163UL) + ((uint64_t)op[7] * 967694480168231465UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 967694480168231465UL) + ((uint64_t)op[1] * 3337648998760853193UL) + ((((uint64_t)op[2] * 7227717734689569545UL) + ((uint64_t)op[3] * 6877501086195889675UL) + ((uint64_t)op[4] * 13966677085885815708UL) + ((uint64_t)op[5] * 16170750529072111897UL) + ((uint64_t)op[6] * 10008440245622709940UL) + ((uint64_t)op[7] * 12033635644162101163UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 12033635644162101163UL) + ((uint64_t)op[1] * 967694480168231465UL) + ((uint64_t)op[2] * 3337648998760853193UL) + ((((uint64_t)op[3] * 7227717734689569545UL) + ((uint64_t)op[4] * 6877501086195889675UL) + ((uint64_t)op[5] * 13966677085885815708UL) + ((uint64_t)op[6] * 16170750529072111897UL) + ((uint64_t)op[7] * 10008440245622709940UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 10008440245622709940UL) + ((uint64_t)op[1] * 12033635644162101163UL) + ((uint64_t)op[2] * 967694480168231465UL) + ((uint64_t)op[3] * 3337648998760853193UL) + ((((uint64_t)op[4] * 7227717734689569545UL) + ((uint64_t)op[5] * 6877501086195889675UL) + ((uint64_t)op[6] * 13966677085885815708UL) + ((uint64_t)op[7] * 16170750529072111897UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 16170750529072111897UL) + ((uint64_t)op[1] * 10008440245622709940UL) + ((uint64_t)op[2] * 12033635644162101163UL) + ((uint64_t)op[3] * 967694480168231465UL) + ((uint64_t)op[4] * 3337648998760853193UL) + ((((uint64_t)op[5] * 7227717734689569545UL) + ((uint64_t)op[6] * 6877501086195889675UL) + ((uint64_t)op[7] * 13966677085885815708UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 13966677085885815708UL) + ((uint64_t)op[1] * 16170750529072111897UL) + ((uint64_t)op[2] * 10008440245622709940UL) + ((uint64_t)op[3] * 12033635644162101163UL) + ((uint64_t)op[4] * 967694480168231465UL) + ((uint64_t)op[5] * 3337648998760853193UL) + ((((uint64_t)op[6] * 7227717734689569545UL) + ((uint64_t)op[7] * 6877501086195889675UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 6877501086195889675UL) + ((uint64_t)op[1] * 13966677085885815708UL) + ((uint64_t)op[2] * 16170750529072111897UL) + ((uint64_t)op[3] * 10008440245622709940UL) + ((uint64_t)op[4] * 12033635644162101163UL) + ((uint64_t)op[5] * 967694480168231465UL) + ((uint64_t)op[6] * 3337648998760853193UL) + ((uint64_t)op[7] * 15965234417321650104UL);
	tmp_q[7] = ((uint64_t)op[0] * 7227717734689569545UL) + ((uint64_t)op[1] * 6877501086195889675UL) + ((uint64_t)op[2] * 13966677085885815708UL) + ((uint64_t)op[3] * 16170750529072111897UL) + ((uint64_t)op[4] * 10008440245622709940UL) + ((uint64_t)op[5] * 12033635644162101163UL) + ((uint64_t)op[6] * 967694480168231465UL) + ((uint64_t)op[7] * 3337648998760853193UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 137664503339319L) - ((-((int128)tmp_q[1] * 23380311911126L) + ((int128)tmp_q[2] * 134296793682691L) + ((int128)tmp_q[3] * 33757614415070L) - ((int128)tmp_q[4] * 117921626065208L) - ((int128)tmp_q[5] * 89981360863737L) - ((int128)tmp_q[6] * 24399926177414L) - ((int128)tmp_q[7] * 95082402874143L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 95082402874143L) + ((int128)tmp_q[1] * 137664503339319L) - ((-((int128)tmp_q[2] * 23380311911126L) + ((int128)tmp_q[3] * 134296793682691L) + ((int128)tmp_q[4] * 33757614415070L) - ((int128)tmp_q[5] * 117921626065208L) - ((int128)tmp_q[6] * 89981360863737L) - ((int128)tmp_q[7] * 24399926177414L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 24399926177414L) - ((int128)tmp_q[1] * 95082402874143L) + ((int128)tmp_q[2] * 137664503339319L) - ((-((int128)tmp_q[3] * 23380311911126L) + ((int128)tmp_q[4] * 134296793682691L) + ((int128)tmp_q[5] * 33757614415070L) - ((int128)tmp_q[6] * 117921626065208L) - ((int128)tmp_q[7] * 89981360863737L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 89981360863737L) - ((int128)tmp_q[1] * 24399926177414L) - ((int128)tmp_q[2] * 95082402874143L) + ((int128)tmp_q[3] * 137664503339319L) - ((-((int128)tmp_q[4] * 23380311911126L) + ((int128)tmp_q[5] * 134296793682691L) + ((int128)tmp_q[6] * 33757614415070L) - ((int128)tmp_q[7] * 117921626065208L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 117921626065208L) - ((int128)tmp_q[1] * 89981360863737L) - ((int128)tmp_q[2] * 24399926177414L) - ((int128)tmp_q[3] * 95082402874143L) + ((int128)tmp_q[4] * 137664503339319L) - ((-((int128)tmp_q[5] * 23380311911126L) + ((int128)tmp_q[6] * 134296793682691L) + ((int128)tmp_q[7] * 33757614415070L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 33757614415070L) - ((int128)tmp_q[1] * 117921626065208L) - ((int128)tmp_q[2] * 89981360863737L) - ((int128)tmp_q[3] * 24399926177414L) - ((int128)tmp_q[4] * 95082402874143L) + ((int128)tmp_q[5] * 137664503339319L) - ((-((int128)tmp_q[6] * 23380311911126L) + ((int128)tmp_q[7] * 134296793682691L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 134296793682691L) + ((int128)tmp_q[1] * 33757614415070L) - ((int128)tmp_q[2] * 117921626065208L) - ((int128)tmp_q[3] * 89981360863737L) - ((int128)tmp_q[4] * 24399926177414L) - ((int128)tmp_q[5] * 95082402874143L) + ((int128)tmp_q[6] * 137664503339319L) + ((int128)tmp_q[7] * 187042495289008L);
	tmp_zero[7] = -((int128)tmp_q[0] * 23380311911126L) + ((int128)tmp_q[1] * 134296793682691L) + ((int128)tmp_q[2] * 33757614415070L) - ((int128)tmp_q[3] * 117921626065208L) - ((int128)tmp_q[4] * 89981360863737L) - ((int128)tmp_q[5] * 24399926177414L) - ((int128)tmp_q[6] * 95082402874143L) + ((int128)tmp_q[7] * 137664503339319L);

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

