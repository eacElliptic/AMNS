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
	tmp_q[0] = ((uint64_t)op[0] * 6747879995015643189UL) + ((((uint64_t)op[1] * 7947569409371862095UL) + ((uint64_t)op[2] * 6566749008003890043UL) + ((uint64_t)op[3] * 16913353855706661744UL) + ((uint64_t)op[4] * 11045776113192476354UL) + ((uint64_t)op[5] * 7008775382811200779UL) + ((uint64_t)op[6] * 4854050199574680321UL) + ((uint64_t)op[7] * 9968235917052480426UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 9968235917052480426UL) + ((uint64_t)op[1] * 6747879995015643189UL) + ((((uint64_t)op[2] * 7947569409371862095UL) + ((uint64_t)op[3] * 6566749008003890043UL) + ((uint64_t)op[4] * 16913353855706661744UL) + ((uint64_t)op[5] * 11045776113192476354UL) + ((uint64_t)op[6] * 7008775382811200779UL) + ((uint64_t)op[7] * 4854050199574680321UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 4854050199574680321UL) + ((uint64_t)op[1] * 9968235917052480426UL) + ((uint64_t)op[2] * 6747879995015643189UL) + ((((uint64_t)op[3] * 7947569409371862095UL) + ((uint64_t)op[4] * 6566749008003890043UL) + ((uint64_t)op[5] * 16913353855706661744UL) + ((uint64_t)op[6] * 11045776113192476354UL) + ((uint64_t)op[7] * 7008775382811200779UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 7008775382811200779UL) + ((uint64_t)op[1] * 4854050199574680321UL) + ((uint64_t)op[2] * 9968235917052480426UL) + ((uint64_t)op[3] * 6747879995015643189UL) + ((((uint64_t)op[4] * 7947569409371862095UL) + ((uint64_t)op[5] * 6566749008003890043UL) + ((uint64_t)op[6] * 16913353855706661744UL) + ((uint64_t)op[7] * 11045776113192476354UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 11045776113192476354UL) + ((uint64_t)op[1] * 7008775382811200779UL) + ((uint64_t)op[2] * 4854050199574680321UL) + ((uint64_t)op[3] * 9968235917052480426UL) + ((uint64_t)op[4] * 6747879995015643189UL) + ((((uint64_t)op[5] * 7947569409371862095UL) + ((uint64_t)op[6] * 6566749008003890043UL) + ((uint64_t)op[7] * 16913353855706661744UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 16913353855706661744UL) + ((uint64_t)op[1] * 11045776113192476354UL) + ((uint64_t)op[2] * 7008775382811200779UL) + ((uint64_t)op[3] * 4854050199574680321UL) + ((uint64_t)op[4] * 9968235917052480426UL) + ((uint64_t)op[5] * 6747879995015643189UL) + ((((uint64_t)op[6] * 7947569409371862095UL) + ((uint64_t)op[7] * 6566749008003890043UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 6566749008003890043UL) + ((uint64_t)op[1] * 16913353855706661744UL) + ((uint64_t)op[2] * 11045776113192476354UL) + ((uint64_t)op[3] * 7008775382811200779UL) + ((uint64_t)op[4] * 4854050199574680321UL) + ((uint64_t)op[5] * 9968235917052480426UL) + ((uint64_t)op[6] * 6747879995015643189UL) + ((uint64_t)op[7] * 13050779919303516947UL);
	tmp_q[7] = ((uint64_t)op[0] * 7947569409371862095UL) + ((uint64_t)op[1] * 6566749008003890043UL) + ((uint64_t)op[2] * 16913353855706661744UL) + ((uint64_t)op[3] * 11045776113192476354UL) + ((uint64_t)op[4] * 7008775382811200779UL) + ((uint64_t)op[5] * 4854050199574680321UL) + ((uint64_t)op[6] * 9968235917052480426UL) + ((uint64_t)op[7] * 6747879995015643189UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 173334457162913L) - ((-((int128)tmp_q[1] * 139421154739237L) - ((int128)tmp_q[2] * 81215401766359L) - ((int128)tmp_q[3] * 33084910600736L) + ((int128)tmp_q[4] * 42171771305404L) - ((int128)tmp_q[5] * 130893386118697L) + ((int128)tmp_q[6] * 113703450776263L) + ((int128)tmp_q[7] * 130348825225994L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 130348825225994L) + ((int128)tmp_q[1] * 173334457162913L) - ((-((int128)tmp_q[2] * 139421154739237L) - ((int128)tmp_q[3] * 81215401766359L) - ((int128)tmp_q[4] * 33084910600736L) + ((int128)tmp_q[5] * 42171771305404L) - ((int128)tmp_q[6] * 130893386118697L) + ((int128)tmp_q[7] * 113703450776263L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 113703450776263L) + ((int128)tmp_q[1] * 130348825225994L) + ((int128)tmp_q[2] * 173334457162913L) - ((-((int128)tmp_q[3] * 139421154739237L) - ((int128)tmp_q[4] * 81215401766359L) - ((int128)tmp_q[5] * 33084910600736L) + ((int128)tmp_q[6] * 42171771305404L) - ((int128)tmp_q[7] * 130893386118697L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 130893386118697L) + ((int128)tmp_q[1] * 113703450776263L) + ((int128)tmp_q[2] * 130348825225994L) + ((int128)tmp_q[3] * 173334457162913L) - ((-((int128)tmp_q[4] * 139421154739237L) - ((int128)tmp_q[5] * 81215401766359L) - ((int128)tmp_q[6] * 33084910600736L) + ((int128)tmp_q[7] * 42171771305404L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 42171771305404L) - ((int128)tmp_q[1] * 130893386118697L) + ((int128)tmp_q[2] * 113703450776263L) + ((int128)tmp_q[3] * 130348825225994L) + ((int128)tmp_q[4] * 173334457162913L) - ((-((int128)tmp_q[5] * 139421154739237L) - ((int128)tmp_q[6] * 81215401766359L) - ((int128)tmp_q[7] * 33084910600736L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 33084910600736L) + ((int128)tmp_q[1] * 42171771305404L) - ((int128)tmp_q[2] * 130893386118697L) + ((int128)tmp_q[3] * 113703450776263L) + ((int128)tmp_q[4] * 130348825225994L) + ((int128)tmp_q[5] * 173334457162913L) - ((-((int128)tmp_q[6] * 139421154739237L) - ((int128)tmp_q[7] * 81215401766359L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 81215401766359L) - ((int128)tmp_q[1] * 33084910600736L) + ((int128)tmp_q[2] * 42171771305404L) - ((int128)tmp_q[3] * 130893386118697L) + ((int128)tmp_q[4] * 113703450776263L) + ((int128)tmp_q[5] * 130348825225994L) + ((int128)tmp_q[6] * 173334457162913L) + ((int128)tmp_q[7] * 418263464217711L);
	tmp_zero[7] = -((int128)tmp_q[0] * 139421154739237L) - ((int128)tmp_q[1] * 81215401766359L) - ((int128)tmp_q[2] * 33084910600736L) + ((int128)tmp_q[3] * 42171771305404L) - ((int128)tmp_q[4] * 130893386118697L) + ((int128)tmp_q[5] * 113703450776263L) + ((int128)tmp_q[6] * 130348825225994L) + ((int128)tmp_q[7] * 173334457162913L);

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

