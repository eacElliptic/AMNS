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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17896400840486198815UL) + ((((uint64_t)op[1] * 12170599126770835310UL) + ((uint64_t)op[2] * 16268455540795920680UL) + ((uint64_t)op[3] * 9165057691112852048UL) + ((uint64_t)op[4] * 5793578204039702154UL) + ((uint64_t)op[5] * 9237580225959537643UL) + ((uint64_t)op[6] * 1158749317795765080UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 1158749317795765080UL) + ((uint64_t)op[1] * 17896400840486198815UL) + ((((uint64_t)op[2] * 12170599126770835310UL) + ((uint64_t)op[3] * 16268455540795920680UL) + ((uint64_t)op[4] * 9165057691112852048UL) + ((uint64_t)op[5] * 5793578204039702154UL) + ((uint64_t)op[6] * 9237580225959537643UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 9237580225959537643UL) + ((uint64_t)op[1] * 1158749317795765080UL) + ((uint64_t)op[2] * 17896400840486198815UL) + ((((uint64_t)op[3] * 12170599126770835310UL) + ((uint64_t)op[4] * 16268455540795920680UL) + ((uint64_t)op[5] * 9165057691112852048UL) + ((uint64_t)op[6] * 5793578204039702154UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 5793578204039702154UL) + ((uint64_t)op[1] * 9237580225959537643UL) + ((uint64_t)op[2] * 1158749317795765080UL) + ((uint64_t)op[3] * 17896400840486198815UL) + ((((uint64_t)op[4] * 12170599126770835310UL) + ((uint64_t)op[5] * 16268455540795920680UL) + ((uint64_t)op[6] * 9165057691112852048UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 9165057691112852048UL) + ((uint64_t)op[1] * 5793578204039702154UL) + ((uint64_t)op[2] * 9237580225959537643UL) + ((uint64_t)op[3] * 1158749317795765080UL) + ((uint64_t)op[4] * 17896400840486198815UL) + ((((uint64_t)op[5] * 12170599126770835310UL) + ((uint64_t)op[6] * 16268455540795920680UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 16268455540795920680UL) + ((uint64_t)op[1] * 9165057691112852048UL) + ((uint64_t)op[2] * 5793578204039702154UL) + ((uint64_t)op[3] * 9237580225959537643UL) + ((uint64_t)op[4] * 1158749317795765080UL) + ((uint64_t)op[5] * 17896400840486198815UL) + ((uint64_t)op[6] * 5131072645618924400UL);
	tmp_q[6] = ((uint64_t)op[0] * 12170599126770835310UL) + ((uint64_t)op[1] * 16268455540795920680UL) + ((uint64_t)op[2] * 9165057691112852048UL) + ((uint64_t)op[3] * 5793578204039702154UL) + ((uint64_t)op[4] * 9237580225959537643UL) + ((uint64_t)op[5] * 1158749317795765080UL) + ((uint64_t)op[6] * 17896400840486198815UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 10331610705L) + ((-((int128)tmp_q[1] * 20211387579L) + ((int128)tmp_q[2] * 33176866932L) - ((int128)tmp_q[3] * 16267621735L) + ((int128)tmp_q[4] * 44188503042L) + ((int128)tmp_q[5] * 61322762091L) + ((int128)tmp_q[6] * 568312480L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 568312480L) + ((int128)tmp_q[1] * 10331610705L) + ((-((int128)tmp_q[2] * 20211387579L) + ((int128)tmp_q[3] * 33176866932L) - ((int128)tmp_q[4] * 16267621735L) + ((int128)tmp_q[5] * 44188503042L) + ((int128)tmp_q[6] * 61322762091L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 61322762091L) + ((int128)tmp_q[1] * 568312480L) + ((int128)tmp_q[2] * 10331610705L) + ((-((int128)tmp_q[3] * 20211387579L) + ((int128)tmp_q[4] * 33176866932L) - ((int128)tmp_q[5] * 16267621735L) + ((int128)tmp_q[6] * 44188503042L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 44188503042L) + ((int128)tmp_q[1] * 61322762091L) + ((int128)tmp_q[2] * 568312480L) + ((int128)tmp_q[3] * 10331610705L) + ((-((int128)tmp_q[4] * 20211387579L) + ((int128)tmp_q[5] * 33176866932L) - ((int128)tmp_q[6] * 16267621735L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 16267621735L) + ((int128)tmp_q[1] * 44188503042L) + ((int128)tmp_q[2] * 61322762091L) + ((int128)tmp_q[3] * 568312480L) + ((int128)tmp_q[4] * 10331610705L) + ((-((int128)tmp_q[5] * 20211387579L) + ((int128)tmp_q[6] * 33176866932L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 33176866932L) - ((int128)tmp_q[1] * 16267621735L) + ((int128)tmp_q[2] * 44188503042L) + ((int128)tmp_q[3] * 61322762091L) + ((int128)tmp_q[4] * 568312480L) + ((int128)tmp_q[5] * 10331610705L) - ((int128)tmp_q[6] * 161691100632L);
	tmp_zero[6] = -((int128)tmp_q[0] * 20211387579L) + ((int128)tmp_q[1] * 33176866932L) - ((int128)tmp_q[2] * 16267621735L) + ((int128)tmp_q[3] * 44188503042L) + ((int128)tmp_q[4] * 61322762091L) + ((int128)tmp_q[5] * 568312480L) + ((int128)tmp_q[6] * 10331610705L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

