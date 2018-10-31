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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8513081643842350781UL) + ((((uint64_t)op[1] * 5372451320155786022UL) + ((uint64_t)op[2] * 12432549350869629074UL) + ((uint64_t)op[3] * 4857259891832412275UL) + ((uint64_t)op[4] * 9820946361970800903UL) + ((uint64_t)op[5] * 1751316598610767100UL) + ((uint64_t)op[6] * 1814839974827018052UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 1814839974827018052UL) + ((uint64_t)op[1] * 8513081643842350781UL) + ((((uint64_t)op[2] * 5372451320155786022UL) + ((uint64_t)op[3] * 12432549350869629074UL) + ((uint64_t)op[4] * 4857259891832412275UL) + ((uint64_t)op[5] * 9820946361970800903UL) + ((uint64_t)op[6] * 1751316598610767100UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 1751316598610767100UL) + ((uint64_t)op[1] * 1814839974827018052UL) + ((uint64_t)op[2] * 8513081643842350781UL) + ((((uint64_t)op[3] * 5372451320155786022UL) + ((uint64_t)op[4] * 12432549350869629074UL) + ((uint64_t)op[5] * 4857259891832412275UL) + ((uint64_t)op[6] * 9820946361970800903UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 9820946361970800903UL) + ((uint64_t)op[1] * 1751316598610767100UL) + ((uint64_t)op[2] * 1814839974827018052UL) + ((uint64_t)op[3] * 8513081643842350781UL) + ((((uint64_t)op[4] * 5372451320155786022UL) + ((uint64_t)op[5] * 12432549350869629074UL) + ((uint64_t)op[6] * 4857259891832412275UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 4857259891832412275UL) + ((uint64_t)op[1] * 9820946361970800903UL) + ((uint64_t)op[2] * 1751316598610767100UL) + ((uint64_t)op[3] * 1814839974827018052UL) + ((uint64_t)op[4] * 8513081643842350781UL) + ((((uint64_t)op[5] * 5372451320155786022UL) + ((uint64_t)op[6] * 12432549350869629074UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 12432549350869629074UL) + ((uint64_t)op[1] * 4857259891832412275UL) + ((uint64_t)op[2] * 9820946361970800903UL) + ((uint64_t)op[3] * 1751316598610767100UL) + ((uint64_t)op[4] * 1814839974827018052UL) + ((uint64_t)op[5] * 8513081643842350781UL) + ((uint64_t)op[6] * 2329390113242193550UL);
	tmp_q[6] = ((uint64_t)op[0] * 5372451320155786022UL) + ((uint64_t)op[1] * 12432549350869629074UL) + ((uint64_t)op[2] * 4857259891832412275UL) + ((uint64_t)op[3] * 9820946361970800903UL) + ((uint64_t)op[4] * 1751316598610767100UL) + ((uint64_t)op[5] * 1814839974827018052UL) + ((uint64_t)op[6] * 8513081643842350781UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 30886969145L) - ((((int128)tmp_q[1] * 27527110376L) - ((int128)tmp_q[2] * 9028221277L) - ((int128)tmp_q[3] * 44537159641L) + ((int128)tmp_q[4] * 2900242031L) - ((int128)tmp_q[5] * 117258710789L) + ((int128)tmp_q[6] * 9436796810L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 9436796810L) - ((int128)tmp_q[1] * 30886969145L) - ((((int128)tmp_q[2] * 27527110376L) - ((int128)tmp_q[3] * 9028221277L) - ((int128)tmp_q[4] * 44537159641L) + ((int128)tmp_q[5] * 2900242031L) - ((int128)tmp_q[6] * 117258710789L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 117258710789L) + ((int128)tmp_q[1] * 9436796810L) - ((int128)tmp_q[2] * 30886969145L) - ((((int128)tmp_q[3] * 27527110376L) - ((int128)tmp_q[4] * 9028221277L) - ((int128)tmp_q[5] * 44537159641L) + ((int128)tmp_q[6] * 2900242031L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 2900242031L) - ((int128)tmp_q[1] * 117258710789L) + ((int128)tmp_q[2] * 9436796810L) - ((int128)tmp_q[3] * 30886969145L) - ((((int128)tmp_q[4] * 27527110376L) - ((int128)tmp_q[5] * 9028221277L) - ((int128)tmp_q[6] * 44537159641L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 44537159641L) + ((int128)tmp_q[1] * 2900242031L) - ((int128)tmp_q[2] * 117258710789L) + ((int128)tmp_q[3] * 9436796810L) - ((int128)tmp_q[4] * 30886969145L) - ((((int128)tmp_q[5] * 27527110376L) - ((int128)tmp_q[6] * 9028221277L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 9028221277L) - ((int128)tmp_q[1] * 44537159641L) + ((int128)tmp_q[2] * 2900242031L) - ((int128)tmp_q[3] * 117258710789L) + ((int128)tmp_q[4] * 9436796810L) - ((int128)tmp_q[5] * 30886969145L) - ((int128)tmp_q[6] * 82581331128L);
	tmp_zero[6] = ((int128)tmp_q[0] * 27527110376L) - ((int128)tmp_q[1] * 9028221277L) - ((int128)tmp_q[2] * 44537159641L) + ((int128)tmp_q[3] * 2900242031L) - ((int128)tmp_q[4] * 117258710789L) + ((int128)tmp_q[5] * 9436796810L) - ((int128)tmp_q[6] * 30886969145L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

