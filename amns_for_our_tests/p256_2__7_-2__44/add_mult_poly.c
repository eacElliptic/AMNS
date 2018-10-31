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
	tmp_q[0] = ((uint64_t)op[0] * 15042178119535375441UL) + ((((uint64_t)op[1] * 18006788646168884603UL) + ((uint64_t)op[2] * 8195403165923820132UL) + ((uint64_t)op[3] * 4832709365872890339UL) + ((uint64_t)op[4] * 734963362961594652UL) + ((uint64_t)op[5] * 11103381719074146542UL) + ((uint64_t)op[6] * 16505576772107286211UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 16505576772107286211UL) + ((uint64_t)op[1] * 15042178119535375441UL) + ((((uint64_t)op[2] * 18006788646168884603UL) + ((uint64_t)op[3] * 8195403165923820132UL) + ((uint64_t)op[4] * 4832709365872890339UL) + ((uint64_t)op[5] * 734963362961594652UL) + ((uint64_t)op[6] * 11103381719074146542UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 11103381719074146542UL) + ((uint64_t)op[1] * 16505576772107286211UL) + ((uint64_t)op[2] * 15042178119535375441UL) + ((((uint64_t)op[3] * 18006788646168884603UL) + ((uint64_t)op[4] * 8195403165923820132UL) + ((uint64_t)op[5] * 4832709365872890339UL) + ((uint64_t)op[6] * 734963362961594652UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 734963362961594652UL) + ((uint64_t)op[1] * 11103381719074146542UL) + ((uint64_t)op[2] * 16505576772107286211UL) + ((uint64_t)op[3] * 15042178119535375441UL) + ((((uint64_t)op[4] * 18006788646168884603UL) + ((uint64_t)op[5] * 8195403165923820132UL) + ((uint64_t)op[6] * 4832709365872890339UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 4832709365872890339UL) + ((uint64_t)op[1] * 734963362961594652UL) + ((uint64_t)op[2] * 11103381719074146542UL) + ((uint64_t)op[3] * 16505576772107286211UL) + ((uint64_t)op[4] * 15042178119535375441UL) + ((((uint64_t)op[5] * 18006788646168884603UL) + ((uint64_t)op[6] * 8195403165923820132UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 8195403165923820132UL) + ((uint64_t)op[1] * 4832709365872890339UL) + ((uint64_t)op[2] * 734963362961594652UL) + ((uint64_t)op[3] * 11103381719074146542UL) + ((uint64_t)op[4] * 16505576772107286211UL) + ((uint64_t)op[5] * 15042178119535375441UL) + ((uint64_t)op[6] * 879910855081334026UL);
	tmp_q[6] = ((uint64_t)op[0] * 18006788646168884603UL) + ((uint64_t)op[1] * 8195403165923820132UL) + ((uint64_t)op[2] * 4832709365872890339UL) + ((uint64_t)op[3] * 734963362961594652UL) + ((uint64_t)op[4] * 11103381719074146542UL) + ((uint64_t)op[5] * 16505576772107286211UL) + ((uint64_t)op[6] * 15042178119535375441UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 45130217155L) - ((-((int128)tmp_q[1] * 32808727853L) + ((int128)tmp_q[2] * 6252050973L) + ((int128)tmp_q[3] * 3382303182L) - ((int128)tmp_q[4] * 51048307547L) - ((int128)tmp_q[5] * 51584921287L) - ((int128)tmp_q[6] * 38628067833L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 38628067833L) - ((int128)tmp_q[1] * 45130217155L) - ((-((int128)tmp_q[2] * 32808727853L) + ((int128)tmp_q[3] * 6252050973L) + ((int128)tmp_q[4] * 3382303182L) - ((int128)tmp_q[5] * 51048307547L) - ((int128)tmp_q[6] * 51584921287L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 51584921287L) - ((int128)tmp_q[1] * 38628067833L) - ((int128)tmp_q[2] * 45130217155L) - ((-((int128)tmp_q[3] * 32808727853L) + ((int128)tmp_q[4] * 6252050973L) + ((int128)tmp_q[5] * 3382303182L) - ((int128)tmp_q[6] * 51048307547L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 51048307547L) - ((int128)tmp_q[1] * 51584921287L) - ((int128)tmp_q[2] * 38628067833L) - ((int128)tmp_q[3] * 45130217155L) - ((-((int128)tmp_q[4] * 32808727853L) + ((int128)tmp_q[5] * 6252050973L) + ((int128)tmp_q[6] * 3382303182L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 3382303182L) - ((int128)tmp_q[1] * 51048307547L) - ((int128)tmp_q[2] * 51584921287L) - ((int128)tmp_q[3] * 38628067833L) - ((int128)tmp_q[4] * 45130217155L) - ((-((int128)tmp_q[5] * 32808727853L) + ((int128)tmp_q[6] * 6252050973L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 6252050973L) + ((int128)tmp_q[1] * 3382303182L) - ((int128)tmp_q[2] * 51048307547L) - ((int128)tmp_q[3] * 51584921287L) - ((int128)tmp_q[4] * 38628067833L) - ((int128)tmp_q[5] * 45130217155L) + ((int128)tmp_q[6] * 65617455706L);
	tmp_zero[6] = -((int128)tmp_q[0] * 32808727853L) + ((int128)tmp_q[1] * 6252050973L) + ((int128)tmp_q[2] * 3382303182L) - ((int128)tmp_q[3] * 51048307547L) - ((int128)tmp_q[4] * 51584921287L) - ((int128)tmp_q[5] * 38628067833L) - ((int128)tmp_q[6] * 45130217155L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

