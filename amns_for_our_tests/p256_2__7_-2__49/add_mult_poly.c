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
	tmp_q[0] = ((uint64_t)op[0] * 7764172757378251353UL) + ((((uint64_t)op[1] * 9856548506353895979UL) + ((uint64_t)op[2] * 2987743363585732343UL) + ((uint64_t)op[3] * 1221047436249972043UL) + ((uint64_t)op[4] * 14369581604987370562UL) + ((uint64_t)op[5] * 2457836853558069972UL) + ((uint64_t)op[6] * 15758938278413300233UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 15758938278413300233UL) + ((uint64_t)op[1] * 7764172757378251353UL) + ((((uint64_t)op[2] * 9856548506353895979UL) + ((uint64_t)op[3] * 2987743363585732343UL) + ((uint64_t)op[4] * 1221047436249972043UL) + ((uint64_t)op[5] * 14369581604987370562UL) + ((uint64_t)op[6] * 2457836853558069972UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2457836853558069972UL) + ((uint64_t)op[1] * 15758938278413300233UL) + ((uint64_t)op[2] * 7764172757378251353UL) + ((((uint64_t)op[3] * 9856548506353895979UL) + ((uint64_t)op[4] * 2987743363585732343UL) + ((uint64_t)op[5] * 1221047436249972043UL) + ((uint64_t)op[6] * 14369581604987370562UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 14369581604987370562UL) + ((uint64_t)op[1] * 2457836853558069972UL) + ((uint64_t)op[2] * 15758938278413300233UL) + ((uint64_t)op[3] * 7764172757378251353UL) + ((((uint64_t)op[4] * 9856548506353895979UL) + ((uint64_t)op[5] * 2987743363585732343UL) + ((uint64_t)op[6] * 1221047436249972043UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 1221047436249972043UL) + ((uint64_t)op[1] * 14369581604987370562UL) + ((uint64_t)op[2] * 2457836853558069972UL) + ((uint64_t)op[3] * 15758938278413300233UL) + ((uint64_t)op[4] * 7764172757378251353UL) + ((((uint64_t)op[5] * 9856548506353895979UL) + ((uint64_t)op[6] * 2987743363585732343UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 2987743363585732343UL) + ((uint64_t)op[1] * 1221047436249972043UL) + ((uint64_t)op[2] * 14369581604987370562UL) + ((uint64_t)op[3] * 2457836853558069972UL) + ((uint64_t)op[4] * 15758938278413300233UL) + ((uint64_t)op[5] * 7764172757378251353UL) + ((uint64_t)op[6] * 17180391134711311274UL);
	tmp_q[6] = ((uint64_t)op[0] * 9856548506353895979UL) + ((uint64_t)op[1] * 2987743363585732343UL) + ((uint64_t)op[2] * 1221047436249972043UL) + ((uint64_t)op[3] * 14369581604987370562UL) + ((uint64_t)op[4] * 2457836853558069972UL) + ((uint64_t)op[5] * 15758938278413300233UL) + ((uint64_t)op[6] * 7764172757378251353UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 6767727611L) - ((((int128)tmp_q[1] * 46514103959L) - ((int128)tmp_q[2] * 51309881248L) + ((int128)tmp_q[3] * 12688644400L) - ((int128)tmp_q[4] * 54670880221L) - ((int128)tmp_q[5] * 33490990291L) - ((int128)tmp_q[6] * 21880581803L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 21880581803L) + ((int128)tmp_q[1] * 6767727611L) - ((((int128)tmp_q[2] * 46514103959L) - ((int128)tmp_q[3] * 51309881248L) + ((int128)tmp_q[4] * 12688644400L) - ((int128)tmp_q[5] * 54670880221L) - ((int128)tmp_q[6] * 33490990291L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 33490990291L) - ((int128)tmp_q[1] * 21880581803L) + ((int128)tmp_q[2] * 6767727611L) - ((((int128)tmp_q[3] * 46514103959L) - ((int128)tmp_q[4] * 51309881248L) + ((int128)tmp_q[5] * 12688644400L) - ((int128)tmp_q[6] * 54670880221L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 54670880221L) - ((int128)tmp_q[1] * 33490990291L) - ((int128)tmp_q[2] * 21880581803L) + ((int128)tmp_q[3] * 6767727611L) - ((((int128)tmp_q[4] * 46514103959L) - ((int128)tmp_q[5] * 51309881248L) + ((int128)tmp_q[6] * 12688644400L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 12688644400L) - ((int128)tmp_q[1] * 54670880221L) - ((int128)tmp_q[2] * 33490990291L) - ((int128)tmp_q[3] * 21880581803L) + ((int128)tmp_q[4] * 6767727611L) - ((((int128)tmp_q[5] * 46514103959L) - ((int128)tmp_q[6] * 51309881248L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 51309881248L) + ((int128)tmp_q[1] * 12688644400L) - ((int128)tmp_q[2] * 54670880221L) - ((int128)tmp_q[3] * 33490990291L) - ((int128)tmp_q[4] * 21880581803L) + ((int128)tmp_q[5] * 6767727611L) - ((int128)tmp_q[6] * 93028207918L);
	tmp_zero[6] = ((int128)tmp_q[0] * 46514103959L) - ((int128)tmp_q[1] * 51309881248L) + ((int128)tmp_q[2] * 12688644400L) - ((int128)tmp_q[3] * 54670880221L) - ((int128)tmp_q[4] * 33490990291L) - ((int128)tmp_q[5] * 21880581803L) + ((int128)tmp_q[6] * 6767727611L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

