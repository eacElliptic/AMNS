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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 511594267145265943UL) + ((((uint64_t)op[1] * 11155720387567473234UL) + ((uint64_t)op[2] * 1905751078203060705UL) + ((uint64_t)op[3] * 9916052860725889900UL) + ((uint64_t)op[4] * 17127985665295331503UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 17127985665295331503UL) + ((uint64_t)op[1] * 511594267145265943UL) + ((((uint64_t)op[2] * 11155720387567473234UL) + ((uint64_t)op[3] * 1905751078203060705UL) + ((uint64_t)op[4] * 9916052860725889900UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 9916052860725889900UL) + ((uint64_t)op[1] * 17127985665295331503UL) + ((uint64_t)op[2] * 511594267145265943UL) + ((((uint64_t)op[3] * 11155720387567473234UL) + ((uint64_t)op[4] * 1905751078203060705UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 1905751078203060705UL) + ((uint64_t)op[1] * 9916052860725889900UL) + ((uint64_t)op[2] * 17127985665295331503UL) + ((uint64_t)op[3] * 511594267145265943UL) + ((uint64_t)op[4] * 2987957268007972208UL);
	tmp_q[4] = ((uint64_t)op[0] * 11155720387567473234UL) + ((uint64_t)op[1] * 1905751078203060705UL) + ((uint64_t)op[2] * 9916052860725889900UL) + ((uint64_t)op[3] * 17127985665295331503UL) + ((uint64_t)op[4] * 511594267145265943UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 96292338031959L) - ((((int128)tmp_q[1] * 1202376957299621L) + ((int128)tmp_q[2] * 1148315896303824L) - ((int128)tmp_q[3] * 143348822507539L) + ((int128)tmp_q[4] * 369488178752959L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 369488178752959L) - ((int128)tmp_q[1] * 96292338031959L) - ((((int128)tmp_q[2] * 1202376957299621L) + ((int128)tmp_q[3] * 1148315896303824L) - ((int128)tmp_q[4] * 143348822507539L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 143348822507539L) + ((int128)tmp_q[1] * 369488178752959L) - ((int128)tmp_q[2] * 96292338031959L) - ((((int128)tmp_q[3] * 1202376957299621L) + ((int128)tmp_q[4] * 1148315896303824L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 1148315896303824L) - ((int128)tmp_q[1] * 143348822507539L) + ((int128)tmp_q[2] * 369488178752959L) - ((int128)tmp_q[3] * 96292338031959L) - ((int128)tmp_q[4] * 9619015658396968L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1202376957299621L) + ((int128)tmp_q[1] * 1148315896303824L) - ((int128)tmp_q[2] * 143348822507539L) + ((int128)tmp_q[3] * 369488178752959L) - ((int128)tmp_q[4] * 96292338031959L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

