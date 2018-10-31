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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5319433351767182970UL) + ((((uint64_t)op[1] * 200795756428594956UL) + ((uint64_t)op[2] * 5340656938932321901UL) + ((uint64_t)op[3] * 2108088056516584927UL) + ((uint64_t)op[4] * 16923978709828004525UL) + ((uint64_t)op[5] * 6153899087331820559UL) + ((uint64_t)op[6] * 2551696665027611700UL) + ((uint64_t)op[7] * 2076454713243147415UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 2076454713243147415UL) + ((uint64_t)op[1] * 5319433351767182970UL) + ((((uint64_t)op[2] * 200795756428594956UL) + ((uint64_t)op[3] * 5340656938932321901UL) + ((uint64_t)op[4] * 2108088056516584927UL) + ((uint64_t)op[5] * 16923978709828004525UL) + ((uint64_t)op[6] * 6153899087331820559UL) + ((uint64_t)op[7] * 2551696665027611700UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 2551696665027611700UL) + ((uint64_t)op[1] * 2076454713243147415UL) + ((uint64_t)op[2] * 5319433351767182970UL) + ((((uint64_t)op[3] * 200795756428594956UL) + ((uint64_t)op[4] * 5340656938932321901UL) + ((uint64_t)op[5] * 2108088056516584927UL) + ((uint64_t)op[6] * 16923978709828004525UL) + ((uint64_t)op[7] * 6153899087331820559UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 6153899087331820559UL) + ((uint64_t)op[1] * 2551696665027611700UL) + ((uint64_t)op[2] * 2076454713243147415UL) + ((uint64_t)op[3] * 5319433351767182970UL) + ((((uint64_t)op[4] * 200795756428594956UL) + ((uint64_t)op[5] * 5340656938932321901UL) + ((uint64_t)op[6] * 2108088056516584927UL) + ((uint64_t)op[7] * 16923978709828004525UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 16923978709828004525UL) + ((uint64_t)op[1] * 6153899087331820559UL) + ((uint64_t)op[2] * 2551696665027611700UL) + ((uint64_t)op[3] * 2076454713243147415UL) + ((uint64_t)op[4] * 5319433351767182970UL) + ((((uint64_t)op[5] * 200795756428594956UL) + ((uint64_t)op[6] * 5340656938932321901UL) + ((uint64_t)op[7] * 2108088056516584927UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 2108088056516584927UL) + ((uint64_t)op[1] * 16923978709828004525UL) + ((uint64_t)op[2] * 6153899087331820559UL) + ((uint64_t)op[3] * 2551696665027611700UL) + ((uint64_t)op[4] * 2076454713243147415UL) + ((uint64_t)op[5] * 5319433351767182970UL) + ((((uint64_t)op[6] * 200795756428594956UL) + ((uint64_t)op[7] * 5340656938932321901UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 5340656938932321901UL) + ((uint64_t)op[1] * 2108088056516584927UL) + ((uint64_t)op[2] * 16923978709828004525UL) + ((uint64_t)op[3] * 6153899087331820559UL) + ((uint64_t)op[4] * 2551696665027611700UL) + ((uint64_t)op[5] * 2076454713243147415UL) + ((uint64_t)op[6] * 5319433351767182970UL) + ((uint64_t)op[7] * 1405570295000164692UL);
	tmp_q[7] = ((uint64_t)op[0] * 200795756428594956UL) + ((uint64_t)op[1] * 5340656938932321901UL) + ((uint64_t)op[2] * 2108088056516584927UL) + ((uint64_t)op[3] * 16923978709828004525UL) + ((uint64_t)op[4] * 6153899087331820559UL) + ((uint64_t)op[5] * 2551696665027611700UL) + ((uint64_t)op[6] * 2076454713243147415UL) + ((uint64_t)op[7] * 5319433351767182970UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 133611329235012L) + ((-((int128)tmp_q[1] * 21667317131276L) - ((int128)tmp_q[2] * 57717660479476L) + ((int128)tmp_q[3] * 76587683174003L) + ((int128)tmp_q[4] * 127231524516633L) - ((int128)tmp_q[5] * 147337720652194L) + ((int128)tmp_q[6] * 51453505400869L) + ((int128)tmp_q[7] * 49039148372544L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 49039148372544L) - ((int128)tmp_q[1] * 133611329235012L) + ((-((int128)tmp_q[2] * 21667317131276L) - ((int128)tmp_q[3] * 57717660479476L) + ((int128)tmp_q[4] * 76587683174003L) + ((int128)tmp_q[5] * 127231524516633L) - ((int128)tmp_q[6] * 147337720652194L) + ((int128)tmp_q[7] * 51453505400869L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 51453505400869L) + ((int128)tmp_q[1] * 49039148372544L) - ((int128)tmp_q[2] * 133611329235012L) + ((-((int128)tmp_q[3] * 21667317131276L) - ((int128)tmp_q[4] * 57717660479476L) + ((int128)tmp_q[5] * 76587683174003L) + ((int128)tmp_q[6] * 127231524516633L) - ((int128)tmp_q[7] * 147337720652194L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 147337720652194L) + ((int128)tmp_q[1] * 51453505400869L) + ((int128)tmp_q[2] * 49039148372544L) - ((int128)tmp_q[3] * 133611329235012L) + ((-((int128)tmp_q[4] * 21667317131276L) - ((int128)tmp_q[5] * 57717660479476L) + ((int128)tmp_q[6] * 76587683174003L) + ((int128)tmp_q[7] * 127231524516633L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 127231524516633L) - ((int128)tmp_q[1] * 147337720652194L) + ((int128)tmp_q[2] * 51453505400869L) + ((int128)tmp_q[3] * 49039148372544L) - ((int128)tmp_q[4] * 133611329235012L) + ((-((int128)tmp_q[5] * 21667317131276L) - ((int128)tmp_q[6] * 57717660479476L) + ((int128)tmp_q[7] * 76587683174003L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 76587683174003L) + ((int128)tmp_q[1] * 127231524516633L) - ((int128)tmp_q[2] * 147337720652194L) + ((int128)tmp_q[3] * 51453505400869L) + ((int128)tmp_q[4] * 49039148372544L) - ((int128)tmp_q[5] * 133611329235012L) + ((-((int128)tmp_q[6] * 21667317131276L) - ((int128)tmp_q[7] * 57717660479476L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 57717660479476L) + ((int128)tmp_q[1] * 76587683174003L) + ((int128)tmp_q[2] * 127231524516633L) - ((int128)tmp_q[3] * 147337720652194L) + ((int128)tmp_q[4] * 51453505400869L) + ((int128)tmp_q[5] * 49039148372544L) - ((int128)tmp_q[6] * 133611329235012L) - ((int128)tmp_q[7] * 151671219918932L);
	tmp_zero[7] = -((int128)tmp_q[0] * 21667317131276L) - ((int128)tmp_q[1] * 57717660479476L) + ((int128)tmp_q[2] * 76587683174003L) + ((int128)tmp_q[3] * 127231524516633L) - ((int128)tmp_q[4] * 147337720652194L) + ((int128)tmp_q[5] * 51453505400869L) + ((int128)tmp_q[6] * 49039148372544L) - ((int128)tmp_q[7] * 133611329235012L);

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

