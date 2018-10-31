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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 968428250787770017UL) + ((((uint64_t)op[1] * 6718859142641885792UL) + ((uint64_t)op[2] * 5436098786265108828UL) + ((uint64_t)op[3] * 16470005894427517179UL) + ((uint64_t)op[4] * 672214114642496590UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 672214114642496590UL) + ((uint64_t)op[1] * 968428250787770017UL) + ((((uint64_t)op[2] * 6718859142641885792UL) + ((uint64_t)op[3] * 5436098786265108828UL) + ((uint64_t)op[4] * 16470005894427517179UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 16470005894427517179UL) + ((uint64_t)op[1] * 672214114642496590UL) + ((uint64_t)op[2] * 968428250787770017UL) + ((((uint64_t)op[3] * 6718859142641885792UL) + ((uint64_t)op[4] * 5436098786265108828UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 5436098786265108828UL) + ((uint64_t)op[1] * 16470005894427517179UL) + ((uint64_t)op[2] * 672214114642496590UL) + ((uint64_t)op[3] * 968428250787770017UL) + ((uint64_t)op[4] * 3419666708432211520UL);
	tmp_q[4] = ((uint64_t)op[0] * 6718859142641885792UL) + ((uint64_t)op[1] * 5436098786265108828UL) + ((uint64_t)op[2] * 16470005894427517179UL) + ((uint64_t)op[3] * 672214114642496590UL) + ((uint64_t)op[4] * 968428250787770017UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 265262313736343L) + ((-((int128)tmp_q[1] * 906801969144533L) + ((int128)tmp_q[2] * 2009230641712378L) + ((int128)tmp_q[3] * 336442880530635L) - ((int128)tmp_q[4] * 985747630022600L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 985747630022600L) + ((int128)tmp_q[1] * 265262313736343L) + ((-((int128)tmp_q[2] * 906801969144533L) + ((int128)tmp_q[3] * 2009230641712378L) + ((int128)tmp_q[4] * 336442880530635L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 336442880530635L) - ((int128)tmp_q[1] * 985747630022600L) + ((int128)tmp_q[2] * 265262313736343L) + ((-((int128)tmp_q[3] * 906801969144533L) + ((int128)tmp_q[4] * 2009230641712378L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 2009230641712378L) + ((int128)tmp_q[1] * 336442880530635L) - ((int128)tmp_q[2] * 985747630022600L) + ((int128)tmp_q[3] * 265262313736343L) - ((int128)tmp_q[4] * 5440811814867198L);
	tmp_zero[4] = -((int128)tmp_q[0] * 906801969144533L) + ((int128)tmp_q[1] * 2009230641712378L) + ((int128)tmp_q[2] * 336442880530635L) - ((int128)tmp_q[3] * 985747630022600L) + ((int128)tmp_q[4] * 265262313736343L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

