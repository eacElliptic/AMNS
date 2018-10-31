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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6724408170832475139UL) + ((((uint64_t)op[1] * 8552135142536821830UL) + ((uint64_t)op[2] * 3585266031550457941UL) + ((uint64_t)op[3] * 97665636649665498UL) + ((uint64_t)op[4] * 149081223204560213UL) + ((uint64_t)op[5] * 12158269255358731044UL) + ((uint64_t)op[6] * 13331117755444085713UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 13331117755444085713UL) + ((uint64_t)op[1] * 6724408170832475139UL) + ((((uint64_t)op[2] * 8552135142536821830UL) + ((uint64_t)op[3] * 3585266031550457941UL) + ((uint64_t)op[4] * 97665636649665498UL) + ((uint64_t)op[5] * 149081223204560213UL) + ((uint64_t)op[6] * 12158269255358731044UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 12158269255358731044UL) + ((uint64_t)op[1] * 13331117755444085713UL) + ((uint64_t)op[2] * 6724408170832475139UL) + ((((uint64_t)op[3] * 8552135142536821830UL) + ((uint64_t)op[4] * 3585266031550457941UL) + ((uint64_t)op[5] * 97665636649665498UL) + ((uint64_t)op[6] * 149081223204560213UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 149081223204560213UL) + ((uint64_t)op[1] * 12158269255358731044UL) + ((uint64_t)op[2] * 13331117755444085713UL) + ((uint64_t)op[3] * 6724408170832475139UL) + ((((uint64_t)op[4] * 8552135142536821830UL) + ((uint64_t)op[5] * 3585266031550457941UL) + ((uint64_t)op[6] * 97665636649665498UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 97665636649665498UL) + ((uint64_t)op[1] * 149081223204560213UL) + ((uint64_t)op[2] * 12158269255358731044UL) + ((uint64_t)op[3] * 13331117755444085713UL) + ((uint64_t)op[4] * 6724408170832475139UL) + ((((uint64_t)op[5] * 8552135142536821830UL) + ((uint64_t)op[6] * 3585266031550457941UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 3585266031550457941UL) + ((uint64_t)op[1] * 97665636649665498UL) + ((uint64_t)op[2] * 149081223204560213UL) + ((uint64_t)op[3] * 12158269255358731044UL) + ((uint64_t)op[4] * 13331117755444085713UL) + ((uint64_t)op[5] * 6724408170832475139UL) + ((uint64_t)op[6] * 17104270285073643660UL);
	tmp_q[6] = ((uint64_t)op[0] * 8552135142536821830UL) + ((uint64_t)op[1] * 3585266031550457941UL) + ((uint64_t)op[2] * 97665636649665498UL) + ((uint64_t)op[3] * 149081223204560213UL) + ((uint64_t)op[4] * 12158269255358731044UL) + ((uint64_t)op[5] * 13331117755444085713UL) + ((uint64_t)op[6] * 6724408170832475139UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 10445419062089851L) + ((-((int128)tmp_q[1] * 14581186199786352L) - ((int128)tmp_q[2] * 2641792348377871L) - ((int128)tmp_q[3] * 4073950979238679L) + ((int128)tmp_q[4] * 11434695089154052L) - ((int128)tmp_q[5] * 6686520524857831L) + ((int128)tmp_q[6] * 3635813147981459L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 3635813147981459L) - ((int128)tmp_q[1] * 10445419062089851L) + ((-((int128)tmp_q[2] * 14581186199786352L) - ((int128)tmp_q[3] * 2641792348377871L) - ((int128)tmp_q[4] * 4073950979238679L) + ((int128)tmp_q[5] * 11434695089154052L) - ((int128)tmp_q[6] * 6686520524857831L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 6686520524857831L) + ((int128)tmp_q[1] * 3635813147981459L) - ((int128)tmp_q[2] * 10445419062089851L) + ((-((int128)tmp_q[3] * 14581186199786352L) - ((int128)tmp_q[4] * 2641792348377871L) - ((int128)tmp_q[5] * 4073950979238679L) + ((int128)tmp_q[6] * 11434695089154052L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 11434695089154052L) - ((int128)tmp_q[1] * 6686520524857831L) + ((int128)tmp_q[2] * 3635813147981459L) - ((int128)tmp_q[3] * 10445419062089851L) + ((-((int128)tmp_q[4] * 14581186199786352L) - ((int128)tmp_q[5] * 2641792348377871L) - ((int128)tmp_q[6] * 4073950979238679L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 4073950979238679L) + ((int128)tmp_q[1] * 11434695089154052L) - ((int128)tmp_q[2] * 6686520524857831L) + ((int128)tmp_q[3] * 3635813147981459L) - ((int128)tmp_q[4] * 10445419062089851L) + ((-((int128)tmp_q[5] * 14581186199786352L) - ((int128)tmp_q[6] * 2641792348377871L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 2641792348377871L) - ((int128)tmp_q[1] * 4073950979238679L) + ((int128)tmp_q[2] * 11434695089154052L) - ((int128)tmp_q[3] * 6686520524857831L) + ((int128)tmp_q[4] * 3635813147981459L) - ((int128)tmp_q[5] * 10445419062089851L) - ((int128)tmp_q[6] * 29162372399572704L);
	tmp_zero[6] = -((int128)tmp_q[0] * 14581186199786352L) - ((int128)tmp_q[1] * 2641792348377871L) - ((int128)tmp_q[2] * 4073950979238679L) + ((int128)tmp_q[3] * 11434695089154052L) - ((int128)tmp_q[4] * 6686520524857831L) + ((int128)tmp_q[5] * 3635813147981459L) - ((int128)tmp_q[6] * 10445419062089851L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

