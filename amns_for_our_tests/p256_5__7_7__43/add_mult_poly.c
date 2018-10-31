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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 436304521616961527UL) + ((((uint64_t)op[1] * 15854511070809174569UL) + ((uint64_t)op[2] * 4411907939086405018UL) + ((uint64_t)op[3] * 6824475643257457942UL) + ((uint64_t)op[4] * 2803642121071829719UL) + ((uint64_t)op[5] * 6478878482505650402UL) + ((uint64_t)op[6] * 3285739491163352768UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 3285739491163352768UL) + ((uint64_t)op[1] * 436304521616961527UL) + ((((uint64_t)op[2] * 15854511070809174569UL) + ((uint64_t)op[3] * 4411907939086405018UL) + ((uint64_t)op[4] * 6824475643257457942UL) + ((uint64_t)op[5] * 2803642121071829719UL) + ((uint64_t)op[6] * 6478878482505650402UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 6478878482505650402UL) + ((uint64_t)op[1] * 3285739491163352768UL) + ((uint64_t)op[2] * 436304521616961527UL) + ((((uint64_t)op[3] * 15854511070809174569UL) + ((uint64_t)op[4] * 4411907939086405018UL) + ((uint64_t)op[5] * 6824475643257457942UL) + ((uint64_t)op[6] * 2803642121071829719UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 2803642121071829719UL) + ((uint64_t)op[1] * 6478878482505650402UL) + ((uint64_t)op[2] * 3285739491163352768UL) + ((uint64_t)op[3] * 436304521616961527UL) + ((((uint64_t)op[4] * 15854511070809174569UL) + ((uint64_t)op[5] * 4411907939086405018UL) + ((uint64_t)op[6] * 6824475643257457942UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 6824475643257457942UL) + ((uint64_t)op[1] * 2803642121071829719UL) + ((uint64_t)op[2] * 6478878482505650402UL) + ((uint64_t)op[3] * 3285739491163352768UL) + ((uint64_t)op[4] * 436304521616961527UL) + ((((uint64_t)op[5] * 15854511070809174569UL) + ((uint64_t)op[6] * 4411907939086405018UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 4411907939086405018UL) + ((uint64_t)op[1] * 6824475643257457942UL) + ((uint64_t)op[2] * 2803642121071829719UL) + ((uint64_t)op[3] * 6478878482505650402UL) + ((uint64_t)op[4] * 3285739491163352768UL) + ((uint64_t)op[5] * 436304521616961527UL) + ((uint64_t)op[6] * 301113053406912287UL);
	tmp_q[6] = ((uint64_t)op[0] * 15854511070809174569UL) + ((uint64_t)op[1] * 4411907939086405018UL) + ((uint64_t)op[2] * 6824475643257457942UL) + ((uint64_t)op[3] * 2803642121071829719UL) + ((uint64_t)op[4] * 6478878482505650402UL) + ((uint64_t)op[5] * 3285739491163352768UL) + ((uint64_t)op[6] * 436304521616961527UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 6430511055L) + ((-((int128)tmp_q[1] * 13203206117L) + ((int128)tmp_q[2] * 22644183692L) + ((int128)tmp_q[3] * 2294359463L) + ((int128)tmp_q[4] * 42955051662L) - ((int128)tmp_q[5] * 84890608715L) - ((int128)tmp_q[6] * 40304816285L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 40304816285L) + ((int128)tmp_q[1] * 6430511055L) + ((-((int128)tmp_q[2] * 13203206117L) + ((int128)tmp_q[3] * 22644183692L) + ((int128)tmp_q[4] * 2294359463L) + ((int128)tmp_q[5] * 42955051662L) - ((int128)tmp_q[6] * 84890608715L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 84890608715L) - ((int128)tmp_q[1] * 40304816285L) + ((int128)tmp_q[2] * 6430511055L) + ((-((int128)tmp_q[3] * 13203206117L) + ((int128)tmp_q[4] * 22644183692L) + ((int128)tmp_q[5] * 2294359463L) + ((int128)tmp_q[6] * 42955051662L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 42955051662L) - ((int128)tmp_q[1] * 84890608715L) - ((int128)tmp_q[2] * 40304816285L) + ((int128)tmp_q[3] * 6430511055L) + ((-((int128)tmp_q[4] * 13203206117L) + ((int128)tmp_q[5] * 22644183692L) + ((int128)tmp_q[6] * 2294359463L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 2294359463L) + ((int128)tmp_q[1] * 42955051662L) - ((int128)tmp_q[2] * 84890608715L) - ((int128)tmp_q[3] * 40304816285L) + ((int128)tmp_q[4] * 6430511055L) + ((-((int128)tmp_q[5] * 13203206117L) + ((int128)tmp_q[6] * 22644183692L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 22644183692L) + ((int128)tmp_q[1] * 2294359463L) + ((int128)tmp_q[2] * 42955051662L) - ((int128)tmp_q[3] * 84890608715L) - ((int128)tmp_q[4] * 40304816285L) + ((int128)tmp_q[5] * 6430511055L) - ((int128)tmp_q[6] * 92422442819L);
	tmp_zero[6] = -((int128)tmp_q[0] * 13203206117L) + ((int128)tmp_q[1] * 22644183692L) + ((int128)tmp_q[2] * 2294359463L) + ((int128)tmp_q[3] * 42955051662L) - ((int128)tmp_q[4] * 84890608715L) - ((int128)tmp_q[5] * 40304816285L) + ((int128)tmp_q[6] * 6430511055L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

