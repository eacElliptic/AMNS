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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8167240352869218398UL) + ((((uint64_t)op[1] * 14452537628054282151UL) + ((uint64_t)op[2] * 17018313908731043829UL) + ((uint64_t)op[3] * 15384007859276848942UL) + ((uint64_t)op[4] * 5257667048907356811UL) + ((uint64_t)op[5] * 8250774309030996012UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 8250774309030996012UL) + ((uint64_t)op[1] * 8167240352869218398UL) + ((((uint64_t)op[2] * 14452537628054282151UL) + ((uint64_t)op[3] * 17018313908731043829UL) + ((uint64_t)op[4] * 15384007859276848942UL) + ((uint64_t)op[5] * 5257667048907356811UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 5257667048907356811UL) + ((uint64_t)op[1] * 8250774309030996012UL) + ((uint64_t)op[2] * 8167240352869218398UL) + ((((uint64_t)op[3] * 14452537628054282151UL) + ((uint64_t)op[4] * 17018313908731043829UL) + ((uint64_t)op[5] * 15384007859276848942UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 15384007859276848942UL) + ((uint64_t)op[1] * 5257667048907356811UL) + ((uint64_t)op[2] * 8250774309030996012UL) + ((uint64_t)op[3] * 8167240352869218398UL) + ((((uint64_t)op[4] * 14452537628054282151UL) + ((uint64_t)op[5] * 17018313908731043829UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 17018313908731043829UL) + ((uint64_t)op[1] * 15384007859276848942UL) + ((uint64_t)op[2] * 5257667048907356811UL) + ((uint64_t)op[3] * 8250774309030996012UL) + ((uint64_t)op[4] * 8167240352869218398UL) + ((uint64_t)op[5] * 16922455919142755907UL);
	tmp_q[5] = ((uint64_t)op[0] * 14452537628054282151UL) + ((uint64_t)op[1] * 17018313908731043829UL) + ((uint64_t)op[2] * 15384007859276848942UL) + ((uint64_t)op[3] * 5257667048907356811UL) + ((uint64_t)op[4] * 8250774309030996012UL) + ((uint64_t)op[5] * 8167240352869218398UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2012105877919L) + ((((int128)tmp_q[1] * 4436422037716L) - ((int128)tmp_q[2] * 2531916779856L) - ((int128)tmp_q[3] * 4395446986933L) + ((int128)tmp_q[4] * 2773266022519L) + ((int128)tmp_q[5] * 1268845965636L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1268845965636L) + ((int128)tmp_q[1] * 2012105877919L) + ((((int128)tmp_q[2] * 4436422037716L) - ((int128)tmp_q[3] * 2531916779856L) - ((int128)tmp_q[4] * 4395446986933L) + ((int128)tmp_q[5] * 2773266022519L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 2773266022519L) + ((int128)tmp_q[1] * 1268845965636L) + ((int128)tmp_q[2] * 2012105877919L) + ((((int128)tmp_q[3] * 4436422037716L) - ((int128)tmp_q[4] * 2531916779856L) - ((int128)tmp_q[5] * 4395446986933L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 4395446986933L) + ((int128)tmp_q[1] * 2773266022519L) + ((int128)tmp_q[2] * 1268845965636L) + ((int128)tmp_q[3] * 2012105877919L) + ((((int128)tmp_q[4] * 4436422037716L) - ((int128)tmp_q[5] * 2531916779856L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 2531916779856L) - ((int128)tmp_q[1] * 4395446986933L) + ((int128)tmp_q[2] * 2773266022519L) + ((int128)tmp_q[3] * 1268845965636L) + ((int128)tmp_q[4] * 2012105877919L) + ((int128)tmp_q[5] * 22182110188580L);
	tmp_zero[5] = ((int128)tmp_q[0] * 4436422037716L) - ((int128)tmp_q[1] * 2531916779856L) - ((int128)tmp_q[2] * 4395446986933L) + ((int128)tmp_q[3] * 2773266022519L) + ((int128)tmp_q[4] * 1268845965636L) + ((int128)tmp_q[5] * 2012105877919L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

