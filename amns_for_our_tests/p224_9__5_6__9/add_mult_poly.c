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
	tmp_q[0] = ((uint64_t)op[0] * 4404954920573523943UL) + ((((uint64_t)op[1] * 15809529719067946918UL) + ((uint64_t)op[2] * 6482987377102316880UL) + ((uint64_t)op[3] * 7670464573074097204UL) + ((uint64_t)op[4] * 13506383178365550079UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 13506383178365550079UL) + ((uint64_t)op[1] * 4404954920573523943UL) + ((((uint64_t)op[2] * 15809529719067946918UL) + ((uint64_t)op[3] * 6482987377102316880UL) + ((uint64_t)op[4] * 7670464573074097204UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 7670464573074097204UL) + ((uint64_t)op[1] * 13506383178365550079UL) + ((uint64_t)op[2] * 4404954920573523943UL) + ((((uint64_t)op[3] * 15809529719067946918UL) + ((uint64_t)op[4] * 6482987377102316880UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 6482987377102316880UL) + ((uint64_t)op[1] * 7670464573074097204UL) + ((uint64_t)op[2] * 13506383178365550079UL) + ((uint64_t)op[3] * 4404954920573523943UL) + ((uint64_t)op[4] * 2623457945859923428UL);
	tmp_q[4] = ((uint64_t)op[0] * 15809529719067946918UL) + ((uint64_t)op[1] * 6482987377102316880UL) + ((uint64_t)op[2] * 7670464573074097204UL) + ((uint64_t)op[3] * 13506383178365550079UL) + ((uint64_t)op[4] * 4404954920573523943UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 12674026108719L) + ((-((int128)tmp_q[1] * 14461808195863L) - ((int128)tmp_q[2] * 27547200706747L) + ((int128)tmp_q[3] * 15357811003355L) + ((int128)tmp_q[4] * 2313306341501L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 2313306341501L) + ((int128)tmp_q[1] * 12674026108719L) + ((-((int128)tmp_q[2] * 14461808195863L) - ((int128)tmp_q[3] * 27547200706747L) + ((int128)tmp_q[4] * 15357811003355L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 15357811003355L) + ((int128)tmp_q[1] * 2313306341501L) + ((int128)tmp_q[2] * 12674026108719L) + ((-((int128)tmp_q[3] * 14461808195863L) - ((int128)tmp_q[4] * 27547200706747L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 27547200706747L) + ((int128)tmp_q[1] * 15357811003355L) + ((int128)tmp_q[2] * 2313306341501L) + ((int128)tmp_q[3] * 12674026108719L) - ((int128)tmp_q[4] * 86770849175178L);
	tmp_zero[4] = -((int128)tmp_q[0] * 14461808195863L) - ((int128)tmp_q[1] * 27547200706747L) + ((int128)tmp_q[2] * 15357811003355L) + ((int128)tmp_q[3] * 2313306341501L) + ((int128)tmp_q[4] * 12674026108719L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

