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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5611418817079275215UL) + ((((uint64_t)op[1] * 2178275878485539606UL) + ((uint64_t)op[2] * 16622380109279607139UL) + ((uint64_t)op[3] * 13319290999106622943UL) + ((uint64_t)op[4] * 6668789490752282506UL) + ((uint64_t)op[5] * 4409162516987859136UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 4409162516987859136UL) + ((uint64_t)op[1] * 5611418817079275215UL) + ((((uint64_t)op[2] * 2178275878485539606UL) + ((uint64_t)op[3] * 16622380109279607139UL) + ((uint64_t)op[4] * 13319290999106622943UL) + ((uint64_t)op[5] * 6668789490752282506UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 6668789490752282506UL) + ((uint64_t)op[1] * 4409162516987859136UL) + ((uint64_t)op[2] * 5611418817079275215UL) + ((((uint64_t)op[3] * 2178275878485539606UL) + ((uint64_t)op[4] * 16622380109279607139UL) + ((uint64_t)op[5] * 13319290999106622943UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 13319290999106622943UL) + ((uint64_t)op[1] * 6668789490752282506UL) + ((uint64_t)op[2] * 4409162516987859136UL) + ((uint64_t)op[3] * 5611418817079275215UL) + ((((uint64_t)op[4] * 2178275878485539606UL) + ((uint64_t)op[5] * 16622380109279607139UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 16622380109279607139UL) + ((uint64_t)op[1] * 13319290999106622943UL) + ((uint64_t)op[2] * 6668789490752282506UL) + ((uint64_t)op[3] * 4409162516987859136UL) + ((uint64_t)op[4] * 5611418817079275215UL) + ((uint64_t)op[5] * 4356551756971079212UL);
	tmp_q[5] = ((uint64_t)op[0] * 2178275878485539606UL) + ((uint64_t)op[1] * 16622380109279607139UL) + ((uint64_t)op[2] * 13319290999106622943UL) + ((uint64_t)op[3] * 6668789490752282506UL) + ((uint64_t)op[4] * 4409162516987859136UL) + ((uint64_t)op[5] * 5611418817079275215UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2260775201139L) + ((-((int128)tmp_q[1] * 4545062891320L) - ((int128)tmp_q[2] * 2059624301427L) - ((int128)tmp_q[3] * 2855358254431L) - ((int128)tmp_q[4] * 3120885660760L) + ((int128)tmp_q[5] * 2323887082668L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 2323887082668L) + ((int128)tmp_q[1] * 2260775201139L) + ((-((int128)tmp_q[2] * 4545062891320L) - ((int128)tmp_q[3] * 2059624301427L) - ((int128)tmp_q[4] * 2855358254431L) - ((int128)tmp_q[5] * 3120885660760L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 3120885660760L) + ((int128)tmp_q[1] * 2323887082668L) + ((int128)tmp_q[2] * 2260775201139L) + ((-((int128)tmp_q[3] * 4545062891320L) - ((int128)tmp_q[4] * 2059624301427L) - ((int128)tmp_q[5] * 2855358254431L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 2855358254431L) - ((int128)tmp_q[1] * 3120885660760L) + ((int128)tmp_q[2] * 2323887082668L) + ((int128)tmp_q[3] * 2260775201139L) + ((-((int128)tmp_q[4] * 4545062891320L) - ((int128)tmp_q[5] * 2059624301427L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 2059624301427L) - ((int128)tmp_q[1] * 2855358254431L) - ((int128)tmp_q[2] * 3120885660760L) + ((int128)tmp_q[3] * 2323887082668L) + ((int128)tmp_q[4] * 2260775201139L) - ((int128)tmp_q[5] * 9090125782640L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4545062891320L) - ((int128)tmp_q[1] * 2059624301427L) - ((int128)tmp_q[2] * 2855358254431L) - ((int128)tmp_q[3] * 3120885660760L) + ((int128)tmp_q[4] * 2323887082668L) + ((int128)tmp_q[5] * 2260775201139L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

