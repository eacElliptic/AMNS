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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2398683323897937787UL) + ((((uint64_t)op[1] * 3325589521812522590UL) + ((uint64_t)op[2] * 1452362891708725607UL) + ((uint64_t)op[3] * 4360519530014489405UL) + ((uint64_t)op[4] * 5930275427978093131UL) + ((uint64_t)op[5] * 3939042262855703057UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 3939042262855703057UL) + ((uint64_t)op[1] * 2398683323897937787UL) + ((((uint64_t)op[2] * 3325589521812522590UL) + ((uint64_t)op[3] * 1452362891708725607UL) + ((uint64_t)op[4] * 4360519530014489405UL) + ((uint64_t)op[5] * 5930275427978093131UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 5930275427978093131UL) + ((uint64_t)op[1] * 3939042262855703057UL) + ((uint64_t)op[2] * 2398683323897937787UL) + ((((uint64_t)op[3] * 3325589521812522590UL) + ((uint64_t)op[4] * 1452362891708725607UL) + ((uint64_t)op[5] * 4360519530014489405UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 4360519530014489405UL) + ((uint64_t)op[1] * 5930275427978093131UL) + ((uint64_t)op[2] * 3939042262855703057UL) + ((uint64_t)op[3] * 2398683323897937787UL) + ((((uint64_t)op[4] * 3325589521812522590UL) + ((uint64_t)op[5] * 1452362891708725607UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 1452362891708725607UL) + ((uint64_t)op[1] * 4360519530014489405UL) + ((uint64_t)op[2] * 5930275427978093131UL) + ((uint64_t)op[3] * 3939042262855703057UL) + ((uint64_t)op[4] * 2398683323897937787UL) + ((uint64_t)op[5] * 11795565030084506436UL);
	tmp_q[5] = ((uint64_t)op[0] * 3325589521812522590UL) + ((uint64_t)op[1] * 1452362891708725607UL) + ((uint64_t)op[2] * 4360519530014489405UL) + ((uint64_t)op[3] * 5930275427978093131UL) + ((uint64_t)op[4] * 3939042262855703057UL) + ((uint64_t)op[5] * 2398683323897937787UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 10689727353L) - ((((int128)tmp_q[1] * 8026763167L) + ((int128)tmp_q[2] * 53851229614L) + ((int128)tmp_q[3] * 46429449824L) + ((int128)tmp_q[4] * 8748860288L) - ((int128)tmp_q[5] * 129340227731L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 129340227731L) - ((int128)tmp_q[1] * 10689727353L) - ((((int128)tmp_q[2] * 8026763167L) + ((int128)tmp_q[3] * 53851229614L) + ((int128)tmp_q[4] * 46429449824L) + ((int128)tmp_q[5] * 8748860288L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 8748860288L) - ((int128)tmp_q[1] * 129340227731L) - ((int128)tmp_q[2] * 10689727353L) - ((((int128)tmp_q[3] * 8026763167L) + ((int128)tmp_q[4] * 53851229614L) + ((int128)tmp_q[5] * 46429449824L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 46429449824L) + ((int128)tmp_q[1] * 8748860288L) - ((int128)tmp_q[2] * 129340227731L) - ((int128)tmp_q[3] * 10689727353L) - ((((int128)tmp_q[4] * 8026763167L) + ((int128)tmp_q[5] * 53851229614L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 53851229614L) + ((int128)tmp_q[1] * 46429449824L) + ((int128)tmp_q[2] * 8748860288L) - ((int128)tmp_q[3] * 129340227731L) - ((int128)tmp_q[4] * 10689727353L) - ((int128)tmp_q[5] * 16053526334L);
	tmp_zero[5] = ((int128)tmp_q[0] * 8026763167L) + ((int128)tmp_q[1] * 53851229614L) + ((int128)tmp_q[2] * 46429449824L) + ((int128)tmp_q[3] * 8748860288L) - ((int128)tmp_q[4] * 129340227731L) - ((int128)tmp_q[5] * 10689727353L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

