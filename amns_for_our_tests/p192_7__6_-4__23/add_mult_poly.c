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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13765450266572920317UL) + ((((uint64_t)op[1] * 8507905803171335591UL) + ((uint64_t)op[2] * 8987529842996372382UL) + ((uint64_t)op[3] * 16366639998270493236UL) + ((uint64_t)op[4] * 7227271588011792585UL) + ((uint64_t)op[5] * 15257808896920616783UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 15257808896920616783UL) + ((uint64_t)op[1] * 13765450266572920317UL) + ((((uint64_t)op[2] * 8507905803171335591UL) + ((uint64_t)op[3] * 8987529842996372382UL) + ((uint64_t)op[4] * 16366639998270493236UL) + ((uint64_t)op[5] * 7227271588011792585UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 7227271588011792585UL) + ((uint64_t)op[1] * 15257808896920616783UL) + ((uint64_t)op[2] * 13765450266572920317UL) + ((((uint64_t)op[3] * 8507905803171335591UL) + ((uint64_t)op[4] * 8987529842996372382UL) + ((uint64_t)op[5] * 16366639998270493236UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 16366639998270493236UL) + ((uint64_t)op[1] * 7227271588011792585UL) + ((uint64_t)op[2] * 15257808896920616783UL) + ((uint64_t)op[3] * 13765450266572920317UL) + ((((uint64_t)op[4] * 8507905803171335591UL) + ((uint64_t)op[5] * 8987529842996372382UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 8987529842996372382UL) + ((uint64_t)op[1] * 16366639998270493236UL) + ((uint64_t)op[2] * 7227271588011792585UL) + ((uint64_t)op[3] * 15257808896920616783UL) + ((uint64_t)op[4] * 13765450266572920317UL) + ((uint64_t)op[5] * 2861864934733760868UL);
	tmp_q[5] = ((uint64_t)op[0] * 8507905803171335591UL) + ((uint64_t)op[1] * 8987529842996372382UL) + ((uint64_t)op[2] * 16366639998270493236UL) + ((uint64_t)op[3] * 7227271588011792585UL) + ((uint64_t)op[4] * 15257808896920616783UL) + ((uint64_t)op[5] * 13765450266572920317UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1062225745L) - ((((int128)tmp_q[1] * 1188203055L) + ((int128)tmp_q[2] * 1460932319L) + ((int128)tmp_q[3] * 1489752985L) - ((int128)tmp_q[4] * 282719460L) - ((int128)tmp_q[5] * 1453960753L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 1453960753L) - ((int128)tmp_q[1] * 1062225745L) - ((((int128)tmp_q[2] * 1188203055L) + ((int128)tmp_q[3] * 1460932319L) + ((int128)tmp_q[4] * 1489752985L) - ((int128)tmp_q[5] * 282719460L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 282719460L) - ((int128)tmp_q[1] * 1453960753L) - ((int128)tmp_q[2] * 1062225745L) - ((((int128)tmp_q[3] * 1188203055L) + ((int128)tmp_q[4] * 1460932319L) + ((int128)tmp_q[5] * 1489752985L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 1489752985L) - ((int128)tmp_q[1] * 282719460L) - ((int128)tmp_q[2] * 1453960753L) - ((int128)tmp_q[3] * 1062225745L) - ((((int128)tmp_q[4] * 1188203055L) + ((int128)tmp_q[5] * 1460932319L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 1460932319L) + ((int128)tmp_q[1] * 1489752985L) - ((int128)tmp_q[2] * 282719460L) - ((int128)tmp_q[3] * 1453960753L) - ((int128)tmp_q[4] * 1062225745L) - ((int128)tmp_q[5] * 4752812220L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1188203055L) + ((int128)tmp_q[1] * 1460932319L) + ((int128)tmp_q[2] * 1489752985L) - ((int128)tmp_q[3] * 282719460L) - ((int128)tmp_q[4] * 1453960753L) - ((int128)tmp_q[5] * 1062225745L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

