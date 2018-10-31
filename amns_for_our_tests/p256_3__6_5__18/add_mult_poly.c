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
	tmp_q[0] = ((uint64_t)op[0] * 15394554892283549752UL) + ((((uint64_t)op[1] * 3125428820988755517UL) + ((uint64_t)op[2] * 1592867291347948557UL) + ((uint64_t)op[3] * 10694720006058051538UL) + ((uint64_t)op[4] * 15563183728542233481UL) + ((uint64_t)op[5] * 14968472398854843700UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 14968472398854843700UL) + ((uint64_t)op[1] * 15394554892283549752UL) + ((((uint64_t)op[2] * 3125428820988755517UL) + ((uint64_t)op[3] * 1592867291347948557UL) + ((uint64_t)op[4] * 10694720006058051538UL) + ((uint64_t)op[5] * 15563183728542233481UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 15563183728542233481UL) + ((uint64_t)op[1] * 14968472398854843700UL) + ((uint64_t)op[2] * 15394554892283549752UL) + ((((uint64_t)op[3] * 3125428820988755517UL) + ((uint64_t)op[4] * 1592867291347948557UL) + ((uint64_t)op[5] * 10694720006058051538UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 10694720006058051538UL) + ((uint64_t)op[1] * 15563183728542233481UL) + ((uint64_t)op[2] * 14968472398854843700UL) + ((uint64_t)op[3] * 15394554892283549752UL) + ((((uint64_t)op[4] * 3125428820988755517UL) + ((uint64_t)op[5] * 1592867291347948557UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 1592867291347948557UL) + ((uint64_t)op[1] * 10694720006058051538UL) + ((uint64_t)op[2] * 15563183728542233481UL) + ((uint64_t)op[3] * 14968472398854843700UL) + ((uint64_t)op[4] * 15394554892283549752UL) + ((uint64_t)op[5] * 15627144104943777585UL);
	tmp_q[5] = ((uint64_t)op[0] * 3125428820988755517UL) + ((uint64_t)op[1] * 1592867291347948557UL) + ((uint64_t)op[2] * 10694720006058051538UL) + ((uint64_t)op[3] * 15563183728542233481UL) + ((uint64_t)op[4] * 14968472398854843700UL) + ((uint64_t)op[5] * 15394554892283549752UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3992928879659L) + ((-((int128)tmp_q[1] * 3716359481576L) - ((int128)tmp_q[2] * 3798964114570L) - ((int128)tmp_q[3] * 1344380224971L) + ((int128)tmp_q[4] * 1954684911267L) + ((int128)tmp_q[5] * 1831752249284L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1831752249284L) - ((int128)tmp_q[1] * 3992928879659L) + ((-((int128)tmp_q[2] * 3716359481576L) - ((int128)tmp_q[3] * 3798964114570L) - ((int128)tmp_q[4] * 1344380224971L) + ((int128)tmp_q[5] * 1954684911267L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 1954684911267L) + ((int128)tmp_q[1] * 1831752249284L) - ((int128)tmp_q[2] * 3992928879659L) + ((-((int128)tmp_q[3] * 3716359481576L) - ((int128)tmp_q[4] * 3798964114570L) - ((int128)tmp_q[5] * 1344380224971L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 1344380224971L) + ((int128)tmp_q[1] * 1954684911267L) + ((int128)tmp_q[2] * 1831752249284L) - ((int128)tmp_q[3] * 3992928879659L) + ((-((int128)tmp_q[4] * 3716359481576L) - ((int128)tmp_q[5] * 3798964114570L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 3798964114570L) - ((int128)tmp_q[1] * 1344380224971L) + ((int128)tmp_q[2] * 1954684911267L) + ((int128)tmp_q[3] * 1831752249284L) - ((int128)tmp_q[4] * 3992928879659L) - ((int128)tmp_q[5] * 18581797407880L);
	tmp_zero[5] = -((int128)tmp_q[0] * 3716359481576L) - ((int128)tmp_q[1] * 3798964114570L) - ((int128)tmp_q[2] * 1344380224971L) + ((int128)tmp_q[3] * 1954684911267L) + ((int128)tmp_q[4] * 1831752249284L) - ((int128)tmp_q[5] * 3992928879659L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

