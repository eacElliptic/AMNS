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
	tmp_q[0] = ((uint64_t)op[0] * 13768276089251351745UL) + ((((uint64_t)op[1] * 11173320383993864848UL) + ((uint64_t)op[2] * 18439672511719094859UL) + ((uint64_t)op[3] * 744394751724029273UL) + ((uint64_t)op[4] * 5617173083304154909UL) + ((uint64_t)op[5] * 8066395316739604753UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 8066395316739604753UL) + ((uint64_t)op[1] * 13768276089251351745UL) + ((((uint64_t)op[2] * 11173320383993864848UL) + ((uint64_t)op[3] * 18439672511719094859UL) + ((uint64_t)op[4] * 744394751724029273UL) + ((uint64_t)op[5] * 5617173083304154909UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 5617173083304154909UL) + ((uint64_t)op[1] * 8066395316739604753UL) + ((uint64_t)op[2] * 13768276089251351745UL) + ((((uint64_t)op[3] * 11173320383993864848UL) + ((uint64_t)op[4] * 18439672511719094859UL) + ((uint64_t)op[5] * 744394751724029273UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 744394751724029273UL) + ((uint64_t)op[1] * 5617173083304154909UL) + ((uint64_t)op[2] * 8066395316739604753UL) + ((uint64_t)op[3] * 13768276089251351745UL) + ((((uint64_t)op[4] * 11173320383993864848UL) + ((uint64_t)op[5] * 18439672511719094859UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 18439672511719094859UL) + ((uint64_t)op[1] * 744394751724029273UL) + ((uint64_t)op[2] * 5617173083304154909UL) + ((uint64_t)op[3] * 8066395316739604753UL) + ((uint64_t)op[4] * 13768276089251351745UL) + ((uint64_t)op[5] * 14546847379431373536UL);
	tmp_q[5] = ((uint64_t)op[0] * 11173320383993864848UL) + ((uint64_t)op[1] * 18439672511719094859UL) + ((uint64_t)op[2] * 744394751724029273UL) + ((uint64_t)op[3] * 5617173083304154909UL) + ((uint64_t)op[4] * 8066395316739604753UL) + ((uint64_t)op[5] * 13768276089251351745UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3177478130855L) - ((((int128)tmp_q[1] * 881676945041L) - ((int128)tmp_q[2] * 2242363787748L) - ((int128)tmp_q[3] * 982025638892L) - ((int128)tmp_q[4] * 3476075588172L) + ((int128)tmp_q[5] * 1404160158177L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1404160158177L) - ((int128)tmp_q[1] * 3177478130855L) - ((((int128)tmp_q[2] * 881676945041L) - ((int128)tmp_q[3] * 2242363787748L) - ((int128)tmp_q[4] * 982025638892L) - ((int128)tmp_q[5] * 3476075588172L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 3476075588172L) + ((int128)tmp_q[1] * 1404160158177L) - ((int128)tmp_q[2] * 3177478130855L) - ((((int128)tmp_q[3] * 881676945041L) - ((int128)tmp_q[4] * 2242363787748L) - ((int128)tmp_q[5] * 982025638892L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 982025638892L) - ((int128)tmp_q[1] * 3476075588172L) + ((int128)tmp_q[2] * 1404160158177L) - ((int128)tmp_q[3] * 3177478130855L) - ((((int128)tmp_q[4] * 881676945041L) - ((int128)tmp_q[5] * 2242363787748L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 2242363787748L) - ((int128)tmp_q[1] * 982025638892L) - ((int128)tmp_q[2] * 3476075588172L) + ((int128)tmp_q[3] * 1404160158177L) - ((int128)tmp_q[4] * 3177478130855L) - ((int128)tmp_q[5] * 1763353890082L);
	tmp_zero[5] = ((int128)tmp_q[0] * 881676945041L) - ((int128)tmp_q[1] * 2242363787748L) - ((int128)tmp_q[2] * 982025638892L) - ((int128)tmp_q[3] * 3476075588172L) + ((int128)tmp_q[4] * 1404160158177L) - ((int128)tmp_q[5] * 3177478130855L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

