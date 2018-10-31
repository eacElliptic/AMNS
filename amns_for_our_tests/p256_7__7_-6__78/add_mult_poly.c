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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16717455452056804911UL) + ((((uint64_t)op[1] * 17480933949076899427UL) + ((uint64_t)op[2] * 5244302127971129157UL) + ((uint64_t)op[3] * 5109464532809255396UL) + ((uint64_t)op[4] * 9804440194430136299UL) + ((uint64_t)op[5] * 13665733474261823351UL) + ((uint64_t)op[6] * 17330850486005049771UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 17330850486005049771UL) + ((uint64_t)op[1] * 16717455452056804911UL) + ((((uint64_t)op[2] * 17480933949076899427UL) + ((uint64_t)op[3] * 5244302127971129157UL) + ((uint64_t)op[4] * 5109464532809255396UL) + ((uint64_t)op[5] * 9804440194430136299UL) + ((uint64_t)op[6] * 13665733474261823351UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 13665733474261823351UL) + ((uint64_t)op[1] * 17330850486005049771UL) + ((uint64_t)op[2] * 16717455452056804911UL) + ((((uint64_t)op[3] * 17480933949076899427UL) + ((uint64_t)op[4] * 5244302127971129157UL) + ((uint64_t)op[5] * 5109464532809255396UL) + ((uint64_t)op[6] * 9804440194430136299UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 9804440194430136299UL) + ((uint64_t)op[1] * 13665733474261823351UL) + ((uint64_t)op[2] * 17330850486005049771UL) + ((uint64_t)op[3] * 16717455452056804911UL) + ((((uint64_t)op[4] * 17480933949076899427UL) + ((uint64_t)op[5] * 5244302127971129157UL) + ((uint64_t)op[6] * 5109464532809255396UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 5109464532809255396UL) + ((uint64_t)op[1] * 9804440194430136299UL) + ((uint64_t)op[2] * 13665733474261823351UL) + ((uint64_t)op[3] * 17330850486005049771UL) + ((uint64_t)op[4] * 16717455452056804911UL) + ((((uint64_t)op[5] * 17480933949076899427UL) + ((uint64_t)op[6] * 5244302127971129157UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 5244302127971129157UL) + ((uint64_t)op[1] * 5109464532809255396UL) + ((uint64_t)op[2] * 9804440194430136299UL) + ((uint64_t)op[3] * 13665733474261823351UL) + ((uint64_t)op[4] * 17330850486005049771UL) + ((uint64_t)op[5] * 16717455452056804911UL) + ((uint64_t)op[6] * 5794860747795913134UL);
	tmp_q[6] = ((uint64_t)op[0] * 17480933949076899427UL) + ((uint64_t)op[1] * 5244302127971129157UL) + ((uint64_t)op[2] * 5109464532809255396UL) + ((uint64_t)op[3] * 9804440194430136299UL) + ((uint64_t)op[4] * 13665733474261823351UL) + ((uint64_t)op[5] * 17330850486005049771UL) + ((uint64_t)op[6] * 16717455452056804911UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 48229084915L) - ((-((int128)tmp_q[1] * 15919486151L) - ((int128)tmp_q[2] * 30282337766L) + ((int128)tmp_q[3] * 8122566903L) - ((int128)tmp_q[4] * 62913494970L) + ((int128)tmp_q[5] * 38286359102L) + ((int128)tmp_q[6] * 23862956859L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 23862956859L) + ((int128)tmp_q[1] * 48229084915L) - ((-((int128)tmp_q[2] * 15919486151L) - ((int128)tmp_q[3] * 30282337766L) + ((int128)tmp_q[4] * 8122566903L) - ((int128)tmp_q[5] * 62913494970L) + ((int128)tmp_q[6] * 38286359102L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 38286359102L) + ((int128)tmp_q[1] * 23862956859L) + ((int128)tmp_q[2] * 48229084915L) - ((-((int128)tmp_q[3] * 15919486151L) - ((int128)tmp_q[4] * 30282337766L) + ((int128)tmp_q[5] * 8122566903L) - ((int128)tmp_q[6] * 62913494970L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 62913494970L) + ((int128)tmp_q[1] * 38286359102L) + ((int128)tmp_q[2] * 23862956859L) + ((int128)tmp_q[3] * 48229084915L) - ((-((int128)tmp_q[4] * 15919486151L) - ((int128)tmp_q[5] * 30282337766L) + ((int128)tmp_q[6] * 8122566903L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 8122566903L) - ((int128)tmp_q[1] * 62913494970L) + ((int128)tmp_q[2] * 38286359102L) + ((int128)tmp_q[3] * 23862956859L) + ((int128)tmp_q[4] * 48229084915L) - ((-((int128)tmp_q[5] * 15919486151L) - ((int128)tmp_q[6] * 30282337766L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 30282337766L) + ((int128)tmp_q[1] * 8122566903L) - ((int128)tmp_q[2] * 62913494970L) + ((int128)tmp_q[3] * 38286359102L) + ((int128)tmp_q[4] * 23862956859L) + ((int128)tmp_q[5] * 48229084915L) + ((int128)tmp_q[6] * 95516916906L);
	tmp_zero[6] = -((int128)tmp_q[0] * 15919486151L) - ((int128)tmp_q[1] * 30282337766L) + ((int128)tmp_q[2] * 8122566903L) - ((int128)tmp_q[3] * 62913494970L) + ((int128)tmp_q[4] * 38286359102L) + ((int128)tmp_q[5] * 23862956859L) + ((int128)tmp_q[6] * 48229084915L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

