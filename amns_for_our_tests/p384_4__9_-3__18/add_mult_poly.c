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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17416945808973117467UL) + ((((uint64_t)op[1] * 15105047561194622506UL) + ((uint64_t)op[2] * 931767206095425046UL) + ((uint64_t)op[3] * 9413610373667466946UL) + ((uint64_t)op[4] * 189560475937767962UL) + ((uint64_t)op[5] * 2403526764199910087UL) + ((uint64_t)op[6] * 12882045241393250821UL) + ((uint64_t)op[7] * 136675690194236212UL) + ((uint64_t)op[8] * 3533001591371269778UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 3533001591371269778UL) + ((uint64_t)op[1] * 17416945808973117467UL) + ((((uint64_t)op[2] * 15105047561194622506UL) + ((uint64_t)op[3] * 931767206095425046UL) + ((uint64_t)op[4] * 9413610373667466946UL) + ((uint64_t)op[5] * 189560475937767962UL) + ((uint64_t)op[6] * 2403526764199910087UL) + ((uint64_t)op[7] * 12882045241393250821UL) + ((uint64_t)op[8] * 136675690194236212UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 136675690194236212UL) + ((uint64_t)op[1] * 3533001591371269778UL) + ((uint64_t)op[2] * 17416945808973117467UL) + ((((uint64_t)op[3] * 15105047561194622506UL) + ((uint64_t)op[4] * 931767206095425046UL) + ((uint64_t)op[5] * 9413610373667466946UL) + ((uint64_t)op[6] * 189560475937767962UL) + ((uint64_t)op[7] * 2403526764199910087UL) + ((uint64_t)op[8] * 12882045241393250821UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 12882045241393250821UL) + ((uint64_t)op[1] * 136675690194236212UL) + ((uint64_t)op[2] * 3533001591371269778UL) + ((uint64_t)op[3] * 17416945808973117467UL) + ((((uint64_t)op[4] * 15105047561194622506UL) + ((uint64_t)op[5] * 931767206095425046UL) + ((uint64_t)op[6] * 9413610373667466946UL) + ((uint64_t)op[7] * 189560475937767962UL) + ((uint64_t)op[8] * 2403526764199910087UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 2403526764199910087UL) + ((uint64_t)op[1] * 12882045241393250821UL) + ((uint64_t)op[2] * 136675690194236212UL) + ((uint64_t)op[3] * 3533001591371269778UL) + ((uint64_t)op[4] * 17416945808973117467UL) + ((((uint64_t)op[5] * 15105047561194622506UL) + ((uint64_t)op[6] * 931767206095425046UL) + ((uint64_t)op[7] * 9413610373667466946UL) + ((uint64_t)op[8] * 189560475937767962UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 189560475937767962UL) + ((uint64_t)op[1] * 2403526764199910087UL) + ((uint64_t)op[2] * 12882045241393250821UL) + ((uint64_t)op[3] * 136675690194236212UL) + ((uint64_t)op[4] * 3533001591371269778UL) + ((uint64_t)op[5] * 17416945808973117467UL) + ((((uint64_t)op[6] * 15105047561194622506UL) + ((uint64_t)op[7] * 931767206095425046UL) + ((uint64_t)op[8] * 9413610373667466946UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 9413610373667466946UL) + ((uint64_t)op[1] * 189560475937767962UL) + ((uint64_t)op[2] * 2403526764199910087UL) + ((uint64_t)op[3] * 12882045241393250821UL) + ((uint64_t)op[4] * 136675690194236212UL) + ((uint64_t)op[5] * 3533001591371269778UL) + ((uint64_t)op[6] * 17416945808973117467UL) + ((((uint64_t)op[7] * 15105047561194622506UL) + ((uint64_t)op[8] * 931767206095425046UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 931767206095425046UL) + ((uint64_t)op[1] * 9413610373667466946UL) + ((uint64_t)op[2] * 189560475937767962UL) + ((uint64_t)op[3] * 2403526764199910087UL) + ((uint64_t)op[4] * 12882045241393250821UL) + ((uint64_t)op[5] * 136675690194236212UL) + ((uint64_t)op[6] * 3533001591371269778UL) + ((uint64_t)op[7] * 17416945808973117467UL) + ((uint64_t)op[8] * 10025089537544787330UL);
	tmp_q[8] = ((uint64_t)op[0] * 15105047561194622506UL) + ((uint64_t)op[1] * 931767206095425046UL) + ((uint64_t)op[2] * 9413610373667466946UL) + ((uint64_t)op[3] * 189560475937767962UL) + ((uint64_t)op[4] * 2403526764199910087UL) + ((uint64_t)op[5] * 12882045241393250821UL) + ((uint64_t)op[6] * 136675690194236212UL) + ((uint64_t)op[7] * 3533001591371269778UL) + ((uint64_t)op[8] * 17416945808973117467UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4537791832463L) - ((-((int128)tmp_q[1] * 734735230194L) + ((int128)tmp_q[2] * 736225682273L) + ((int128)tmp_q[3] * 150935145204L) - ((int128)tmp_q[4] * 1473718638666L) - ((int128)tmp_q[5] * 3397755734282L) - ((int128)tmp_q[6] * 4589056541373L) - ((int128)tmp_q[7] * 3366693268747L) - ((int128)tmp_q[8] * 1655139725081L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 1655139725081L) - ((int128)tmp_q[1] * 4537791832463L) - ((-((int128)tmp_q[2] * 734735230194L) + ((int128)tmp_q[3] * 736225682273L) + ((int128)tmp_q[4] * 150935145204L) - ((int128)tmp_q[5] * 1473718638666L) - ((int128)tmp_q[6] * 3397755734282L) - ((int128)tmp_q[7] * 4589056541373L) - ((int128)tmp_q[8] * 3366693268747L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 3366693268747L) - ((int128)tmp_q[1] * 1655139725081L) - ((int128)tmp_q[2] * 4537791832463L) - ((-((int128)tmp_q[3] * 734735230194L) + ((int128)tmp_q[4] * 736225682273L) + ((int128)tmp_q[5] * 150935145204L) - ((int128)tmp_q[6] * 1473718638666L) - ((int128)tmp_q[7] * 3397755734282L) - ((int128)tmp_q[8] * 4589056541373L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 4589056541373L) - ((int128)tmp_q[1] * 3366693268747L) - ((int128)tmp_q[2] * 1655139725081L) - ((int128)tmp_q[3] * 4537791832463L) - ((-((int128)tmp_q[4] * 734735230194L) + ((int128)tmp_q[5] * 736225682273L) + ((int128)tmp_q[6] * 150935145204L) - ((int128)tmp_q[7] * 1473718638666L) - ((int128)tmp_q[8] * 3397755734282L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 3397755734282L) - ((int128)tmp_q[1] * 4589056541373L) - ((int128)tmp_q[2] * 3366693268747L) - ((int128)tmp_q[3] * 1655139725081L) - ((int128)tmp_q[4] * 4537791832463L) - ((-((int128)tmp_q[5] * 734735230194L) + ((int128)tmp_q[6] * 736225682273L) + ((int128)tmp_q[7] * 150935145204L) - ((int128)tmp_q[8] * 1473718638666L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 1473718638666L) - ((int128)tmp_q[1] * 3397755734282L) - ((int128)tmp_q[2] * 4589056541373L) - ((int128)tmp_q[3] * 3366693268747L) - ((int128)tmp_q[4] * 1655139725081L) - ((int128)tmp_q[5] * 4537791832463L) - ((-((int128)tmp_q[6] * 734735230194L) + ((int128)tmp_q[7] * 736225682273L) + ((int128)tmp_q[8] * 150935145204L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 150935145204L) - ((int128)tmp_q[1] * 1473718638666L) - ((int128)tmp_q[2] * 3397755734282L) - ((int128)tmp_q[3] * 4589056541373L) - ((int128)tmp_q[4] * 3366693268747L) - ((int128)tmp_q[5] * 1655139725081L) - ((int128)tmp_q[6] * 4537791832463L) - ((-((int128)tmp_q[7] * 734735230194L) + ((int128)tmp_q[8] * 736225682273L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 736225682273L) + ((int128)tmp_q[1] * 150935145204L) - ((int128)tmp_q[2] * 1473718638666L) - ((int128)tmp_q[3] * 3397755734282L) - ((int128)tmp_q[4] * 4589056541373L) - ((int128)tmp_q[5] * 3366693268747L) - ((int128)tmp_q[6] * 1655139725081L) - ((int128)tmp_q[7] * 4537791832463L) + ((int128)tmp_q[8] * 2204205690582L);
	tmp_zero[8] = -((int128)tmp_q[0] * 734735230194L) + ((int128)tmp_q[1] * 736225682273L) + ((int128)tmp_q[2] * 150935145204L) - ((int128)tmp_q[3] * 1473718638666L) - ((int128)tmp_q[4] * 3397755734282L) - ((int128)tmp_q[5] * 4589056541373L) - ((int128)tmp_q[6] * 3366693268747L) - ((int128)tmp_q[7] * 1655139725081L) - ((int128)tmp_q[8] * 4537791832463L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
}

