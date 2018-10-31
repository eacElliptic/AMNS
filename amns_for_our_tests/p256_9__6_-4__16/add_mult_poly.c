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
	tmp_q[0] = ((uint64_t)op[0] * 4743391980215843099UL) + ((((uint64_t)op[1] * 10151408009908531440UL) + ((uint64_t)op[2] * 496552594208016366UL) + ((uint64_t)op[3] * 11913621570039014168UL) + ((uint64_t)op[4] * 3713849246019306254UL) + ((uint64_t)op[5] * 12054912410056150093UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 12054912410056150093UL) + ((uint64_t)op[1] * 4743391980215843099UL) + ((((uint64_t)op[2] * 10151408009908531440UL) + ((uint64_t)op[3] * 496552594208016366UL) + ((uint64_t)op[4] * 11913621570039014168UL) + ((uint64_t)op[5] * 3713849246019306254UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 3713849246019306254UL) + ((uint64_t)op[1] * 12054912410056150093UL) + ((uint64_t)op[2] * 4743391980215843099UL) + ((((uint64_t)op[3] * 10151408009908531440UL) + ((uint64_t)op[4] * 496552594208016366UL) + ((uint64_t)op[5] * 11913621570039014168UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 11913621570039014168UL) + ((uint64_t)op[1] * 3713849246019306254UL) + ((uint64_t)op[2] * 12054912410056150093UL) + ((uint64_t)op[3] * 4743391980215843099UL) + ((((uint64_t)op[4] * 10151408009908531440UL) + ((uint64_t)op[5] * 496552594208016366UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 496552594208016366UL) + ((uint64_t)op[1] * 11913621570039014168UL) + ((uint64_t)op[2] * 3713849246019306254UL) + ((uint64_t)op[3] * 12054912410056150093UL) + ((uint64_t)op[4] * 4743391980215843099UL) + ((uint64_t)op[5] * 14734600181494529088UL);
	tmp_q[5] = ((uint64_t)op[0] * 10151408009908531440UL) + ((uint64_t)op[1] * 496552594208016366UL) + ((uint64_t)op[2] * 11913621570039014168UL) + ((uint64_t)op[3] * 3713849246019306254UL) + ((uint64_t)op[4] * 12054912410056150093UL) + ((uint64_t)op[5] * 4743391980215843099UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3160378482711L) - ((-((int128)tmp_q[1] * 1080772629975L) - ((int128)tmp_q[2] * 3271798142555L) - ((int128)tmp_q[3] * 388090873371L) - ((int128)tmp_q[4] * 3225245152169L) - ((int128)tmp_q[5] * 3136359513215L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 3136359513215L) - ((int128)tmp_q[1] * 3160378482711L) - ((-((int128)tmp_q[2] * 1080772629975L) - ((int128)tmp_q[3] * 3271798142555L) - ((int128)tmp_q[4] * 388090873371L) - ((int128)tmp_q[5] * 3225245152169L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 3225245152169L) - ((int128)tmp_q[1] * 3136359513215L) - ((int128)tmp_q[2] * 3160378482711L) - ((-((int128)tmp_q[3] * 1080772629975L) - ((int128)tmp_q[4] * 3271798142555L) - ((int128)tmp_q[5] * 388090873371L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 388090873371L) - ((int128)tmp_q[1] * 3225245152169L) - ((int128)tmp_q[2] * 3136359513215L) - ((int128)tmp_q[3] * 3160378482711L) - ((-((int128)tmp_q[4] * 1080772629975L) - ((int128)tmp_q[5] * 3271798142555L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 3271798142555L) - ((int128)tmp_q[1] * 388090873371L) - ((int128)tmp_q[2] * 3225245152169L) - ((int128)tmp_q[3] * 3136359513215L) - ((int128)tmp_q[4] * 3160378482711L) + ((int128)tmp_q[5] * 4323090519900L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1080772629975L) - ((int128)tmp_q[1] * 3271798142555L) - ((int128)tmp_q[2] * 388090873371L) - ((int128)tmp_q[3] * 3225245152169L) - ((int128)tmp_q[4] * 3136359513215L) - ((int128)tmp_q[5] * 3160378482711L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

