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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5616854483819164721UL) + ((((uint64_t)op[1] * 14371001463891355814UL) + ((uint64_t)op[2] * 17624769020852639252UL) + ((uint64_t)op[3] * 15924730722243145632UL) + ((uint64_t)op[4] * 16699243538379992068UL) + ((uint64_t)op[5] * 771100805996217020UL) + ((uint64_t)op[6] * 16072810767031056530UL) + ((uint64_t)op[7] * 13403188892413167857UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 13403188892413167857UL) + ((uint64_t)op[1] * 5616854483819164721UL) + ((((uint64_t)op[2] * 14371001463891355814UL) + ((uint64_t)op[3] * 17624769020852639252UL) + ((uint64_t)op[4] * 15924730722243145632UL) + ((uint64_t)op[5] * 16699243538379992068UL) + ((uint64_t)op[6] * 771100805996217020UL) + ((uint64_t)op[7] * 16072810767031056530UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 16072810767031056530UL) + ((uint64_t)op[1] * 13403188892413167857UL) + ((uint64_t)op[2] * 5616854483819164721UL) + ((((uint64_t)op[3] * 14371001463891355814UL) + ((uint64_t)op[4] * 17624769020852639252UL) + ((uint64_t)op[5] * 15924730722243145632UL) + ((uint64_t)op[6] * 16699243538379992068UL) + ((uint64_t)op[7] * 771100805996217020UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 771100805996217020UL) + ((uint64_t)op[1] * 16072810767031056530UL) + ((uint64_t)op[2] * 13403188892413167857UL) + ((uint64_t)op[3] * 5616854483819164721UL) + ((((uint64_t)op[4] * 14371001463891355814UL) + ((uint64_t)op[5] * 17624769020852639252UL) + ((uint64_t)op[6] * 15924730722243145632UL) + ((uint64_t)op[7] * 16699243538379992068UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 16699243538379992068UL) + ((uint64_t)op[1] * 771100805996217020UL) + ((uint64_t)op[2] * 16072810767031056530UL) + ((uint64_t)op[3] * 13403188892413167857UL) + ((uint64_t)op[4] * 5616854483819164721UL) + ((((uint64_t)op[5] * 14371001463891355814UL) + ((uint64_t)op[6] * 17624769020852639252UL) + ((uint64_t)op[7] * 15924730722243145632UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 15924730722243145632UL) + ((uint64_t)op[1] * 16699243538379992068UL) + ((uint64_t)op[2] * 771100805996217020UL) + ((uint64_t)op[3] * 16072810767031056530UL) + ((uint64_t)op[4] * 13403188892413167857UL) + ((uint64_t)op[5] * 5616854483819164721UL) + ((((uint64_t)op[6] * 14371001463891355814UL) + ((uint64_t)op[7] * 17624769020852639252UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 17624769020852639252UL) + ((uint64_t)op[1] * 15924730722243145632UL) + ((uint64_t)op[2] * 16699243538379992068UL) + ((uint64_t)op[3] * 771100805996217020UL) + ((uint64_t)op[4] * 16072810767031056530UL) + ((uint64_t)op[5] * 13403188892413167857UL) + ((uint64_t)op[6] * 5616854483819164721UL) + ((uint64_t)op[7] * 8151485219636391604UL);
	tmp_q[7] = ((uint64_t)op[0] * 14371001463891355814UL) + ((uint64_t)op[1] * 17624769020852639252UL) + ((uint64_t)op[2] * 15924730722243145632UL) + ((uint64_t)op[3] * 16699243538379992068UL) + ((uint64_t)op[4] * 771100805996217020UL) + ((uint64_t)op[5] * 16072810767031056530UL) + ((uint64_t)op[6] * 13403188892413167857UL) + ((uint64_t)op[7] * 5616854483819164721UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 41932947737201L) - ((-((int128)tmp_q[1] * 65582662583835L) + ((int128)tmp_q[2] * 9576573096827L) + ((int128)tmp_q[3] * 31505170155727L) + ((int128)tmp_q[4] * 28038124384511L) - ((int128)tmp_q[5] * 199096681272833L) + ((int128)tmp_q[6] * 10577831459947L) - ((int128)tmp_q[7] * 20610535805057L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 20610535805057L) + ((int128)tmp_q[1] * 41932947737201L) - ((-((int128)tmp_q[2] * 65582662583835L) + ((int128)tmp_q[3] * 9576573096827L) + ((int128)tmp_q[4] * 31505170155727L) + ((int128)tmp_q[5] * 28038124384511L) - ((int128)tmp_q[6] * 199096681272833L) + ((int128)tmp_q[7] * 10577831459947L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 10577831459947L) - ((int128)tmp_q[1] * 20610535805057L) + ((int128)tmp_q[2] * 41932947737201L) - ((-((int128)tmp_q[3] * 65582662583835L) + ((int128)tmp_q[4] * 9576573096827L) + ((int128)tmp_q[5] * 31505170155727L) + ((int128)tmp_q[6] * 28038124384511L) - ((int128)tmp_q[7] * 199096681272833L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 199096681272833L) + ((int128)tmp_q[1] * 10577831459947L) - ((int128)tmp_q[2] * 20610535805057L) + ((int128)tmp_q[3] * 41932947737201L) - ((-((int128)tmp_q[4] * 65582662583835L) + ((int128)tmp_q[5] * 9576573096827L) + ((int128)tmp_q[6] * 31505170155727L) + ((int128)tmp_q[7] * 28038124384511L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 28038124384511L) - ((int128)tmp_q[1] * 199096681272833L) + ((int128)tmp_q[2] * 10577831459947L) - ((int128)tmp_q[3] * 20610535805057L) + ((int128)tmp_q[4] * 41932947737201L) - ((-((int128)tmp_q[5] * 65582662583835L) + ((int128)tmp_q[6] * 9576573096827L) + ((int128)tmp_q[7] * 31505170155727L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 31505170155727L) + ((int128)tmp_q[1] * 28038124384511L) - ((int128)tmp_q[2] * 199096681272833L) + ((int128)tmp_q[3] * 10577831459947L) - ((int128)tmp_q[4] * 20610535805057L) + ((int128)tmp_q[5] * 41932947737201L) - ((-((int128)tmp_q[6] * 65582662583835L) + ((int128)tmp_q[7] * 9576573096827L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 9576573096827L) + ((int128)tmp_q[1] * 31505170155727L) + ((int128)tmp_q[2] * 28038124384511L) - ((int128)tmp_q[3] * 199096681272833L) + ((int128)tmp_q[4] * 10577831459947L) - ((int128)tmp_q[5] * 20610535805057L) + ((int128)tmp_q[6] * 41932947737201L) + ((int128)tmp_q[7] * 131165325167670L);
	tmp_zero[7] = -((int128)tmp_q[0] * 65582662583835L) + ((int128)tmp_q[1] * 9576573096827L) + ((int128)tmp_q[2] * 31505170155727L) + ((int128)tmp_q[3] * 28038124384511L) - ((int128)tmp_q[4] * 199096681272833L) + ((int128)tmp_q[5] * 10577831459947L) - ((int128)tmp_q[6] * 20610535805057L) + ((int128)tmp_q[7] * 41932947737201L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

