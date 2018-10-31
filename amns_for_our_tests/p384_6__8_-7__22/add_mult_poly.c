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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14251391383375371719UL) + ((((uint64_t)op[1] * 812040082186469399UL) + ((uint64_t)op[2] * 4991660447168437441UL) + ((uint64_t)op[3] * 12286438486169380410UL) + ((uint64_t)op[4] * 15981930321093064857UL) + ((uint64_t)op[5] * 16279053083065573703UL) + ((uint64_t)op[6] * 10192932819930287901UL) + ((uint64_t)op[7] * 13609877424277678863UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 13609877424277678863UL) + ((uint64_t)op[1] * 14251391383375371719UL) + ((((uint64_t)op[2] * 812040082186469399UL) + ((uint64_t)op[3] * 4991660447168437441UL) + ((uint64_t)op[4] * 12286438486169380410UL) + ((uint64_t)op[5] * 15981930321093064857UL) + ((uint64_t)op[6] * 16279053083065573703UL) + ((uint64_t)op[7] * 10192932819930287901UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 10192932819930287901UL) + ((uint64_t)op[1] * 13609877424277678863UL) + ((uint64_t)op[2] * 14251391383375371719UL) + ((((uint64_t)op[3] * 812040082186469399UL) + ((uint64_t)op[4] * 4991660447168437441UL) + ((uint64_t)op[5] * 12286438486169380410UL) + ((uint64_t)op[6] * 15981930321093064857UL) + ((uint64_t)op[7] * 16279053083065573703UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 16279053083065573703UL) + ((uint64_t)op[1] * 10192932819930287901UL) + ((uint64_t)op[2] * 13609877424277678863UL) + ((uint64_t)op[3] * 14251391383375371719UL) + ((((uint64_t)op[4] * 812040082186469399UL) + ((uint64_t)op[5] * 4991660447168437441UL) + ((uint64_t)op[6] * 12286438486169380410UL) + ((uint64_t)op[7] * 15981930321093064857UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 15981930321093064857UL) + ((uint64_t)op[1] * 16279053083065573703UL) + ((uint64_t)op[2] * 10192932819930287901UL) + ((uint64_t)op[3] * 13609877424277678863UL) + ((uint64_t)op[4] * 14251391383375371719UL) + ((((uint64_t)op[5] * 812040082186469399UL) + ((uint64_t)op[6] * 4991660447168437441UL) + ((uint64_t)op[7] * 12286438486169380410UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 12286438486169380410UL) + ((uint64_t)op[1] * 15981930321093064857UL) + ((uint64_t)op[2] * 16279053083065573703UL) + ((uint64_t)op[3] * 10192932819930287901UL) + ((uint64_t)op[4] * 13609877424277678863UL) + ((uint64_t)op[5] * 14251391383375371719UL) + ((((uint64_t)op[6] * 812040082186469399UL) + ((uint64_t)op[7] * 4991660447168437441UL)) * 18446744073709551609);
	tmp_q[6] = ((uint64_t)op[0] * 4991660447168437441UL) + ((uint64_t)op[1] * 12286438486169380410UL) + ((uint64_t)op[2] * 15981930321093064857UL) + ((uint64_t)op[3] * 16279053083065573703UL) + ((uint64_t)op[4] * 10192932819930287901UL) + ((uint64_t)op[5] * 13609877424277678863UL) + ((uint64_t)op[6] * 14251391383375371719UL) + ((uint64_t)op[7] * 12762463498404265823UL);
	tmp_q[7] = ((uint64_t)op[0] * 812040082186469399UL) + ((uint64_t)op[1] * 4991660447168437441UL) + ((uint64_t)op[2] * 12286438486169380410UL) + ((uint64_t)op[3] * 15981930321093064857UL) + ((uint64_t)op[4] * 16279053083065573703UL) + ((uint64_t)op[5] * 10192932819930287901UL) + ((uint64_t)op[6] * 13609877424277678863UL) + ((uint64_t)op[7] * 14251391383375371719UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 59111330024477L) - ((-((int128)tmp_q[1] * 162023228565281L) + ((int128)tmp_q[2] * 26956893739903L) - ((int128)tmp_q[3] * 92008790496257L) + ((int128)tmp_q[4] * 12367957105593L) + ((int128)tmp_q[5] * 47449261742860L) - ((int128)tmp_q[6] * 144867274125047L) + ((int128)tmp_q[7] * 115311016485619L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 115311016485619L) + ((int128)tmp_q[1] * 59111330024477L) - ((-((int128)tmp_q[2] * 162023228565281L) + ((int128)tmp_q[3] * 26956893739903L) - ((int128)tmp_q[4] * 92008790496257L) + ((int128)tmp_q[5] * 12367957105593L) + ((int128)tmp_q[6] * 47449261742860L) - ((int128)tmp_q[7] * 144867274125047L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 144867274125047L) + ((int128)tmp_q[1] * 115311016485619L) + ((int128)tmp_q[2] * 59111330024477L) - ((-((int128)tmp_q[3] * 162023228565281L) + ((int128)tmp_q[4] * 26956893739903L) - ((int128)tmp_q[5] * 92008790496257L) + ((int128)tmp_q[6] * 12367957105593L) + ((int128)tmp_q[7] * 47449261742860L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 47449261742860L) - ((int128)tmp_q[1] * 144867274125047L) + ((int128)tmp_q[2] * 115311016485619L) + ((int128)tmp_q[3] * 59111330024477L) - ((-((int128)tmp_q[4] * 162023228565281L) + ((int128)tmp_q[5] * 26956893739903L) - ((int128)tmp_q[6] * 92008790496257L) + ((int128)tmp_q[7] * 12367957105593L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 12367957105593L) + ((int128)tmp_q[1] * 47449261742860L) - ((int128)tmp_q[2] * 144867274125047L) + ((int128)tmp_q[3] * 115311016485619L) + ((int128)tmp_q[4] * 59111330024477L) - ((-((int128)tmp_q[5] * 162023228565281L) + ((int128)tmp_q[6] * 26956893739903L) - ((int128)tmp_q[7] * 92008790496257L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 92008790496257L) + ((int128)tmp_q[1] * 12367957105593L) + ((int128)tmp_q[2] * 47449261742860L) - ((int128)tmp_q[3] * 144867274125047L) + ((int128)tmp_q[4] * 115311016485619L) + ((int128)tmp_q[5] * 59111330024477L) - ((-((int128)tmp_q[6] * 162023228565281L) + ((int128)tmp_q[7] * 26956893739903L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 26956893739903L) - ((int128)tmp_q[1] * 92008790496257L) + ((int128)tmp_q[2] * 12367957105593L) + ((int128)tmp_q[3] * 47449261742860L) - ((int128)tmp_q[4] * 144867274125047L) + ((int128)tmp_q[5] * 115311016485619L) + ((int128)tmp_q[6] * 59111330024477L) + ((int128)tmp_q[7] * 1134162599956967L);
	tmp_zero[7] = -((int128)tmp_q[0] * 162023228565281L) + ((int128)tmp_q[1] * 26956893739903L) - ((int128)tmp_q[2] * 92008790496257L) + ((int128)tmp_q[3] * 12367957105593L) + ((int128)tmp_q[4] * 47449261742860L) - ((int128)tmp_q[5] * 144867274125047L) + ((int128)tmp_q[6] * 115311016485619L) + ((int128)tmp_q[7] * 59111330024477L);

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

