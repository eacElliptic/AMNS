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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) << 1);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) << 1);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12244769951407379845UL) + ((((uint64_t)op[1] * 6576809034633757256UL) + ((uint64_t)op[2] * 17511026182575434339UL) + ((uint64_t)op[3] * 13190968276066342230UL) + ((uint64_t)op[4] * 3895805748254427670UL) + ((uint64_t)op[5] * 7659444158275050714UL) + ((uint64_t)op[6] * 12882217396591551391UL) + ((uint64_t)op[7] * 6528724723375687858UL) + ((uint64_t)op[8] * 1582981382933837901UL) + ((uint64_t)op[9] * 10131860461472180963UL) + ((uint64_t)op[10] * 14847956638722417591UL) + ((uint64_t)op[11] * 3453704736416508062UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 3453704736416508062UL) + ((uint64_t)op[1] * 12244769951407379845UL) + ((((uint64_t)op[2] * 6576809034633757256UL) + ((uint64_t)op[3] * 17511026182575434339UL) + ((uint64_t)op[4] * 13190968276066342230UL) + ((uint64_t)op[5] * 3895805748254427670UL) + ((uint64_t)op[6] * 7659444158275050714UL) + ((uint64_t)op[7] * 12882217396591551391UL) + ((uint64_t)op[8] * 6528724723375687858UL) + ((uint64_t)op[9] * 1582981382933837901UL) + ((uint64_t)op[10] * 10131860461472180963UL) + ((uint64_t)op[11] * 14847956638722417591UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 14847956638722417591UL) + ((uint64_t)op[1] * 3453704736416508062UL) + ((uint64_t)op[2] * 12244769951407379845UL) + ((((uint64_t)op[3] * 6576809034633757256UL) + ((uint64_t)op[4] * 17511026182575434339UL) + ((uint64_t)op[5] * 13190968276066342230UL) + ((uint64_t)op[6] * 3895805748254427670UL) + ((uint64_t)op[7] * 7659444158275050714UL) + ((uint64_t)op[8] * 12882217396591551391UL) + ((uint64_t)op[9] * 6528724723375687858UL) + ((uint64_t)op[10] * 1582981382933837901UL) + ((uint64_t)op[11] * 10131860461472180963UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 10131860461472180963UL) + ((uint64_t)op[1] * 14847956638722417591UL) + ((uint64_t)op[2] * 3453704736416508062UL) + ((uint64_t)op[3] * 12244769951407379845UL) + ((((uint64_t)op[4] * 6576809034633757256UL) + ((uint64_t)op[5] * 17511026182575434339UL) + ((uint64_t)op[6] * 13190968276066342230UL) + ((uint64_t)op[7] * 3895805748254427670UL) + ((uint64_t)op[8] * 7659444158275050714UL) + ((uint64_t)op[9] * 12882217396591551391UL) + ((uint64_t)op[10] * 6528724723375687858UL) + ((uint64_t)op[11] * 1582981382933837901UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 1582981382933837901UL) + ((uint64_t)op[1] * 10131860461472180963UL) + ((uint64_t)op[2] * 14847956638722417591UL) + ((uint64_t)op[3] * 3453704736416508062UL) + ((uint64_t)op[4] * 12244769951407379845UL) + ((((uint64_t)op[5] * 6576809034633757256UL) + ((uint64_t)op[6] * 17511026182575434339UL) + ((uint64_t)op[7] * 13190968276066342230UL) + ((uint64_t)op[8] * 3895805748254427670UL) + ((uint64_t)op[9] * 7659444158275050714UL) + ((uint64_t)op[10] * 12882217396591551391UL) + ((uint64_t)op[11] * 6528724723375687858UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 6528724723375687858UL) + ((uint64_t)op[1] * 1582981382933837901UL) + ((uint64_t)op[2] * 10131860461472180963UL) + ((uint64_t)op[3] * 14847956638722417591UL) + ((uint64_t)op[4] * 3453704736416508062UL) + ((uint64_t)op[5] * 12244769951407379845UL) + ((((uint64_t)op[6] * 6576809034633757256UL) + ((uint64_t)op[7] * 17511026182575434339UL) + ((uint64_t)op[8] * 13190968276066342230UL) + ((uint64_t)op[9] * 3895805748254427670UL) + ((uint64_t)op[10] * 7659444158275050714UL) + ((uint64_t)op[11] * 12882217396591551391UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 12882217396591551391UL) + ((uint64_t)op[1] * 6528724723375687858UL) + ((uint64_t)op[2] * 1582981382933837901UL) + ((uint64_t)op[3] * 10131860461472180963UL) + ((uint64_t)op[4] * 14847956638722417591UL) + ((uint64_t)op[5] * 3453704736416508062UL) + ((uint64_t)op[6] * 12244769951407379845UL) + ((((uint64_t)op[7] * 6576809034633757256UL) + ((uint64_t)op[8] * 17511026182575434339UL) + ((uint64_t)op[9] * 13190968276066342230UL) + ((uint64_t)op[10] * 3895805748254427670UL) + ((uint64_t)op[11] * 7659444158275050714UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 7659444158275050714UL) + ((uint64_t)op[1] * 12882217396591551391UL) + ((uint64_t)op[2] * 6528724723375687858UL) + ((uint64_t)op[3] * 1582981382933837901UL) + ((uint64_t)op[4] * 10131860461472180963UL) + ((uint64_t)op[5] * 14847956638722417591UL) + ((uint64_t)op[6] * 3453704736416508062UL) + ((uint64_t)op[7] * 12244769951407379845UL) + ((((uint64_t)op[8] * 6576809034633757256UL) + ((uint64_t)op[9] * 17511026182575434339UL) + ((uint64_t)op[10] * 13190968276066342230UL) + ((uint64_t)op[11] * 3895805748254427670UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 3895805748254427670UL) + ((uint64_t)op[1] * 7659444158275050714UL) + ((uint64_t)op[2] * 12882217396591551391UL) + ((uint64_t)op[3] * 6528724723375687858UL) + ((uint64_t)op[4] * 1582981382933837901UL) + ((uint64_t)op[5] * 10131860461472180963UL) + ((uint64_t)op[6] * 14847956638722417591UL) + ((uint64_t)op[7] * 3453704736416508062UL) + ((uint64_t)op[8] * 12244769951407379845UL) + ((((uint64_t)op[9] * 6576809034633757256UL) + ((uint64_t)op[10] * 17511026182575434339UL) + ((uint64_t)op[11] * 13190968276066342230UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 13190968276066342230UL) + ((uint64_t)op[1] * 3895805748254427670UL) + ((uint64_t)op[2] * 7659444158275050714UL) + ((uint64_t)op[3] * 12882217396591551391UL) + ((uint64_t)op[4] * 6528724723375687858UL) + ((uint64_t)op[5] * 1582981382933837901UL) + ((uint64_t)op[6] * 10131860461472180963UL) + ((uint64_t)op[7] * 14847956638722417591UL) + ((uint64_t)op[8] * 3453704736416508062UL) + ((uint64_t)op[9] * 12244769951407379845UL) + ((((uint64_t)op[10] * 6576809034633757256UL) + ((uint64_t)op[11] * 17511026182575434339UL)) * 2);
	tmp_q[10] = ((uint64_t)op[0] * 17511026182575434339UL) + ((uint64_t)op[1] * 13190968276066342230UL) + ((uint64_t)op[2] * 3895805748254427670UL) + ((uint64_t)op[3] * 7659444158275050714UL) + ((uint64_t)op[4] * 12882217396591551391UL) + ((uint64_t)op[5] * 6528724723375687858UL) + ((uint64_t)op[6] * 1582981382933837901UL) + ((uint64_t)op[7] * 10131860461472180963UL) + ((uint64_t)op[8] * 14847956638722417591UL) + ((uint64_t)op[9] * 3453704736416508062UL) + ((uint64_t)op[10] * 12244769951407379845UL) + ((uint64_t)op[11] * 13153618069267514512UL);
	tmp_q[11] = ((uint64_t)op[0] * 6576809034633757256UL) + ((uint64_t)op[1] * 17511026182575434339UL) + ((uint64_t)op[2] * 13190968276066342230UL) + ((uint64_t)op[3] * 3895805748254427670UL) + ((uint64_t)op[4] * 7659444158275050714UL) + ((uint64_t)op[5] * 12882217396591551391UL) + ((uint64_t)op[6] * 6528724723375687858UL) + ((uint64_t)op[7] * 1582981382933837901UL) + ((uint64_t)op[8] * 10131860461472180963UL) + ((uint64_t)op[9] * 14847956638722417591UL) + ((uint64_t)op[10] * 3453704736416508062UL) + ((uint64_t)op[11] * 12244769951407379845UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2654914074067L) + ((((int128)tmp_q[1] * 488147865300L) + ((int128)tmp_q[2] * 1621270044345L) + ((int128)tmp_q[3] * 5586838634425L) + ((int128)tmp_q[4] * 2629211345828L) + ((int128)tmp_q[5] * 4710521375939L) + ((int128)tmp_q[6] * 6708231327849L) + ((int128)tmp_q[7] * 1193847570444L) - ((int128)tmp_q[8] * 1991633098806L) + ((int128)tmp_q[9] * 4835611276965L) - ((int128)tmp_q[10] * 2911267132887L) + ((int128)tmp_q[11] * 2791698083058L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 2791698083058L) + ((int128)tmp_q[1] * 2654914074067L) + ((((int128)tmp_q[2] * 488147865300L) + ((int128)tmp_q[3] * 1621270044345L) + ((int128)tmp_q[4] * 5586838634425L) + ((int128)tmp_q[5] * 2629211345828L) + ((int128)tmp_q[6] * 4710521375939L) + ((int128)tmp_q[7] * 6708231327849L) + ((int128)tmp_q[8] * 1193847570444L) - ((int128)tmp_q[9] * 1991633098806L) + ((int128)tmp_q[10] * 4835611276965L) - ((int128)tmp_q[11] * 2911267132887L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 2911267132887L) + ((int128)tmp_q[1] * 2791698083058L) + ((int128)tmp_q[2] * 2654914074067L) + ((((int128)tmp_q[3] * 488147865300L) + ((int128)tmp_q[4] * 1621270044345L) + ((int128)tmp_q[5] * 5586838634425L) + ((int128)tmp_q[6] * 2629211345828L) + ((int128)tmp_q[7] * 4710521375939L) + ((int128)tmp_q[8] * 6708231327849L) + ((int128)tmp_q[9] * 1193847570444L) - ((int128)tmp_q[10] * 1991633098806L) + ((int128)tmp_q[11] * 4835611276965L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 4835611276965L) - ((int128)tmp_q[1] * 2911267132887L) + ((int128)tmp_q[2] * 2791698083058L) + ((int128)tmp_q[3] * 2654914074067L) + ((((int128)tmp_q[4] * 488147865300L) + ((int128)tmp_q[5] * 1621270044345L) + ((int128)tmp_q[6] * 5586838634425L) + ((int128)tmp_q[7] * 2629211345828L) + ((int128)tmp_q[8] * 4710521375939L) + ((int128)tmp_q[9] * 6708231327849L) + ((int128)tmp_q[10] * 1193847570444L) - ((int128)tmp_q[11] * 1991633098806L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1991633098806L) + ((int128)tmp_q[1] * 4835611276965L) - ((int128)tmp_q[2] * 2911267132887L) + ((int128)tmp_q[3] * 2791698083058L) + ((int128)tmp_q[4] * 2654914074067L) + ((((int128)tmp_q[5] * 488147865300L) + ((int128)tmp_q[6] * 1621270044345L) + ((int128)tmp_q[7] * 5586838634425L) + ((int128)tmp_q[8] * 2629211345828L) + ((int128)tmp_q[9] * 4710521375939L) + ((int128)tmp_q[10] * 6708231327849L) + ((int128)tmp_q[11] * 1193847570444L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 1193847570444L) - ((int128)tmp_q[1] * 1991633098806L) + ((int128)tmp_q[2] * 4835611276965L) - ((int128)tmp_q[3] * 2911267132887L) + ((int128)tmp_q[4] * 2791698083058L) + ((int128)tmp_q[5] * 2654914074067L) + ((((int128)tmp_q[6] * 488147865300L) + ((int128)tmp_q[7] * 1621270044345L) + ((int128)tmp_q[8] * 5586838634425L) + ((int128)tmp_q[9] * 2629211345828L) + ((int128)tmp_q[10] * 4710521375939L) + ((int128)tmp_q[11] * 6708231327849L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 6708231327849L) + ((int128)tmp_q[1] * 1193847570444L) - ((int128)tmp_q[2] * 1991633098806L) + ((int128)tmp_q[3] * 4835611276965L) - ((int128)tmp_q[4] * 2911267132887L) + ((int128)tmp_q[5] * 2791698083058L) + ((int128)tmp_q[6] * 2654914074067L) + ((((int128)tmp_q[7] * 488147865300L) + ((int128)tmp_q[8] * 1621270044345L) + ((int128)tmp_q[9] * 5586838634425L) + ((int128)tmp_q[10] * 2629211345828L) + ((int128)tmp_q[11] * 4710521375939L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 4710521375939L) + ((int128)tmp_q[1] * 6708231327849L) + ((int128)tmp_q[2] * 1193847570444L) - ((int128)tmp_q[3] * 1991633098806L) + ((int128)tmp_q[4] * 4835611276965L) - ((int128)tmp_q[5] * 2911267132887L) + ((int128)tmp_q[6] * 2791698083058L) + ((int128)tmp_q[7] * 2654914074067L) + ((((int128)tmp_q[8] * 488147865300L) + ((int128)tmp_q[9] * 1621270044345L) + ((int128)tmp_q[10] * 5586838634425L) + ((int128)tmp_q[11] * 2629211345828L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 2629211345828L) + ((int128)tmp_q[1] * 4710521375939L) + ((int128)tmp_q[2] * 6708231327849L) + ((int128)tmp_q[3] * 1193847570444L) - ((int128)tmp_q[4] * 1991633098806L) + ((int128)tmp_q[5] * 4835611276965L) - ((int128)tmp_q[6] * 2911267132887L) + ((int128)tmp_q[7] * 2791698083058L) + ((int128)tmp_q[8] * 2654914074067L) + ((((int128)tmp_q[9] * 488147865300L) + ((int128)tmp_q[10] * 1621270044345L) + ((int128)tmp_q[11] * 5586838634425L)) * 2);
	tmp_zero[9] = ((int128)tmp_q[0] * 5586838634425L) + ((int128)tmp_q[1] * 2629211345828L) + ((int128)tmp_q[2] * 4710521375939L) + ((int128)tmp_q[3] * 6708231327849L) + ((int128)tmp_q[4] * 1193847570444L) - ((int128)tmp_q[5] * 1991633098806L) + ((int128)tmp_q[6] * 4835611276965L) - ((int128)tmp_q[7] * 2911267132887L) + ((int128)tmp_q[8] * 2791698083058L) + ((int128)tmp_q[9] * 2654914074067L) + ((((int128)tmp_q[10] * 488147865300L) + ((int128)tmp_q[11] * 1621270044345L)) * 2);
	tmp_zero[10] = ((int128)tmp_q[0] * 1621270044345L) + ((int128)tmp_q[1] * 5586838634425L) + ((int128)tmp_q[2] * 2629211345828L) + ((int128)tmp_q[3] * 4710521375939L) + ((int128)tmp_q[4] * 6708231327849L) + ((int128)tmp_q[5] * 1193847570444L) - ((int128)tmp_q[6] * 1991633098806L) + ((int128)tmp_q[7] * 4835611276965L) - ((int128)tmp_q[8] * 2911267132887L) + ((int128)tmp_q[9] * 2791698083058L) + ((int128)tmp_q[10] * 2654914074067L) + ((int128)tmp_q[11] * 976295730600L);
	tmp_zero[11] = ((int128)tmp_q[0] * 488147865300L) + ((int128)tmp_q[1] * 1621270044345L) + ((int128)tmp_q[2] * 5586838634425L) + ((int128)tmp_q[3] * 2629211345828L) + ((int128)tmp_q[4] * 4710521375939L) + ((int128)tmp_q[5] * 6708231327849L) + ((int128)tmp_q[6] * 1193847570444L) - ((int128)tmp_q[7] * 1991633098806L) + ((int128)tmp_q[8] * 4835611276965L) - ((int128)tmp_q[9] * 2911267132887L) + ((int128)tmp_q[10] * 2791698083058L) + ((int128)tmp_q[11] * 2654914074067L);

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
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

