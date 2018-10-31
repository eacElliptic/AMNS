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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16660642165031265267UL) + ((((uint64_t)op[1] * 12261502418162021185UL) + ((uint64_t)op[2] * 7931594085772667453UL) + ((uint64_t)op[3] * 14411819141517247957UL) + ((uint64_t)op[4] * 16177109392092187347UL) + ((uint64_t)op[5] * 7019050322539760260UL) + ((uint64_t)op[6] * 12673696007589134148UL) + ((uint64_t)op[7] * 13320741179269790914UL) + ((uint64_t)op[8] * 16019455570029768599UL) + ((uint64_t)op[9] * 949490978269605801UL) + ((uint64_t)op[10] * 18427021819100861506UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 18427021819100861506UL) + ((uint64_t)op[1] * 16660642165031265267UL) + ((((uint64_t)op[2] * 12261502418162021185UL) + ((uint64_t)op[3] * 7931594085772667453UL) + ((uint64_t)op[4] * 14411819141517247957UL) + ((uint64_t)op[5] * 16177109392092187347UL) + ((uint64_t)op[6] * 7019050322539760260UL) + ((uint64_t)op[7] * 12673696007589134148UL) + ((uint64_t)op[8] * 13320741179269790914UL) + ((uint64_t)op[9] * 16019455570029768599UL) + ((uint64_t)op[10] * 949490978269605801UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 949490978269605801UL) + ((uint64_t)op[1] * 18427021819100861506UL) + ((uint64_t)op[2] * 16660642165031265267UL) + ((((uint64_t)op[3] * 12261502418162021185UL) + ((uint64_t)op[4] * 7931594085772667453UL) + ((uint64_t)op[5] * 14411819141517247957UL) + ((uint64_t)op[6] * 16177109392092187347UL) + ((uint64_t)op[7] * 7019050322539760260UL) + ((uint64_t)op[8] * 12673696007589134148UL) + ((uint64_t)op[9] * 13320741179269790914UL) + ((uint64_t)op[10] * 16019455570029768599UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 16019455570029768599UL) + ((uint64_t)op[1] * 949490978269605801UL) + ((uint64_t)op[2] * 18427021819100861506UL) + ((uint64_t)op[3] * 16660642165031265267UL) + ((((uint64_t)op[4] * 12261502418162021185UL) + ((uint64_t)op[5] * 7931594085772667453UL) + ((uint64_t)op[6] * 14411819141517247957UL) + ((uint64_t)op[7] * 16177109392092187347UL) + ((uint64_t)op[8] * 7019050322539760260UL) + ((uint64_t)op[9] * 12673696007589134148UL) + ((uint64_t)op[10] * 13320741179269790914UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 13320741179269790914UL) + ((uint64_t)op[1] * 16019455570029768599UL) + ((uint64_t)op[2] * 949490978269605801UL) + ((uint64_t)op[3] * 18427021819100861506UL) + ((uint64_t)op[4] * 16660642165031265267UL) + ((((uint64_t)op[5] * 12261502418162021185UL) + ((uint64_t)op[6] * 7931594085772667453UL) + ((uint64_t)op[7] * 14411819141517247957UL) + ((uint64_t)op[8] * 16177109392092187347UL) + ((uint64_t)op[9] * 7019050322539760260UL) + ((uint64_t)op[10] * 12673696007589134148UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 12673696007589134148UL) + ((uint64_t)op[1] * 13320741179269790914UL) + ((uint64_t)op[2] * 16019455570029768599UL) + ((uint64_t)op[3] * 949490978269605801UL) + ((uint64_t)op[4] * 18427021819100861506UL) + ((uint64_t)op[5] * 16660642165031265267UL) + ((((uint64_t)op[6] * 12261502418162021185UL) + ((uint64_t)op[7] * 7931594085772667453UL) + ((uint64_t)op[8] * 14411819141517247957UL) + ((uint64_t)op[9] * 16177109392092187347UL) + ((uint64_t)op[10] * 7019050322539760260UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 7019050322539760260UL) + ((uint64_t)op[1] * 12673696007589134148UL) + ((uint64_t)op[2] * 13320741179269790914UL) + ((uint64_t)op[3] * 16019455570029768599UL) + ((uint64_t)op[4] * 949490978269605801UL) + ((uint64_t)op[5] * 18427021819100861506UL) + ((uint64_t)op[6] * 16660642165031265267UL) + ((((uint64_t)op[7] * 12261502418162021185UL) + ((uint64_t)op[8] * 7931594085772667453UL) + ((uint64_t)op[9] * 14411819141517247957UL) + ((uint64_t)op[10] * 16177109392092187347UL)) * 18446744073709551608);
	tmp_q[7] = ((uint64_t)op[0] * 16177109392092187347UL) + ((uint64_t)op[1] * 7019050322539760260UL) + ((uint64_t)op[2] * 12673696007589134148UL) + ((uint64_t)op[3] * 13320741179269790914UL) + ((uint64_t)op[4] * 16019455570029768599UL) + ((uint64_t)op[5] * 949490978269605801UL) + ((uint64_t)op[6] * 18427021819100861506UL) + ((uint64_t)op[7] * 16660642165031265267UL) + ((((uint64_t)op[8] * 12261502418162021185UL) + ((uint64_t)op[9] * 7931594085772667453UL) + ((uint64_t)op[10] * 14411819141517247957UL)) * 18446744073709551608);
	tmp_q[8] = ((uint64_t)op[0] * 14411819141517247957UL) + ((uint64_t)op[1] * 16177109392092187347UL) + ((uint64_t)op[2] * 7019050322539760260UL) + ((uint64_t)op[3] * 12673696007589134148UL) + ((uint64_t)op[4] * 13320741179269790914UL) + ((uint64_t)op[5] * 16019455570029768599UL) + ((uint64_t)op[6] * 949490978269605801UL) + ((uint64_t)op[7] * 18427021819100861506UL) + ((uint64_t)op[8] * 16660642165031265267UL) + ((((uint64_t)op[9] * 12261502418162021185UL) + ((uint64_t)op[10] * 7931594085772667453UL)) * 18446744073709551608);
	tmp_q[9] = ((uint64_t)op[0] * 7931594085772667453UL) + ((uint64_t)op[1] * 14411819141517247957UL) + ((uint64_t)op[2] * 16177109392092187347UL) + ((uint64_t)op[3] * 7019050322539760260UL) + ((uint64_t)op[4] * 12673696007589134148UL) + ((uint64_t)op[5] * 13320741179269790914UL) + ((uint64_t)op[6] * 16019455570029768599UL) + ((uint64_t)op[7] * 949490978269605801UL) + ((uint64_t)op[8] * 18427021819100861506UL) + ((uint64_t)op[9] * 16660642165031265267UL) + ((uint64_t)op[10] * 12588445096961140216UL);
	tmp_q[10] = ((uint64_t)op[0] * 12261502418162021185UL) + ((uint64_t)op[1] * 7931594085772667453UL) + ((uint64_t)op[2] * 14411819141517247957UL) + ((uint64_t)op[3] * 16177109392092187347UL) + ((uint64_t)op[4] * 7019050322539760260UL) + ((uint64_t)op[5] * 12673696007589134148UL) + ((uint64_t)op[6] * 13320741179269790914UL) + ((uint64_t)op[7] * 16019455570029768599UL) + ((uint64_t)op[8] * 949490978269605801UL) + ((uint64_t)op[9] * 18427021819100861506UL) + ((uint64_t)op[10] * 16660642165031265267UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 27246888244635L) - ((((int128)tmp_q[1] * 4230847855138L) - ((int128)tmp_q[2] * 78919496091016L) - ((int128)tmp_q[3] * 40931149472149L) - ((int128)tmp_q[4] * 57560049703766L) + ((int128)tmp_q[5] * 7643911035690L) + ((int128)tmp_q[6] * 16715813767668L) - ((int128)tmp_q[7] * 90156090999721L) - ((int128)tmp_q[8] * 52445005164181L) + ((int128)tmp_q[9] * 31356111334101L) - ((int128)tmp_q[10] * 104587576292422L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 104587576292422L) - ((int128)tmp_q[1] * 27246888244635L) - ((((int128)tmp_q[2] * 4230847855138L) - ((int128)tmp_q[3] * 78919496091016L) - ((int128)tmp_q[4] * 40931149472149L) - ((int128)tmp_q[5] * 57560049703766L) + ((int128)tmp_q[6] * 7643911035690L) + ((int128)tmp_q[7] * 16715813767668L) - ((int128)tmp_q[8] * 90156090999721L) - ((int128)tmp_q[9] * 52445005164181L) + ((int128)tmp_q[10] * 31356111334101L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 31356111334101L) - ((int128)tmp_q[1] * 104587576292422L) - ((int128)tmp_q[2] * 27246888244635L) - ((((int128)tmp_q[3] * 4230847855138L) - ((int128)tmp_q[4] * 78919496091016L) - ((int128)tmp_q[5] * 40931149472149L) - ((int128)tmp_q[6] * 57560049703766L) + ((int128)tmp_q[7] * 7643911035690L) + ((int128)tmp_q[8] * 16715813767668L) - ((int128)tmp_q[9] * 90156090999721L) - ((int128)tmp_q[10] * 52445005164181L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 52445005164181L) + ((int128)tmp_q[1] * 31356111334101L) - ((int128)tmp_q[2] * 104587576292422L) - ((int128)tmp_q[3] * 27246888244635L) - ((((int128)tmp_q[4] * 4230847855138L) - ((int128)tmp_q[5] * 78919496091016L) - ((int128)tmp_q[6] * 40931149472149L) - ((int128)tmp_q[7] * 57560049703766L) + ((int128)tmp_q[8] * 7643911035690L) + ((int128)tmp_q[9] * 16715813767668L) - ((int128)tmp_q[10] * 90156090999721L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 90156090999721L) - ((int128)tmp_q[1] * 52445005164181L) + ((int128)tmp_q[2] * 31356111334101L) - ((int128)tmp_q[3] * 104587576292422L) - ((int128)tmp_q[4] * 27246888244635L) - ((((int128)tmp_q[5] * 4230847855138L) - ((int128)tmp_q[6] * 78919496091016L) - ((int128)tmp_q[7] * 40931149472149L) - ((int128)tmp_q[8] * 57560049703766L) + ((int128)tmp_q[9] * 7643911035690L) + ((int128)tmp_q[10] * 16715813767668L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 16715813767668L) - ((int128)tmp_q[1] * 90156090999721L) - ((int128)tmp_q[2] * 52445005164181L) + ((int128)tmp_q[3] * 31356111334101L) - ((int128)tmp_q[4] * 104587576292422L) - ((int128)tmp_q[5] * 27246888244635L) - ((((int128)tmp_q[6] * 4230847855138L) - ((int128)tmp_q[7] * 78919496091016L) - ((int128)tmp_q[8] * 40931149472149L) - ((int128)tmp_q[9] * 57560049703766L) + ((int128)tmp_q[10] * 7643911035690L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 7643911035690L) + ((int128)tmp_q[1] * 16715813767668L) - ((int128)tmp_q[2] * 90156090999721L) - ((int128)tmp_q[3] * 52445005164181L) + ((int128)tmp_q[4] * 31356111334101L) - ((int128)tmp_q[5] * 104587576292422L) - ((int128)tmp_q[6] * 27246888244635L) - ((((int128)tmp_q[7] * 4230847855138L) - ((int128)tmp_q[8] * 78919496091016L) - ((int128)tmp_q[9] * 40931149472149L) - ((int128)tmp_q[10] * 57560049703766L)) * 8);
	tmp_zero[7] = -((int128)tmp_q[0] * 57560049703766L) + ((int128)tmp_q[1] * 7643911035690L) + ((int128)tmp_q[2] * 16715813767668L) - ((int128)tmp_q[3] * 90156090999721L) - ((int128)tmp_q[4] * 52445005164181L) + ((int128)tmp_q[5] * 31356111334101L) - ((int128)tmp_q[6] * 104587576292422L) - ((int128)tmp_q[7] * 27246888244635L) - ((((int128)tmp_q[8] * 4230847855138L) - ((int128)tmp_q[9] * 78919496091016L) - ((int128)tmp_q[10] * 40931149472149L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 40931149472149L) - ((int128)tmp_q[1] * 57560049703766L) + ((int128)tmp_q[2] * 7643911035690L) + ((int128)tmp_q[3] * 16715813767668L) - ((int128)tmp_q[4] * 90156090999721L) - ((int128)tmp_q[5] * 52445005164181L) + ((int128)tmp_q[6] * 31356111334101L) - ((int128)tmp_q[7] * 104587576292422L) - ((int128)tmp_q[8] * 27246888244635L) - ((((int128)tmp_q[9] * 4230847855138L) - ((int128)tmp_q[10] * 78919496091016L)) * 8);
	tmp_zero[9] = -((int128)tmp_q[0] * 78919496091016L) - ((int128)tmp_q[1] * 40931149472149L) - ((int128)tmp_q[2] * 57560049703766L) + ((int128)tmp_q[3] * 7643911035690L) + ((int128)tmp_q[4] * 16715813767668L) - ((int128)tmp_q[5] * 90156090999721L) - ((int128)tmp_q[6] * 52445005164181L) + ((int128)tmp_q[7] * 31356111334101L) - ((int128)tmp_q[8] * 104587576292422L) - ((int128)tmp_q[9] * 27246888244635L) - ((int128)tmp_q[10] * 33846782841104L);
	tmp_zero[10] = ((int128)tmp_q[0] * 4230847855138L) - ((int128)tmp_q[1] * 78919496091016L) - ((int128)tmp_q[2] * 40931149472149L) - ((int128)tmp_q[3] * 57560049703766L) + ((int128)tmp_q[4] * 7643911035690L) + ((int128)tmp_q[5] * 16715813767668L) - ((int128)tmp_q[6] * 90156090999721L) - ((int128)tmp_q[7] * 52445005164181L) + ((int128)tmp_q[8] * 31356111334101L) - ((int128)tmp_q[9] * 104587576292422L) - ((int128)tmp_q[10] * 27246888244635L);

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
}

