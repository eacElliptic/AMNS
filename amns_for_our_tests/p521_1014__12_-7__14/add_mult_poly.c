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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 7);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 14);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 14);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 7);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3155881917008880357UL) + ((((uint64_t)op[1] * 2502975817053249358UL) + ((uint64_t)op[2] * 2783709102958595942UL) + ((uint64_t)op[3] * 4303861093280866496UL) + ((uint64_t)op[4] * 12442483447733475858UL) + ((uint64_t)op[5] * 17574392635564835732UL) + ((uint64_t)op[6] * 13618301421045475485UL) + ((uint64_t)op[7] * 11562123458551789729UL) + ((uint64_t)op[8] * 16387303807734353009UL) + ((uint64_t)op[9] * 2443114270154127711UL) + ((uint64_t)op[10] * 18403820962515117585UL) + ((uint64_t)op[11] * 17674925139268631843UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 17674925139268631843UL) + ((uint64_t)op[1] * 3155881917008880357UL) + ((((uint64_t)op[2] * 2502975817053249358UL) + ((uint64_t)op[3] * 2783709102958595942UL) + ((uint64_t)op[4] * 4303861093280866496UL) + ((uint64_t)op[5] * 12442483447733475858UL) + ((uint64_t)op[6] * 17574392635564835732UL) + ((uint64_t)op[7] * 13618301421045475485UL) + ((uint64_t)op[8] * 11562123458551789729UL) + ((uint64_t)op[9] * 16387303807734353009UL) + ((uint64_t)op[10] * 2443114270154127711UL) + ((uint64_t)op[11] * 18403820962515117585UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 18403820962515117585UL) + ((uint64_t)op[1] * 17674925139268631843UL) + ((uint64_t)op[2] * 3155881917008880357UL) + ((((uint64_t)op[3] * 2502975817053249358UL) + ((uint64_t)op[4] * 2783709102958595942UL) + ((uint64_t)op[5] * 4303861093280866496UL) + ((uint64_t)op[6] * 12442483447733475858UL) + ((uint64_t)op[7] * 17574392635564835732UL) + ((uint64_t)op[8] * 13618301421045475485UL) + ((uint64_t)op[9] * 11562123458551789729UL) + ((uint64_t)op[10] * 16387303807734353009UL) + ((uint64_t)op[11] * 2443114270154127711UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 2443114270154127711UL) + ((uint64_t)op[1] * 18403820962515117585UL) + ((uint64_t)op[2] * 17674925139268631843UL) + ((uint64_t)op[3] * 3155881917008880357UL) + ((((uint64_t)op[4] * 2502975817053249358UL) + ((uint64_t)op[5] * 2783709102958595942UL) + ((uint64_t)op[6] * 4303861093280866496UL) + ((uint64_t)op[7] * 12442483447733475858UL) + ((uint64_t)op[8] * 17574392635564835732UL) + ((uint64_t)op[9] * 13618301421045475485UL) + ((uint64_t)op[10] * 11562123458551789729UL) + ((uint64_t)op[11] * 16387303807734353009UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 16387303807734353009UL) + ((uint64_t)op[1] * 2443114270154127711UL) + ((uint64_t)op[2] * 18403820962515117585UL) + ((uint64_t)op[3] * 17674925139268631843UL) + ((uint64_t)op[4] * 3155881917008880357UL) + ((((uint64_t)op[5] * 2502975817053249358UL) + ((uint64_t)op[6] * 2783709102958595942UL) + ((uint64_t)op[7] * 4303861093280866496UL) + ((uint64_t)op[8] * 12442483447733475858UL) + ((uint64_t)op[9] * 17574392635564835732UL) + ((uint64_t)op[10] * 13618301421045475485UL) + ((uint64_t)op[11] * 11562123458551789729UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 11562123458551789729UL) + ((uint64_t)op[1] * 16387303807734353009UL) + ((uint64_t)op[2] * 2443114270154127711UL) + ((uint64_t)op[3] * 18403820962515117585UL) + ((uint64_t)op[4] * 17674925139268631843UL) + ((uint64_t)op[5] * 3155881917008880357UL) + ((((uint64_t)op[6] * 2502975817053249358UL) + ((uint64_t)op[7] * 2783709102958595942UL) + ((uint64_t)op[8] * 4303861093280866496UL) + ((uint64_t)op[9] * 12442483447733475858UL) + ((uint64_t)op[10] * 17574392635564835732UL) + ((uint64_t)op[11] * 13618301421045475485UL)) * 18446744073709551609);
	tmp_q[6] = ((uint64_t)op[0] * 13618301421045475485UL) + ((uint64_t)op[1] * 11562123458551789729UL) + ((uint64_t)op[2] * 16387303807734353009UL) + ((uint64_t)op[3] * 2443114270154127711UL) + ((uint64_t)op[4] * 18403820962515117585UL) + ((uint64_t)op[5] * 17674925139268631843UL) + ((uint64_t)op[6] * 3155881917008880357UL) + ((((uint64_t)op[7] * 2502975817053249358UL) + ((uint64_t)op[8] * 2783709102958595942UL) + ((uint64_t)op[9] * 4303861093280866496UL) + ((uint64_t)op[10] * 12442483447733475858UL) + ((uint64_t)op[11] * 17574392635564835732UL)) * 18446744073709551609);
	tmp_q[7] = ((uint64_t)op[0] * 17574392635564835732UL) + ((uint64_t)op[1] * 13618301421045475485UL) + ((uint64_t)op[2] * 11562123458551789729UL) + ((uint64_t)op[3] * 16387303807734353009UL) + ((uint64_t)op[4] * 2443114270154127711UL) + ((uint64_t)op[5] * 18403820962515117585UL) + ((uint64_t)op[6] * 17674925139268631843UL) + ((uint64_t)op[7] * 3155881917008880357UL) + ((((uint64_t)op[8] * 2502975817053249358UL) + ((uint64_t)op[9] * 2783709102958595942UL) + ((uint64_t)op[10] * 4303861093280866496UL) + ((uint64_t)op[11] * 12442483447733475858UL)) * 18446744073709551609);
	tmp_q[8] = ((uint64_t)op[0] * 12442483447733475858UL) + ((uint64_t)op[1] * 17574392635564835732UL) + ((uint64_t)op[2] * 13618301421045475485UL) + ((uint64_t)op[3] * 11562123458551789729UL) + ((uint64_t)op[4] * 16387303807734353009UL) + ((uint64_t)op[5] * 2443114270154127711UL) + ((uint64_t)op[6] * 18403820962515117585UL) + ((uint64_t)op[7] * 17674925139268631843UL) + ((uint64_t)op[8] * 3155881917008880357UL) + ((((uint64_t)op[9] * 2502975817053249358UL) + ((uint64_t)op[10] * 2783709102958595942UL) + ((uint64_t)op[11] * 4303861093280866496UL)) * 18446744073709551609);
	tmp_q[9] = ((uint64_t)op[0] * 4303861093280866496UL) + ((uint64_t)op[1] * 12442483447733475858UL) + ((uint64_t)op[2] * 17574392635564835732UL) + ((uint64_t)op[3] * 13618301421045475485UL) + ((uint64_t)op[4] * 11562123458551789729UL) + ((uint64_t)op[5] * 16387303807734353009UL) + ((uint64_t)op[6] * 2443114270154127711UL) + ((uint64_t)op[7] * 18403820962515117585UL) + ((uint64_t)op[8] * 17674925139268631843UL) + ((uint64_t)op[9] * 3155881917008880357UL) + ((((uint64_t)op[10] * 2502975817053249358UL) + ((uint64_t)op[11] * 2783709102958595942UL)) * 18446744073709551609);
	tmp_q[10] = ((uint64_t)op[0] * 2783709102958595942UL) + ((uint64_t)op[1] * 4303861093280866496UL) + ((uint64_t)op[2] * 12442483447733475858UL) + ((uint64_t)op[3] * 17574392635564835732UL) + ((uint64_t)op[4] * 13618301421045475485UL) + ((uint64_t)op[5] * 11562123458551789729UL) + ((uint64_t)op[6] * 16387303807734353009UL) + ((uint64_t)op[7] * 2443114270154127711UL) + ((uint64_t)op[8] * 18403820962515117585UL) + ((uint64_t)op[9] * 17674925139268631843UL) + ((uint64_t)op[10] * 3155881917008880357UL) + ((uint64_t)op[11] * 925913354336806110UL);
	tmp_q[11] = ((uint64_t)op[0] * 2502975817053249358UL) + ((uint64_t)op[1] * 2783709102958595942UL) + ((uint64_t)op[2] * 4303861093280866496UL) + ((uint64_t)op[3] * 12442483447733475858UL) + ((uint64_t)op[4] * 17574392635564835732UL) + ((uint64_t)op[5] * 13618301421045475485UL) + ((uint64_t)op[6] * 11562123458551789729UL) + ((uint64_t)op[7] * 16387303807734353009UL) + ((uint64_t)op[8] * 2443114270154127711UL) + ((uint64_t)op[9] * 18403820962515117585UL) + ((uint64_t)op[10] * 17674925139268631843UL) + ((uint64_t)op[11] * 3155881917008880357UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1892237631903L) - ((-((int128)tmp_q[1] * 11772970253065L) + ((int128)tmp_q[2] * 6005532563278L) - ((int128)tmp_q[3] * 409069016295L) + ((int128)tmp_q[4] * 2244623580546L) - ((int128)tmp_q[5] * 3548493976523L) - ((int128)tmp_q[6] * 4569324298057L) - ((int128)tmp_q[7] * 2679799331518L) - ((int128)tmp_q[8] * 1018598267101L) - ((int128)tmp_q[9] * 4648091293512L) - ((int128)tmp_q[10] * 274693612973L) - ((int128)tmp_q[11] * 2718728308516L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 2718728308516L) + ((int128)tmp_q[1] * 1892237631903L) - ((-((int128)tmp_q[2] * 11772970253065L) + ((int128)tmp_q[3] * 6005532563278L) - ((int128)tmp_q[4] * 409069016295L) + ((int128)tmp_q[5] * 2244623580546L) - ((int128)tmp_q[6] * 3548493976523L) - ((int128)tmp_q[7] * 4569324298057L) - ((int128)tmp_q[8] * 2679799331518L) - ((int128)tmp_q[9] * 1018598267101L) - ((int128)tmp_q[10] * 4648091293512L) - ((int128)tmp_q[11] * 274693612973L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 274693612973L) - ((int128)tmp_q[1] * 2718728308516L) + ((int128)tmp_q[2] * 1892237631903L) - ((-((int128)tmp_q[3] * 11772970253065L) + ((int128)tmp_q[4] * 6005532563278L) - ((int128)tmp_q[5] * 409069016295L) + ((int128)tmp_q[6] * 2244623580546L) - ((int128)tmp_q[7] * 3548493976523L) - ((int128)tmp_q[8] * 4569324298057L) - ((int128)tmp_q[9] * 2679799331518L) - ((int128)tmp_q[10] * 1018598267101L) - ((int128)tmp_q[11] * 4648091293512L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 4648091293512L) - ((int128)tmp_q[1] * 274693612973L) - ((int128)tmp_q[2] * 2718728308516L) + ((int128)tmp_q[3] * 1892237631903L) - ((-((int128)tmp_q[4] * 11772970253065L) + ((int128)tmp_q[5] * 6005532563278L) - ((int128)tmp_q[6] * 409069016295L) + ((int128)tmp_q[7] * 2244623580546L) - ((int128)tmp_q[8] * 3548493976523L) - ((int128)tmp_q[9] * 4569324298057L) - ((int128)tmp_q[10] * 2679799331518L) - ((int128)tmp_q[11] * 1018598267101L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 1018598267101L) - ((int128)tmp_q[1] * 4648091293512L) - ((int128)tmp_q[2] * 274693612973L) - ((int128)tmp_q[3] * 2718728308516L) + ((int128)tmp_q[4] * 1892237631903L) - ((-((int128)tmp_q[5] * 11772970253065L) + ((int128)tmp_q[6] * 6005532563278L) - ((int128)tmp_q[7] * 409069016295L) + ((int128)tmp_q[8] * 2244623580546L) - ((int128)tmp_q[9] * 3548493976523L) - ((int128)tmp_q[10] * 4569324298057L) - ((int128)tmp_q[11] * 2679799331518L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 2679799331518L) - ((int128)tmp_q[1] * 1018598267101L) - ((int128)tmp_q[2] * 4648091293512L) - ((int128)tmp_q[3] * 274693612973L) - ((int128)tmp_q[4] * 2718728308516L) + ((int128)tmp_q[5] * 1892237631903L) - ((-((int128)tmp_q[6] * 11772970253065L) + ((int128)tmp_q[7] * 6005532563278L) - ((int128)tmp_q[8] * 409069016295L) + ((int128)tmp_q[9] * 2244623580546L) - ((int128)tmp_q[10] * 3548493976523L) - ((int128)tmp_q[11] * 4569324298057L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 4569324298057L) - ((int128)tmp_q[1] * 2679799331518L) - ((int128)tmp_q[2] * 1018598267101L) - ((int128)tmp_q[3] * 4648091293512L) - ((int128)tmp_q[4] * 274693612973L) - ((int128)tmp_q[5] * 2718728308516L) + ((int128)tmp_q[6] * 1892237631903L) - ((-((int128)tmp_q[7] * 11772970253065L) + ((int128)tmp_q[8] * 6005532563278L) - ((int128)tmp_q[9] * 409069016295L) + ((int128)tmp_q[10] * 2244623580546L) - ((int128)tmp_q[11] * 3548493976523L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 3548493976523L) - ((int128)tmp_q[1] * 4569324298057L) - ((int128)tmp_q[2] * 2679799331518L) - ((int128)tmp_q[3] * 1018598267101L) - ((int128)tmp_q[4] * 4648091293512L) - ((int128)tmp_q[5] * 274693612973L) - ((int128)tmp_q[6] * 2718728308516L) + ((int128)tmp_q[7] * 1892237631903L) - ((-((int128)tmp_q[8] * 11772970253065L) + ((int128)tmp_q[9] * 6005532563278L) - ((int128)tmp_q[10] * 409069016295L) + ((int128)tmp_q[11] * 2244623580546L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 2244623580546L) - ((int128)tmp_q[1] * 3548493976523L) - ((int128)tmp_q[2] * 4569324298057L) - ((int128)tmp_q[3] * 2679799331518L) - ((int128)tmp_q[4] * 1018598267101L) - ((int128)tmp_q[5] * 4648091293512L) - ((int128)tmp_q[6] * 274693612973L) - ((int128)tmp_q[7] * 2718728308516L) + ((int128)tmp_q[8] * 1892237631903L) - ((-((int128)tmp_q[9] * 11772970253065L) + ((int128)tmp_q[10] * 6005532563278L) - ((int128)tmp_q[11] * 409069016295L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 409069016295L) + ((int128)tmp_q[1] * 2244623580546L) - ((int128)tmp_q[2] * 3548493976523L) - ((int128)tmp_q[3] * 4569324298057L) - ((int128)tmp_q[4] * 2679799331518L) - ((int128)tmp_q[5] * 1018598267101L) - ((int128)tmp_q[6] * 4648091293512L) - ((int128)tmp_q[7] * 274693612973L) - ((int128)tmp_q[8] * 2718728308516L) + ((int128)tmp_q[9] * 1892237631903L) - ((-((int128)tmp_q[10] * 11772970253065L) + ((int128)tmp_q[11] * 6005532563278L)) * 7);
	tmp_zero[10] = ((int128)tmp_q[0] * 6005532563278L) - ((int128)tmp_q[1] * 409069016295L) + ((int128)tmp_q[2] * 2244623580546L) - ((int128)tmp_q[3] * 3548493976523L) - ((int128)tmp_q[4] * 4569324298057L) - ((int128)tmp_q[5] * 2679799331518L) - ((int128)tmp_q[6] * 1018598267101L) - ((int128)tmp_q[7] * 4648091293512L) - ((int128)tmp_q[8] * 274693612973L) - ((int128)tmp_q[9] * 2718728308516L) + ((int128)tmp_q[10] * 1892237631903L) + ((int128)tmp_q[11] * 82410791771455L);
	tmp_zero[11] = -((int128)tmp_q[0] * 11772970253065L) + ((int128)tmp_q[1] * 6005532563278L) - ((int128)tmp_q[2] * 409069016295L) + ((int128)tmp_q[3] * 2244623580546L) - ((int128)tmp_q[4] * 3548493976523L) - ((int128)tmp_q[5] * 4569324298057L) - ((int128)tmp_q[6] * 2679799331518L) - ((int128)tmp_q[7] * 1018598267101L) - ((int128)tmp_q[8] * 4648091293512L) - ((int128)tmp_q[9] * 274693612973L) - ((int128)tmp_q[10] * 2718728308516L) + ((int128)tmp_q[11] * 1892237631903L);

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

