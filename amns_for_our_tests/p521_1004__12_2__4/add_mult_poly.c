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
	tmp_q[0] = ((uint64_t)op[0] * 2956371582337097261UL) + ((((uint64_t)op[1] * 14609026992064133095UL) + ((uint64_t)op[2] * 12635142134531874667UL) + ((uint64_t)op[3] * 6625369405568193854UL) + ((uint64_t)op[4] * 17021983163601314319UL) + ((uint64_t)op[5] * 1429927824372220368UL) + ((uint64_t)op[6] * 6960363359420198973UL) + ((uint64_t)op[7] * 1665394132745982568UL) + ((uint64_t)op[8] * 13158778625855637090UL) + ((uint64_t)op[9] * 6539210436102533557UL) + ((uint64_t)op[10] * 7619642665193010426UL) + ((uint64_t)op[11] * 16187228200983445539UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 16187228200983445539UL) + ((uint64_t)op[1] * 2956371582337097261UL) + ((((uint64_t)op[2] * 14609026992064133095UL) + ((uint64_t)op[3] * 12635142134531874667UL) + ((uint64_t)op[4] * 6625369405568193854UL) + ((uint64_t)op[5] * 17021983163601314319UL) + ((uint64_t)op[6] * 1429927824372220368UL) + ((uint64_t)op[7] * 6960363359420198973UL) + ((uint64_t)op[8] * 1665394132745982568UL) + ((uint64_t)op[9] * 13158778625855637090UL) + ((uint64_t)op[10] * 6539210436102533557UL) + ((uint64_t)op[11] * 7619642665193010426UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 7619642665193010426UL) + ((uint64_t)op[1] * 16187228200983445539UL) + ((uint64_t)op[2] * 2956371582337097261UL) + ((((uint64_t)op[3] * 14609026992064133095UL) + ((uint64_t)op[4] * 12635142134531874667UL) + ((uint64_t)op[5] * 6625369405568193854UL) + ((uint64_t)op[6] * 17021983163601314319UL) + ((uint64_t)op[7] * 1429927824372220368UL) + ((uint64_t)op[8] * 6960363359420198973UL) + ((uint64_t)op[9] * 1665394132745982568UL) + ((uint64_t)op[10] * 13158778625855637090UL) + ((uint64_t)op[11] * 6539210436102533557UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 6539210436102533557UL) + ((uint64_t)op[1] * 7619642665193010426UL) + ((uint64_t)op[2] * 16187228200983445539UL) + ((uint64_t)op[3] * 2956371582337097261UL) + ((((uint64_t)op[4] * 14609026992064133095UL) + ((uint64_t)op[5] * 12635142134531874667UL) + ((uint64_t)op[6] * 6625369405568193854UL) + ((uint64_t)op[7] * 17021983163601314319UL) + ((uint64_t)op[8] * 1429927824372220368UL) + ((uint64_t)op[9] * 6960363359420198973UL) + ((uint64_t)op[10] * 1665394132745982568UL) + ((uint64_t)op[11] * 13158778625855637090UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 13158778625855637090UL) + ((uint64_t)op[1] * 6539210436102533557UL) + ((uint64_t)op[2] * 7619642665193010426UL) + ((uint64_t)op[3] * 16187228200983445539UL) + ((uint64_t)op[4] * 2956371582337097261UL) + ((((uint64_t)op[5] * 14609026992064133095UL) + ((uint64_t)op[6] * 12635142134531874667UL) + ((uint64_t)op[7] * 6625369405568193854UL) + ((uint64_t)op[8] * 17021983163601314319UL) + ((uint64_t)op[9] * 1429927824372220368UL) + ((uint64_t)op[10] * 6960363359420198973UL) + ((uint64_t)op[11] * 1665394132745982568UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 1665394132745982568UL) + ((uint64_t)op[1] * 13158778625855637090UL) + ((uint64_t)op[2] * 6539210436102533557UL) + ((uint64_t)op[3] * 7619642665193010426UL) + ((uint64_t)op[4] * 16187228200983445539UL) + ((uint64_t)op[5] * 2956371582337097261UL) + ((((uint64_t)op[6] * 14609026992064133095UL) + ((uint64_t)op[7] * 12635142134531874667UL) + ((uint64_t)op[8] * 6625369405568193854UL) + ((uint64_t)op[9] * 17021983163601314319UL) + ((uint64_t)op[10] * 1429927824372220368UL) + ((uint64_t)op[11] * 6960363359420198973UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 6960363359420198973UL) + ((uint64_t)op[1] * 1665394132745982568UL) + ((uint64_t)op[2] * 13158778625855637090UL) + ((uint64_t)op[3] * 6539210436102533557UL) + ((uint64_t)op[4] * 7619642665193010426UL) + ((uint64_t)op[5] * 16187228200983445539UL) + ((uint64_t)op[6] * 2956371582337097261UL) + ((((uint64_t)op[7] * 14609026992064133095UL) + ((uint64_t)op[8] * 12635142134531874667UL) + ((uint64_t)op[9] * 6625369405568193854UL) + ((uint64_t)op[10] * 17021983163601314319UL) + ((uint64_t)op[11] * 1429927824372220368UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 1429927824372220368UL) + ((uint64_t)op[1] * 6960363359420198973UL) + ((uint64_t)op[2] * 1665394132745982568UL) + ((uint64_t)op[3] * 13158778625855637090UL) + ((uint64_t)op[4] * 6539210436102533557UL) + ((uint64_t)op[5] * 7619642665193010426UL) + ((uint64_t)op[6] * 16187228200983445539UL) + ((uint64_t)op[7] * 2956371582337097261UL) + ((((uint64_t)op[8] * 14609026992064133095UL) + ((uint64_t)op[9] * 12635142134531874667UL) + ((uint64_t)op[10] * 6625369405568193854UL) + ((uint64_t)op[11] * 17021983163601314319UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 17021983163601314319UL) + ((uint64_t)op[1] * 1429927824372220368UL) + ((uint64_t)op[2] * 6960363359420198973UL) + ((uint64_t)op[3] * 1665394132745982568UL) + ((uint64_t)op[4] * 13158778625855637090UL) + ((uint64_t)op[5] * 6539210436102533557UL) + ((uint64_t)op[6] * 7619642665193010426UL) + ((uint64_t)op[7] * 16187228200983445539UL) + ((uint64_t)op[8] * 2956371582337097261UL) + ((((uint64_t)op[9] * 14609026992064133095UL) + ((uint64_t)op[10] * 12635142134531874667UL) + ((uint64_t)op[11] * 6625369405568193854UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 6625369405568193854UL) + ((uint64_t)op[1] * 17021983163601314319UL) + ((uint64_t)op[2] * 1429927824372220368UL) + ((uint64_t)op[3] * 6960363359420198973UL) + ((uint64_t)op[4] * 1665394132745982568UL) + ((uint64_t)op[5] * 13158778625855637090UL) + ((uint64_t)op[6] * 6539210436102533557UL) + ((uint64_t)op[7] * 7619642665193010426UL) + ((uint64_t)op[8] * 16187228200983445539UL) + ((uint64_t)op[9] * 2956371582337097261UL) + ((((uint64_t)op[10] * 14609026992064133095UL) + ((uint64_t)op[11] * 12635142134531874667UL)) * 2);
	tmp_q[10] = ((uint64_t)op[0] * 12635142134531874667UL) + ((uint64_t)op[1] * 6625369405568193854UL) + ((uint64_t)op[2] * 17021983163601314319UL) + ((uint64_t)op[3] * 1429927824372220368UL) + ((uint64_t)op[4] * 6960363359420198973UL) + ((uint64_t)op[5] * 1665394132745982568UL) + ((uint64_t)op[6] * 13158778625855637090UL) + ((uint64_t)op[7] * 6539210436102533557UL) + ((uint64_t)op[8] * 7619642665193010426UL) + ((uint64_t)op[9] * 16187228200983445539UL) + ((uint64_t)op[10] * 2956371582337097261UL) + ((uint64_t)op[11] * 10771309910418714574UL);
	tmp_q[11] = ((uint64_t)op[0] * 14609026992064133095UL) + ((uint64_t)op[1] * 12635142134531874667UL) + ((uint64_t)op[2] * 6625369405568193854UL) + ((uint64_t)op[3] * 17021983163601314319UL) + ((uint64_t)op[4] * 1429927824372220368UL) + ((uint64_t)op[5] * 6960363359420198973UL) + ((uint64_t)op[6] * 1665394132745982568UL) + ((uint64_t)op[7] * 13158778625855637090UL) + ((uint64_t)op[8] * 6539210436102533557UL) + ((uint64_t)op[9] * 7619642665193010426UL) + ((uint64_t)op[10] * 16187228200983445539UL) + ((uint64_t)op[11] * 2956371582337097261UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 996356631689L) + ((-((int128)tmp_q[1] * 986147862756L) + ((int128)tmp_q[2] * 1463342671939L) - ((int128)tmp_q[3] * 335588447957L) - ((int128)tmp_q[4] * 18203900003L) - ((int128)tmp_q[5] * 4495508596719L) - ((int128)tmp_q[6] * 6644650019967L) + ((int128)tmp_q[7] * 2749281266010L) - ((int128)tmp_q[8] * 1834900190451L) + ((int128)tmp_q[9] * 3344209626198L) + ((int128)tmp_q[10] * 739522334343L) - ((int128)tmp_q[11] * 589257295569L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 589257295569L) + ((int128)tmp_q[1] * 996356631689L) + ((-((int128)tmp_q[2] * 986147862756L) + ((int128)tmp_q[3] * 1463342671939L) - ((int128)tmp_q[4] * 335588447957L) - ((int128)tmp_q[5] * 18203900003L) - ((int128)tmp_q[6] * 4495508596719L) - ((int128)tmp_q[7] * 6644650019967L) + ((int128)tmp_q[8] * 2749281266010L) - ((int128)tmp_q[9] * 1834900190451L) + ((int128)tmp_q[10] * 3344209626198L) + ((int128)tmp_q[11] * 739522334343L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 739522334343L) - ((int128)tmp_q[1] * 589257295569L) + ((int128)tmp_q[2] * 996356631689L) + ((-((int128)tmp_q[3] * 986147862756L) + ((int128)tmp_q[4] * 1463342671939L) - ((int128)tmp_q[5] * 335588447957L) - ((int128)tmp_q[6] * 18203900003L) - ((int128)tmp_q[7] * 4495508596719L) - ((int128)tmp_q[8] * 6644650019967L) + ((int128)tmp_q[9] * 2749281266010L) - ((int128)tmp_q[10] * 1834900190451L) + ((int128)tmp_q[11] * 3344209626198L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 3344209626198L) + ((int128)tmp_q[1] * 739522334343L) - ((int128)tmp_q[2] * 589257295569L) + ((int128)tmp_q[3] * 996356631689L) + ((-((int128)tmp_q[4] * 986147862756L) + ((int128)tmp_q[5] * 1463342671939L) - ((int128)tmp_q[6] * 335588447957L) - ((int128)tmp_q[7] * 18203900003L) - ((int128)tmp_q[8] * 4495508596719L) - ((int128)tmp_q[9] * 6644650019967L) + ((int128)tmp_q[10] * 2749281266010L) - ((int128)tmp_q[11] * 1834900190451L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1834900190451L) + ((int128)tmp_q[1] * 3344209626198L) + ((int128)tmp_q[2] * 739522334343L) - ((int128)tmp_q[3] * 589257295569L) + ((int128)tmp_q[4] * 996356631689L) + ((-((int128)tmp_q[5] * 986147862756L) + ((int128)tmp_q[6] * 1463342671939L) - ((int128)tmp_q[7] * 335588447957L) - ((int128)tmp_q[8] * 18203900003L) - ((int128)tmp_q[9] * 4495508596719L) - ((int128)tmp_q[10] * 6644650019967L) + ((int128)tmp_q[11] * 2749281266010L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 2749281266010L) - ((int128)tmp_q[1] * 1834900190451L) + ((int128)tmp_q[2] * 3344209626198L) + ((int128)tmp_q[3] * 739522334343L) - ((int128)tmp_q[4] * 589257295569L) + ((int128)tmp_q[5] * 996356631689L) + ((-((int128)tmp_q[6] * 986147862756L) + ((int128)tmp_q[7] * 1463342671939L) - ((int128)tmp_q[8] * 335588447957L) - ((int128)tmp_q[9] * 18203900003L) - ((int128)tmp_q[10] * 4495508596719L) - ((int128)tmp_q[11] * 6644650019967L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 6644650019967L) + ((int128)tmp_q[1] * 2749281266010L) - ((int128)tmp_q[2] * 1834900190451L) + ((int128)tmp_q[3] * 3344209626198L) + ((int128)tmp_q[4] * 739522334343L) - ((int128)tmp_q[5] * 589257295569L) + ((int128)tmp_q[6] * 996356631689L) + ((-((int128)tmp_q[7] * 986147862756L) + ((int128)tmp_q[8] * 1463342671939L) - ((int128)tmp_q[9] * 335588447957L) - ((int128)tmp_q[10] * 18203900003L) - ((int128)tmp_q[11] * 4495508596719L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 4495508596719L) - ((int128)tmp_q[1] * 6644650019967L) + ((int128)tmp_q[2] * 2749281266010L) - ((int128)tmp_q[3] * 1834900190451L) + ((int128)tmp_q[4] * 3344209626198L) + ((int128)tmp_q[5] * 739522334343L) - ((int128)tmp_q[6] * 589257295569L) + ((int128)tmp_q[7] * 996356631689L) + ((-((int128)tmp_q[8] * 986147862756L) + ((int128)tmp_q[9] * 1463342671939L) - ((int128)tmp_q[10] * 335588447957L) - ((int128)tmp_q[11] * 18203900003L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 18203900003L) - ((int128)tmp_q[1] * 4495508596719L) - ((int128)tmp_q[2] * 6644650019967L) + ((int128)tmp_q[3] * 2749281266010L) - ((int128)tmp_q[4] * 1834900190451L) + ((int128)tmp_q[5] * 3344209626198L) + ((int128)tmp_q[6] * 739522334343L) - ((int128)tmp_q[7] * 589257295569L) + ((int128)tmp_q[8] * 996356631689L) + ((-((int128)tmp_q[9] * 986147862756L) + ((int128)tmp_q[10] * 1463342671939L) - ((int128)tmp_q[11] * 335588447957L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 335588447957L) - ((int128)tmp_q[1] * 18203900003L) - ((int128)tmp_q[2] * 4495508596719L) - ((int128)tmp_q[3] * 6644650019967L) + ((int128)tmp_q[4] * 2749281266010L) - ((int128)tmp_q[5] * 1834900190451L) + ((int128)tmp_q[6] * 3344209626198L) + ((int128)tmp_q[7] * 739522334343L) - ((int128)tmp_q[8] * 589257295569L) + ((int128)tmp_q[9] * 996356631689L) + ((-((int128)tmp_q[10] * 986147862756L) + ((int128)tmp_q[11] * 1463342671939L)) * 2);
	tmp_zero[10] = ((int128)tmp_q[0] * 1463342671939L) - ((int128)tmp_q[1] * 335588447957L) - ((int128)tmp_q[2] * 18203900003L) - ((int128)tmp_q[3] * 4495508596719L) - ((int128)tmp_q[4] * 6644650019967L) + ((int128)tmp_q[5] * 2749281266010L) - ((int128)tmp_q[6] * 1834900190451L) + ((int128)tmp_q[7] * 3344209626198L) + ((int128)tmp_q[8] * 739522334343L) - ((int128)tmp_q[9] * 589257295569L) + ((int128)tmp_q[10] * 996356631689L) - ((int128)tmp_q[11] * 1972295725512L);
	tmp_zero[11] = -((int128)tmp_q[0] * 986147862756L) + ((int128)tmp_q[1] * 1463342671939L) - ((int128)tmp_q[2] * 335588447957L) - ((int128)tmp_q[3] * 18203900003L) - ((int128)tmp_q[4] * 4495508596719L) - ((int128)tmp_q[5] * 6644650019967L) + ((int128)tmp_q[6] * 2749281266010L) - ((int128)tmp_q[7] * 1834900190451L) + ((int128)tmp_q[8] * 3344209626198L) + ((int128)tmp_q[9] * 739522334343L) - ((int128)tmp_q[10] * 589257295569L) + ((int128)tmp_q[11] * 996356631689L);

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

