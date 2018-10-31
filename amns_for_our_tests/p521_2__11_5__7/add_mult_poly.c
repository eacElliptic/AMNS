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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4643019197462740508UL) + ((((uint64_t)op[1] * 3045064833138833494UL) + ((uint64_t)op[2] * 14015256880692556711UL) + ((uint64_t)op[3] * 1509847002612132585UL) + ((uint64_t)op[4] * 4990943503534666822UL) + ((uint64_t)op[5] * 8279085017251641886UL) + ((uint64_t)op[6] * 1479629907338073194UL) + ((uint64_t)op[7] * 9615137692728010420UL) + ((uint64_t)op[8] * 15439977747366564160UL) + ((uint64_t)op[9] * 2866864371958093822UL) + ((uint64_t)op[10] * 17978139908939332339UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 17978139908939332339UL) + ((uint64_t)op[1] * 4643019197462740508UL) + ((((uint64_t)op[2] * 3045064833138833494UL) + ((uint64_t)op[3] * 14015256880692556711UL) + ((uint64_t)op[4] * 1509847002612132585UL) + ((uint64_t)op[5] * 4990943503534666822UL) + ((uint64_t)op[6] * 8279085017251641886UL) + ((uint64_t)op[7] * 1479629907338073194UL) + ((uint64_t)op[8] * 9615137692728010420UL) + ((uint64_t)op[9] * 15439977747366564160UL) + ((uint64_t)op[10] * 2866864371958093822UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 2866864371958093822UL) + ((uint64_t)op[1] * 17978139908939332339UL) + ((uint64_t)op[2] * 4643019197462740508UL) + ((((uint64_t)op[3] * 3045064833138833494UL) + ((uint64_t)op[4] * 14015256880692556711UL) + ((uint64_t)op[5] * 1509847002612132585UL) + ((uint64_t)op[6] * 4990943503534666822UL) + ((uint64_t)op[7] * 8279085017251641886UL) + ((uint64_t)op[8] * 1479629907338073194UL) + ((uint64_t)op[9] * 9615137692728010420UL) + ((uint64_t)op[10] * 15439977747366564160UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 15439977747366564160UL) + ((uint64_t)op[1] * 2866864371958093822UL) + ((uint64_t)op[2] * 17978139908939332339UL) + ((uint64_t)op[3] * 4643019197462740508UL) + ((((uint64_t)op[4] * 3045064833138833494UL) + ((uint64_t)op[5] * 14015256880692556711UL) + ((uint64_t)op[6] * 1509847002612132585UL) + ((uint64_t)op[7] * 4990943503534666822UL) + ((uint64_t)op[8] * 8279085017251641886UL) + ((uint64_t)op[9] * 1479629907338073194UL) + ((uint64_t)op[10] * 9615137692728010420UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 9615137692728010420UL) + ((uint64_t)op[1] * 15439977747366564160UL) + ((uint64_t)op[2] * 2866864371958093822UL) + ((uint64_t)op[3] * 17978139908939332339UL) + ((uint64_t)op[4] * 4643019197462740508UL) + ((((uint64_t)op[5] * 3045064833138833494UL) + ((uint64_t)op[6] * 14015256880692556711UL) + ((uint64_t)op[7] * 1509847002612132585UL) + ((uint64_t)op[8] * 4990943503534666822UL) + ((uint64_t)op[9] * 8279085017251641886UL) + ((uint64_t)op[10] * 1479629907338073194UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 1479629907338073194UL) + ((uint64_t)op[1] * 9615137692728010420UL) + ((uint64_t)op[2] * 15439977747366564160UL) + ((uint64_t)op[3] * 2866864371958093822UL) + ((uint64_t)op[4] * 17978139908939332339UL) + ((uint64_t)op[5] * 4643019197462740508UL) + ((((uint64_t)op[6] * 3045064833138833494UL) + ((uint64_t)op[7] * 14015256880692556711UL) + ((uint64_t)op[8] * 1509847002612132585UL) + ((uint64_t)op[9] * 4990943503534666822UL) + ((uint64_t)op[10] * 8279085017251641886UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 8279085017251641886UL) + ((uint64_t)op[1] * 1479629907338073194UL) + ((uint64_t)op[2] * 9615137692728010420UL) + ((uint64_t)op[3] * 15439977747366564160UL) + ((uint64_t)op[4] * 2866864371958093822UL) + ((uint64_t)op[5] * 17978139908939332339UL) + ((uint64_t)op[6] * 4643019197462740508UL) + ((((uint64_t)op[7] * 3045064833138833494UL) + ((uint64_t)op[8] * 14015256880692556711UL) + ((uint64_t)op[9] * 1509847002612132585UL) + ((uint64_t)op[10] * 4990943503534666822UL)) * 5);
	tmp_q[7] = ((uint64_t)op[0] * 4990943503534666822UL) + ((uint64_t)op[1] * 8279085017251641886UL) + ((uint64_t)op[2] * 1479629907338073194UL) + ((uint64_t)op[3] * 9615137692728010420UL) + ((uint64_t)op[4] * 15439977747366564160UL) + ((uint64_t)op[5] * 2866864371958093822UL) + ((uint64_t)op[6] * 17978139908939332339UL) + ((uint64_t)op[7] * 4643019197462740508UL) + ((((uint64_t)op[8] * 3045064833138833494UL) + ((uint64_t)op[9] * 14015256880692556711UL) + ((uint64_t)op[10] * 1509847002612132585UL)) * 5);
	tmp_q[8] = ((uint64_t)op[0] * 1509847002612132585UL) + ((uint64_t)op[1] * 4990943503534666822UL) + ((uint64_t)op[2] * 8279085017251641886UL) + ((uint64_t)op[3] * 1479629907338073194UL) + ((uint64_t)op[4] * 9615137692728010420UL) + ((uint64_t)op[5] * 15439977747366564160UL) + ((uint64_t)op[6] * 2866864371958093822UL) + ((uint64_t)op[7] * 17978139908939332339UL) + ((uint64_t)op[8] * 4643019197462740508UL) + ((((uint64_t)op[9] * 3045064833138833494UL) + ((uint64_t)op[10] * 14015256880692556711UL)) * 5);
	tmp_q[9] = ((uint64_t)op[0] * 14015256880692556711UL) + ((uint64_t)op[1] * 1509847002612132585UL) + ((uint64_t)op[2] * 4990943503534666822UL) + ((uint64_t)op[3] * 8279085017251641886UL) + ((uint64_t)op[4] * 1479629907338073194UL) + ((uint64_t)op[5] * 9615137692728010420UL) + ((uint64_t)op[6] * 15439977747366564160UL) + ((uint64_t)op[7] * 2866864371958093822UL) + ((uint64_t)op[8] * 17978139908939332339UL) + ((uint64_t)op[9] * 4643019197462740508UL) + ((uint64_t)op[10] * 15225324165694167470UL);
	tmp_q[10] = ((uint64_t)op[0] * 3045064833138833494UL) + ((uint64_t)op[1] * 14015256880692556711UL) + ((uint64_t)op[2] * 1509847002612132585UL) + ((uint64_t)op[3] * 4990943503534666822UL) + ((uint64_t)op[4] * 8279085017251641886UL) + ((uint64_t)op[5] * 1479629907338073194UL) + ((uint64_t)op[6] * 9615137692728010420UL) + ((uint64_t)op[7] * 15439977747366564160UL) + ((uint64_t)op[8] * 2866864371958093822UL) + ((uint64_t)op[9] * 17978139908939332339UL) + ((uint64_t)op[10] * 4643019197462740508UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 67758780336604L) + ((-((int128)tmp_q[1] * 66115779007995L) + ((int128)tmp_q[2] * 88964214607219L) + ((int128)tmp_q[3] * 4215989045762L) + ((int128)tmp_q[4] * 26766804843989L) - ((int128)tmp_q[5] * 10601243006566L) + ((int128)tmp_q[6] * 5531213049491L) + ((int128)tmp_q[7] * 24388825508183L) - ((int128)tmp_q[8] * 16043182575793L) + ((int128)tmp_q[9] * 64964093907225L) + ((int128)tmp_q[10] * 64797792294508L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 64797792294508L) + ((int128)tmp_q[1] * 67758780336604L) + ((-((int128)tmp_q[2] * 66115779007995L) + ((int128)tmp_q[3] * 88964214607219L) + ((int128)tmp_q[4] * 4215989045762L) + ((int128)tmp_q[5] * 26766804843989L) - ((int128)tmp_q[6] * 10601243006566L) + ((int128)tmp_q[7] * 5531213049491L) + ((int128)tmp_q[8] * 24388825508183L) - ((int128)tmp_q[9] * 16043182575793L) + ((int128)tmp_q[10] * 64964093907225L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 64964093907225L) + ((int128)tmp_q[1] * 64797792294508L) + ((int128)tmp_q[2] * 67758780336604L) + ((-((int128)tmp_q[3] * 66115779007995L) + ((int128)tmp_q[4] * 88964214607219L) + ((int128)tmp_q[5] * 4215989045762L) + ((int128)tmp_q[6] * 26766804843989L) - ((int128)tmp_q[7] * 10601243006566L) + ((int128)tmp_q[8] * 5531213049491L) + ((int128)tmp_q[9] * 24388825508183L) - ((int128)tmp_q[10] * 16043182575793L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 16043182575793L) + ((int128)tmp_q[1] * 64964093907225L) + ((int128)tmp_q[2] * 64797792294508L) + ((int128)tmp_q[3] * 67758780336604L) + ((-((int128)tmp_q[4] * 66115779007995L) + ((int128)tmp_q[5] * 88964214607219L) + ((int128)tmp_q[6] * 4215989045762L) + ((int128)tmp_q[7] * 26766804843989L) - ((int128)tmp_q[8] * 10601243006566L) + ((int128)tmp_q[9] * 5531213049491L) + ((int128)tmp_q[10] * 24388825508183L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 24388825508183L) - ((int128)tmp_q[1] * 16043182575793L) + ((int128)tmp_q[2] * 64964093907225L) + ((int128)tmp_q[3] * 64797792294508L) + ((int128)tmp_q[4] * 67758780336604L) + ((-((int128)tmp_q[5] * 66115779007995L) + ((int128)tmp_q[6] * 88964214607219L) + ((int128)tmp_q[7] * 4215989045762L) + ((int128)tmp_q[8] * 26766804843989L) - ((int128)tmp_q[9] * 10601243006566L) + ((int128)tmp_q[10] * 5531213049491L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 5531213049491L) + ((int128)tmp_q[1] * 24388825508183L) - ((int128)tmp_q[2] * 16043182575793L) + ((int128)tmp_q[3] * 64964093907225L) + ((int128)tmp_q[4] * 64797792294508L) + ((int128)tmp_q[5] * 67758780336604L) + ((-((int128)tmp_q[6] * 66115779007995L) + ((int128)tmp_q[7] * 88964214607219L) + ((int128)tmp_q[8] * 4215989045762L) + ((int128)tmp_q[9] * 26766804843989L) - ((int128)tmp_q[10] * 10601243006566L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 10601243006566L) + ((int128)tmp_q[1] * 5531213049491L) + ((int128)tmp_q[2] * 24388825508183L) - ((int128)tmp_q[3] * 16043182575793L) + ((int128)tmp_q[4] * 64964093907225L) + ((int128)tmp_q[5] * 64797792294508L) + ((int128)tmp_q[6] * 67758780336604L) + ((-((int128)tmp_q[7] * 66115779007995L) + ((int128)tmp_q[8] * 88964214607219L) + ((int128)tmp_q[9] * 4215989045762L) + ((int128)tmp_q[10] * 26766804843989L)) * 5);
	tmp_zero[7] = ((int128)tmp_q[0] * 26766804843989L) - ((int128)tmp_q[1] * 10601243006566L) + ((int128)tmp_q[2] * 5531213049491L) + ((int128)tmp_q[3] * 24388825508183L) - ((int128)tmp_q[4] * 16043182575793L) + ((int128)tmp_q[5] * 64964093907225L) + ((int128)tmp_q[6] * 64797792294508L) + ((int128)tmp_q[7] * 67758780336604L) + ((-((int128)tmp_q[8] * 66115779007995L) + ((int128)tmp_q[9] * 88964214607219L) + ((int128)tmp_q[10] * 4215989045762L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 4215989045762L) + ((int128)tmp_q[1] * 26766804843989L) - ((int128)tmp_q[2] * 10601243006566L) + ((int128)tmp_q[3] * 5531213049491L) + ((int128)tmp_q[4] * 24388825508183L) - ((int128)tmp_q[5] * 16043182575793L) + ((int128)tmp_q[6] * 64964093907225L) + ((int128)tmp_q[7] * 64797792294508L) + ((int128)tmp_q[8] * 67758780336604L) + ((-((int128)tmp_q[9] * 66115779007995L) + ((int128)tmp_q[10] * 88964214607219L)) * 5);
	tmp_zero[9] = ((int128)tmp_q[0] * 88964214607219L) + ((int128)tmp_q[1] * 4215989045762L) + ((int128)tmp_q[2] * 26766804843989L) - ((int128)tmp_q[3] * 10601243006566L) + ((int128)tmp_q[4] * 5531213049491L) + ((int128)tmp_q[5] * 24388825508183L) - ((int128)tmp_q[6] * 16043182575793L) + ((int128)tmp_q[7] * 64964093907225L) + ((int128)tmp_q[8] * 64797792294508L) + ((int128)tmp_q[9] * 67758780336604L) - ((int128)tmp_q[10] * 330578895039975L);
	tmp_zero[10] = -((int128)tmp_q[0] * 66115779007995L) + ((int128)tmp_q[1] * 88964214607219L) + ((int128)tmp_q[2] * 4215989045762L) + ((int128)tmp_q[3] * 26766804843989L) - ((int128)tmp_q[4] * 10601243006566L) + ((int128)tmp_q[5] * 5531213049491L) + ((int128)tmp_q[6] * 24388825508183L) - ((int128)tmp_q[7] * 16043182575793L) + ((int128)tmp_q[8] * 64964093907225L) + ((int128)tmp_q[9] * 64797792294508L) + ((int128)tmp_q[10] * 67758780336604L);

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

