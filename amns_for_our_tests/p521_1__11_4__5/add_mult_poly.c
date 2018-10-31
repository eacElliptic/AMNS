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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5851854850093054879UL) + ((((uint64_t)op[1] * 1986626803087918625UL) + ((uint64_t)op[2] * 1170269630229199703UL) + ((uint64_t)op[3] * 15062026177312334876UL) + ((uint64_t)op[4] * 12555398568625425479UL) + ((uint64_t)op[5] * 9226291044724333206UL) + ((uint64_t)op[6] * 12070639868163715111UL) + ((uint64_t)op[7] * 1093162736320069068UL) + ((uint64_t)op[8] * 8887331613627784179UL) + ((uint64_t)op[9] * 8402347405076345566UL) + ((uint64_t)op[10] * 9664495871165975587UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 9664495871165975587UL) + ((uint64_t)op[1] * 5851854850093054879UL) + ((((uint64_t)op[2] * 1986626803087918625UL) + ((uint64_t)op[3] * 1170269630229199703UL) + ((uint64_t)op[4] * 15062026177312334876UL) + ((uint64_t)op[5] * 12555398568625425479UL) + ((uint64_t)op[6] * 9226291044724333206UL) + ((uint64_t)op[7] * 12070639868163715111UL) + ((uint64_t)op[8] * 1093162736320069068UL) + ((uint64_t)op[9] * 8887331613627784179UL) + ((uint64_t)op[10] * 8402347405076345566UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 8402347405076345566UL) + ((uint64_t)op[1] * 9664495871165975587UL) + ((uint64_t)op[2] * 5851854850093054879UL) + ((((uint64_t)op[3] * 1986626803087918625UL) + ((uint64_t)op[4] * 1170269630229199703UL) + ((uint64_t)op[5] * 15062026177312334876UL) + ((uint64_t)op[6] * 12555398568625425479UL) + ((uint64_t)op[7] * 9226291044724333206UL) + ((uint64_t)op[8] * 12070639868163715111UL) + ((uint64_t)op[9] * 1093162736320069068UL) + ((uint64_t)op[10] * 8887331613627784179UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 8887331613627784179UL) + ((uint64_t)op[1] * 8402347405076345566UL) + ((uint64_t)op[2] * 9664495871165975587UL) + ((uint64_t)op[3] * 5851854850093054879UL) + ((((uint64_t)op[4] * 1986626803087918625UL) + ((uint64_t)op[5] * 1170269630229199703UL) + ((uint64_t)op[6] * 15062026177312334876UL) + ((uint64_t)op[7] * 12555398568625425479UL) + ((uint64_t)op[8] * 9226291044724333206UL) + ((uint64_t)op[9] * 12070639868163715111UL) + ((uint64_t)op[10] * 1093162736320069068UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 1093162736320069068UL) + ((uint64_t)op[1] * 8887331613627784179UL) + ((uint64_t)op[2] * 8402347405076345566UL) + ((uint64_t)op[3] * 9664495871165975587UL) + ((uint64_t)op[4] * 5851854850093054879UL) + ((((uint64_t)op[5] * 1986626803087918625UL) + ((uint64_t)op[6] * 1170269630229199703UL) + ((uint64_t)op[7] * 15062026177312334876UL) + ((uint64_t)op[8] * 12555398568625425479UL) + ((uint64_t)op[9] * 9226291044724333206UL) + ((uint64_t)op[10] * 12070639868163715111UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 12070639868163715111UL) + ((uint64_t)op[1] * 1093162736320069068UL) + ((uint64_t)op[2] * 8887331613627784179UL) + ((uint64_t)op[3] * 8402347405076345566UL) + ((uint64_t)op[4] * 9664495871165975587UL) + ((uint64_t)op[5] * 5851854850093054879UL) + ((((uint64_t)op[6] * 1986626803087918625UL) + ((uint64_t)op[7] * 1170269630229199703UL) + ((uint64_t)op[8] * 15062026177312334876UL) + ((uint64_t)op[9] * 12555398568625425479UL) + ((uint64_t)op[10] * 9226291044724333206UL)) * 4);
	tmp_q[6] = ((uint64_t)op[0] * 9226291044724333206UL) + ((uint64_t)op[1] * 12070639868163715111UL) + ((uint64_t)op[2] * 1093162736320069068UL) + ((uint64_t)op[3] * 8887331613627784179UL) + ((uint64_t)op[4] * 8402347405076345566UL) + ((uint64_t)op[5] * 9664495871165975587UL) + ((uint64_t)op[6] * 5851854850093054879UL) + ((((uint64_t)op[7] * 1986626803087918625UL) + ((uint64_t)op[8] * 1170269630229199703UL) + ((uint64_t)op[9] * 15062026177312334876UL) + ((uint64_t)op[10] * 12555398568625425479UL)) * 4);
	tmp_q[7] = ((uint64_t)op[0] * 12555398568625425479UL) + ((uint64_t)op[1] * 9226291044724333206UL) + ((uint64_t)op[2] * 12070639868163715111UL) + ((uint64_t)op[3] * 1093162736320069068UL) + ((uint64_t)op[4] * 8887331613627784179UL) + ((uint64_t)op[5] * 8402347405076345566UL) + ((uint64_t)op[6] * 9664495871165975587UL) + ((uint64_t)op[7] * 5851854850093054879UL) + ((((uint64_t)op[8] * 1986626803087918625UL) + ((uint64_t)op[9] * 1170269630229199703UL) + ((uint64_t)op[10] * 15062026177312334876UL)) * 4);
	tmp_q[8] = ((uint64_t)op[0] * 15062026177312334876UL) + ((uint64_t)op[1] * 12555398568625425479UL) + ((uint64_t)op[2] * 9226291044724333206UL) + ((uint64_t)op[3] * 12070639868163715111UL) + ((uint64_t)op[4] * 1093162736320069068UL) + ((uint64_t)op[5] * 8887331613627784179UL) + ((uint64_t)op[6] * 8402347405076345566UL) + ((uint64_t)op[7] * 9664495871165975587UL) + ((uint64_t)op[8] * 5851854850093054879UL) + ((((uint64_t)op[9] * 1986626803087918625UL) + ((uint64_t)op[10] * 1170269630229199703UL)) * 4);
	tmp_q[9] = ((uint64_t)op[0] * 1170269630229199703UL) + ((uint64_t)op[1] * 15062026177312334876UL) + ((uint64_t)op[2] * 12555398568625425479UL) + ((uint64_t)op[3] * 9226291044724333206UL) + ((uint64_t)op[4] * 12070639868163715111UL) + ((uint64_t)op[5] * 1093162736320069068UL) + ((uint64_t)op[6] * 8887331613627784179UL) + ((uint64_t)op[7] * 8402347405076345566UL) + ((uint64_t)op[8] * 9664495871165975587UL) + ((uint64_t)op[9] * 5851854850093054879UL) + ((uint64_t)op[10] * 7946507212351674500UL);
	tmp_q[10] = ((uint64_t)op[0] * 1986626803087918625UL) + ((uint64_t)op[1] * 1170269630229199703UL) + ((uint64_t)op[2] * 15062026177312334876UL) + ((uint64_t)op[3] * 12555398568625425479UL) + ((uint64_t)op[4] * 9226291044724333206UL) + ((uint64_t)op[5] * 12070639868163715111UL) + ((uint64_t)op[6] * 1093162736320069068UL) + ((uint64_t)op[7] * 8887331613627784179UL) + ((uint64_t)op[8] * 8402347405076345566UL) + ((uint64_t)op[9] * 9664495871165975587UL) + ((uint64_t)op[10] * 5851854850093054879UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 64488342116119L) + ((-((int128)tmp_q[1] * 52260534817252L) - ((int128)tmp_q[2] * 20367165628150L) - ((int128)tmp_q[3] * 27779921859753L) + ((int128)tmp_q[4] * 29666677969291L) - ((int128)tmp_q[5] * 80567076598036L) + ((int128)tmp_q[6] * 12622655601359L) - ((int128)tmp_q[7] * 95265916203583L) - ((int128)tmp_q[8] * 43849324668678L) - ((int128)tmp_q[9] * 70495152739237L) - ((int128)tmp_q[10] * 6417196460721L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 6417196460721L) - ((int128)tmp_q[1] * 64488342116119L) + ((-((int128)tmp_q[2] * 52260534817252L) - ((int128)tmp_q[3] * 20367165628150L) - ((int128)tmp_q[4] * 27779921859753L) + ((int128)tmp_q[5] * 29666677969291L) - ((int128)tmp_q[6] * 80567076598036L) + ((int128)tmp_q[7] * 12622655601359L) - ((int128)tmp_q[8] * 95265916203583L) - ((int128)tmp_q[9] * 43849324668678L) - ((int128)tmp_q[10] * 70495152739237L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 70495152739237L) - ((int128)tmp_q[1] * 6417196460721L) - ((int128)tmp_q[2] * 64488342116119L) + ((-((int128)tmp_q[3] * 52260534817252L) - ((int128)tmp_q[4] * 20367165628150L) - ((int128)tmp_q[5] * 27779921859753L) + ((int128)tmp_q[6] * 29666677969291L) - ((int128)tmp_q[7] * 80567076598036L) + ((int128)tmp_q[8] * 12622655601359L) - ((int128)tmp_q[9] * 95265916203583L) - ((int128)tmp_q[10] * 43849324668678L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 43849324668678L) - ((int128)tmp_q[1] * 70495152739237L) - ((int128)tmp_q[2] * 6417196460721L) - ((int128)tmp_q[3] * 64488342116119L) + ((-((int128)tmp_q[4] * 52260534817252L) - ((int128)tmp_q[5] * 20367165628150L) - ((int128)tmp_q[6] * 27779921859753L) + ((int128)tmp_q[7] * 29666677969291L) - ((int128)tmp_q[8] * 80567076598036L) + ((int128)tmp_q[9] * 12622655601359L) - ((int128)tmp_q[10] * 95265916203583L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 95265916203583L) - ((int128)tmp_q[1] * 43849324668678L) - ((int128)tmp_q[2] * 70495152739237L) - ((int128)tmp_q[3] * 6417196460721L) - ((int128)tmp_q[4] * 64488342116119L) + ((-((int128)tmp_q[5] * 52260534817252L) - ((int128)tmp_q[6] * 20367165628150L) - ((int128)tmp_q[7] * 27779921859753L) + ((int128)tmp_q[8] * 29666677969291L) - ((int128)tmp_q[9] * 80567076598036L) + ((int128)tmp_q[10] * 12622655601359L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 12622655601359L) - ((int128)tmp_q[1] * 95265916203583L) - ((int128)tmp_q[2] * 43849324668678L) - ((int128)tmp_q[3] * 70495152739237L) - ((int128)tmp_q[4] * 6417196460721L) - ((int128)tmp_q[5] * 64488342116119L) + ((-((int128)tmp_q[6] * 52260534817252L) - ((int128)tmp_q[7] * 20367165628150L) - ((int128)tmp_q[8] * 27779921859753L) + ((int128)tmp_q[9] * 29666677969291L) - ((int128)tmp_q[10] * 80567076598036L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 80567076598036L) + ((int128)tmp_q[1] * 12622655601359L) - ((int128)tmp_q[2] * 95265916203583L) - ((int128)tmp_q[3] * 43849324668678L) - ((int128)tmp_q[4] * 70495152739237L) - ((int128)tmp_q[5] * 6417196460721L) - ((int128)tmp_q[6] * 64488342116119L) + ((-((int128)tmp_q[7] * 52260534817252L) - ((int128)tmp_q[8] * 20367165628150L) - ((int128)tmp_q[9] * 27779921859753L) + ((int128)tmp_q[10] * 29666677969291L)) * 4);
	tmp_zero[7] = ((int128)tmp_q[0] * 29666677969291L) - ((int128)tmp_q[1] * 80567076598036L) + ((int128)tmp_q[2] * 12622655601359L) - ((int128)tmp_q[3] * 95265916203583L) - ((int128)tmp_q[4] * 43849324668678L) - ((int128)tmp_q[5] * 70495152739237L) - ((int128)tmp_q[6] * 6417196460721L) - ((int128)tmp_q[7] * 64488342116119L) + ((-((int128)tmp_q[8] * 52260534817252L) - ((int128)tmp_q[9] * 20367165628150L) - ((int128)tmp_q[10] * 27779921859753L)) * 4);
	tmp_zero[8] = -((int128)tmp_q[0] * 27779921859753L) + ((int128)tmp_q[1] * 29666677969291L) - ((int128)tmp_q[2] * 80567076598036L) + ((int128)tmp_q[3] * 12622655601359L) - ((int128)tmp_q[4] * 95265916203583L) - ((int128)tmp_q[5] * 43849324668678L) - ((int128)tmp_q[6] * 70495152739237L) - ((int128)tmp_q[7] * 6417196460721L) - ((int128)tmp_q[8] * 64488342116119L) + ((-((int128)tmp_q[9] * 52260534817252L) - ((int128)tmp_q[10] * 20367165628150L)) * 4);
	tmp_zero[9] = -((int128)tmp_q[0] * 20367165628150L) - ((int128)tmp_q[1] * 27779921859753L) + ((int128)tmp_q[2] * 29666677969291L) - ((int128)tmp_q[3] * 80567076598036L) + ((int128)tmp_q[4] * 12622655601359L) - ((int128)tmp_q[5] * 95265916203583L) - ((int128)tmp_q[6] * 43849324668678L) - ((int128)tmp_q[7] * 70495152739237L) - ((int128)tmp_q[8] * 6417196460721L) - ((int128)tmp_q[9] * 64488342116119L) - ((int128)tmp_q[10] * 209042139269008L);
	tmp_zero[10] = -((int128)tmp_q[0] * 52260534817252L) - ((int128)tmp_q[1] * 20367165628150L) - ((int128)tmp_q[2] * 27779921859753L) + ((int128)tmp_q[3] * 29666677969291L) - ((int128)tmp_q[4] * 80567076598036L) + ((int128)tmp_q[5] * 12622655601359L) - ((int128)tmp_q[6] * 95265916203583L) - ((int128)tmp_q[7] * 43849324668678L) - ((int128)tmp_q[8] * 70495152739237L) - ((int128)tmp_q[9] * 6417196460721L) - ((int128)tmp_q[10] * 64488342116119L);

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

