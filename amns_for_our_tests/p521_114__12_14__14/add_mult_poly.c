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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11036042466249130213UL) + ((((uint64_t)op[1] * 16717715905366450109UL) + ((uint64_t)op[2] * 7509539440426322358UL) + ((uint64_t)op[3] * 13263611064373458816UL) + ((uint64_t)op[4] * 1280029432987368483UL) + ((uint64_t)op[5] * 4348454407688982625UL) + ((uint64_t)op[6] * 1201602343619193026UL) + ((uint64_t)op[7] * 7256010295119599817UL) + ((uint64_t)op[8] * 3775972711501086697UL) + ((uint64_t)op[9] * 17049971660829750810UL) + ((uint64_t)op[10] * 189323481733675003UL) + ((uint64_t)op[11] * 9762236051368343767UL)) * 14);
	tmp_q[1] = ((uint64_t)op[0] * 9762236051368343767UL) + ((uint64_t)op[1] * 11036042466249130213UL) + ((((uint64_t)op[2] * 16717715905366450109UL) + ((uint64_t)op[3] * 7509539440426322358UL) + ((uint64_t)op[4] * 13263611064373458816UL) + ((uint64_t)op[5] * 1280029432987368483UL) + ((uint64_t)op[6] * 4348454407688982625UL) + ((uint64_t)op[7] * 1201602343619193026UL) + ((uint64_t)op[8] * 7256010295119599817UL) + ((uint64_t)op[9] * 3775972711501086697UL) + ((uint64_t)op[10] * 17049971660829750810UL) + ((uint64_t)op[11] * 189323481733675003UL)) * 14);
	tmp_q[2] = ((uint64_t)op[0] * 189323481733675003UL) + ((uint64_t)op[1] * 9762236051368343767UL) + ((uint64_t)op[2] * 11036042466249130213UL) + ((((uint64_t)op[3] * 16717715905366450109UL) + ((uint64_t)op[4] * 7509539440426322358UL) + ((uint64_t)op[5] * 13263611064373458816UL) + ((uint64_t)op[6] * 1280029432987368483UL) + ((uint64_t)op[7] * 4348454407688982625UL) + ((uint64_t)op[8] * 1201602343619193026UL) + ((uint64_t)op[9] * 7256010295119599817UL) + ((uint64_t)op[10] * 3775972711501086697UL) + ((uint64_t)op[11] * 17049971660829750810UL)) * 14);
	tmp_q[3] = ((uint64_t)op[0] * 17049971660829750810UL) + ((uint64_t)op[1] * 189323481733675003UL) + ((uint64_t)op[2] * 9762236051368343767UL) + ((uint64_t)op[3] * 11036042466249130213UL) + ((((uint64_t)op[4] * 16717715905366450109UL) + ((uint64_t)op[5] * 7509539440426322358UL) + ((uint64_t)op[6] * 13263611064373458816UL) + ((uint64_t)op[7] * 1280029432987368483UL) + ((uint64_t)op[8] * 4348454407688982625UL) + ((uint64_t)op[9] * 1201602343619193026UL) + ((uint64_t)op[10] * 7256010295119599817UL) + ((uint64_t)op[11] * 3775972711501086697UL)) * 14);
	tmp_q[4] = ((uint64_t)op[0] * 3775972711501086697UL) + ((uint64_t)op[1] * 17049971660829750810UL) + ((uint64_t)op[2] * 189323481733675003UL) + ((uint64_t)op[3] * 9762236051368343767UL) + ((uint64_t)op[4] * 11036042466249130213UL) + ((((uint64_t)op[5] * 16717715905366450109UL) + ((uint64_t)op[6] * 7509539440426322358UL) + ((uint64_t)op[7] * 13263611064373458816UL) + ((uint64_t)op[8] * 1280029432987368483UL) + ((uint64_t)op[9] * 4348454407688982625UL) + ((uint64_t)op[10] * 1201602343619193026UL) + ((uint64_t)op[11] * 7256010295119599817UL)) * 14);
	tmp_q[5] = ((uint64_t)op[0] * 7256010295119599817UL) + ((uint64_t)op[1] * 3775972711501086697UL) + ((uint64_t)op[2] * 17049971660829750810UL) + ((uint64_t)op[3] * 189323481733675003UL) + ((uint64_t)op[4] * 9762236051368343767UL) + ((uint64_t)op[5] * 11036042466249130213UL) + ((((uint64_t)op[6] * 16717715905366450109UL) + ((uint64_t)op[7] * 7509539440426322358UL) + ((uint64_t)op[8] * 13263611064373458816UL) + ((uint64_t)op[9] * 1280029432987368483UL) + ((uint64_t)op[10] * 4348454407688982625UL) + ((uint64_t)op[11] * 1201602343619193026UL)) * 14);
	tmp_q[6] = ((uint64_t)op[0] * 1201602343619193026UL) + ((uint64_t)op[1] * 7256010295119599817UL) + ((uint64_t)op[2] * 3775972711501086697UL) + ((uint64_t)op[3] * 17049971660829750810UL) + ((uint64_t)op[4] * 189323481733675003UL) + ((uint64_t)op[5] * 9762236051368343767UL) + ((uint64_t)op[6] * 11036042466249130213UL) + ((((uint64_t)op[7] * 16717715905366450109UL) + ((uint64_t)op[8] * 7509539440426322358UL) + ((uint64_t)op[9] * 13263611064373458816UL) + ((uint64_t)op[10] * 1280029432987368483UL) + ((uint64_t)op[11] * 4348454407688982625UL)) * 14);
	tmp_q[7] = ((uint64_t)op[0] * 4348454407688982625UL) + ((uint64_t)op[1] * 1201602343619193026UL) + ((uint64_t)op[2] * 7256010295119599817UL) + ((uint64_t)op[3] * 3775972711501086697UL) + ((uint64_t)op[4] * 17049971660829750810UL) + ((uint64_t)op[5] * 189323481733675003UL) + ((uint64_t)op[6] * 9762236051368343767UL) + ((uint64_t)op[7] * 11036042466249130213UL) + ((((uint64_t)op[8] * 16717715905366450109UL) + ((uint64_t)op[9] * 7509539440426322358UL) + ((uint64_t)op[10] * 13263611064373458816UL) + ((uint64_t)op[11] * 1280029432987368483UL)) * 14);
	tmp_q[8] = ((uint64_t)op[0] * 1280029432987368483UL) + ((uint64_t)op[1] * 4348454407688982625UL) + ((uint64_t)op[2] * 1201602343619193026UL) + ((uint64_t)op[3] * 7256010295119599817UL) + ((uint64_t)op[4] * 3775972711501086697UL) + ((uint64_t)op[5] * 17049971660829750810UL) + ((uint64_t)op[6] * 189323481733675003UL) + ((uint64_t)op[7] * 9762236051368343767UL) + ((uint64_t)op[8] * 11036042466249130213UL) + ((((uint64_t)op[9] * 16717715905366450109UL) + ((uint64_t)op[10] * 7509539440426322358UL) + ((uint64_t)op[11] * 13263611064373458816UL)) * 14);
	tmp_q[9] = ((uint64_t)op[0] * 13263611064373458816UL) + ((uint64_t)op[1] * 1280029432987368483UL) + ((uint64_t)op[2] * 4348454407688982625UL) + ((uint64_t)op[3] * 1201602343619193026UL) + ((uint64_t)op[4] * 7256010295119599817UL) + ((uint64_t)op[5] * 3775972711501086697UL) + ((uint64_t)op[6] * 17049971660829750810UL) + ((uint64_t)op[7] * 189323481733675003UL) + ((uint64_t)op[8] * 9762236051368343767UL) + ((uint64_t)op[9] * 11036042466249130213UL) + ((((uint64_t)op[10] * 16717715905366450109UL) + ((uint64_t)op[11] * 7509539440426322358UL)) * 14);
	tmp_q[10] = ((uint64_t)op[0] * 7509539440426322358UL) + ((uint64_t)op[1] * 13263611064373458816UL) + ((uint64_t)op[2] * 1280029432987368483UL) + ((uint64_t)op[3] * 4348454407688982625UL) + ((uint64_t)op[4] * 1201602343619193026UL) + ((uint64_t)op[5] * 7256010295119599817UL) + ((uint64_t)op[6] * 3775972711501086697UL) + ((uint64_t)op[7] * 17049971660829750810UL) + ((uint64_t)op[8] * 189323481733675003UL) + ((uint64_t)op[9] * 9762236051368343767UL) + ((uint64_t)op[10] * 11036042466249130213UL) + ((uint64_t)op[11] * 12687093790615682134UL);
	tmp_q[11] = ((uint64_t)op[0] * 16717715905366450109UL) + ((uint64_t)op[1] * 7509539440426322358UL) + ((uint64_t)op[2] * 13263611064373458816UL) + ((uint64_t)op[3] * 1280029432987368483UL) + ((uint64_t)op[4] * 4348454407688982625UL) + ((uint64_t)op[5] * 1201602343619193026UL) + ((uint64_t)op[6] * 7256010295119599817UL) + ((uint64_t)op[7] * 3775972711501086697UL) + ((uint64_t)op[8] * 17049971660829750810UL) + ((uint64_t)op[9] * 189323481733675003UL) + ((uint64_t)op[10] * 9762236051368343767UL) + ((uint64_t)op[11] * 11036042466249130213UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2865464289689L) + ((-((int128)tmp_q[1] * 2334629458781L) - ((int128)tmp_q[2] * 504042314465L) - ((int128)tmp_q[3] * 3024994324511L) + ((int128)tmp_q[4] * 4325420975166L) - ((int128)tmp_q[5] * 57768311181L) - ((int128)tmp_q[6] * 6299679492494L) - ((int128)tmp_q[7] * 2992185645775L) + ((int128)tmp_q[8] * 1950565702210L) + ((int128)tmp_q[9] * 5354231837075L) - ((int128)tmp_q[10] * 2598956588428L) + ((int128)tmp_q[11] * 1451346432163L)) * 14);
	tmp_zero[1] = ((int128)tmp_q[0] * 1451346432163L) + ((int128)tmp_q[1] * 2865464289689L) + ((-((int128)tmp_q[2] * 2334629458781L) - ((int128)tmp_q[3] * 504042314465L) - ((int128)tmp_q[4] * 3024994324511L) + ((int128)tmp_q[5] * 4325420975166L) - ((int128)tmp_q[6] * 57768311181L) - ((int128)tmp_q[7] * 6299679492494L) - ((int128)tmp_q[8] * 2992185645775L) + ((int128)tmp_q[9] * 1950565702210L) + ((int128)tmp_q[10] * 5354231837075L) - ((int128)tmp_q[11] * 2598956588428L)) * 14);
	tmp_zero[2] = -((int128)tmp_q[0] * 2598956588428L) + ((int128)tmp_q[1] * 1451346432163L) + ((int128)tmp_q[2] * 2865464289689L) + ((-((int128)tmp_q[3] * 2334629458781L) - ((int128)tmp_q[4] * 504042314465L) - ((int128)tmp_q[5] * 3024994324511L) + ((int128)tmp_q[6] * 4325420975166L) - ((int128)tmp_q[7] * 57768311181L) - ((int128)tmp_q[8] * 6299679492494L) - ((int128)tmp_q[9] * 2992185645775L) + ((int128)tmp_q[10] * 1950565702210L) + ((int128)tmp_q[11] * 5354231837075L)) * 14);
	tmp_zero[3] = ((int128)tmp_q[0] * 5354231837075L) - ((int128)tmp_q[1] * 2598956588428L) + ((int128)tmp_q[2] * 1451346432163L) + ((int128)tmp_q[3] * 2865464289689L) + ((-((int128)tmp_q[4] * 2334629458781L) - ((int128)tmp_q[5] * 504042314465L) - ((int128)tmp_q[6] * 3024994324511L) + ((int128)tmp_q[7] * 4325420975166L) - ((int128)tmp_q[8] * 57768311181L) - ((int128)tmp_q[9] * 6299679492494L) - ((int128)tmp_q[10] * 2992185645775L) + ((int128)tmp_q[11] * 1950565702210L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 1950565702210L) + ((int128)tmp_q[1] * 5354231837075L) - ((int128)tmp_q[2] * 2598956588428L) + ((int128)tmp_q[3] * 1451346432163L) + ((int128)tmp_q[4] * 2865464289689L) + ((-((int128)tmp_q[5] * 2334629458781L) - ((int128)tmp_q[6] * 504042314465L) - ((int128)tmp_q[7] * 3024994324511L) + ((int128)tmp_q[8] * 4325420975166L) - ((int128)tmp_q[9] * 57768311181L) - ((int128)tmp_q[10] * 6299679492494L) - ((int128)tmp_q[11] * 2992185645775L)) * 14);
	tmp_zero[5] = -((int128)tmp_q[0] * 2992185645775L) + ((int128)tmp_q[1] * 1950565702210L) + ((int128)tmp_q[2] * 5354231837075L) - ((int128)tmp_q[3] * 2598956588428L) + ((int128)tmp_q[4] * 1451346432163L) + ((int128)tmp_q[5] * 2865464289689L) + ((-((int128)tmp_q[6] * 2334629458781L) - ((int128)tmp_q[7] * 504042314465L) - ((int128)tmp_q[8] * 3024994324511L) + ((int128)tmp_q[9] * 4325420975166L) - ((int128)tmp_q[10] * 57768311181L) - ((int128)tmp_q[11] * 6299679492494L)) * 14);
	tmp_zero[6] = -((int128)tmp_q[0] * 6299679492494L) - ((int128)tmp_q[1] * 2992185645775L) + ((int128)tmp_q[2] * 1950565702210L) + ((int128)tmp_q[3] * 5354231837075L) - ((int128)tmp_q[4] * 2598956588428L) + ((int128)tmp_q[5] * 1451346432163L) + ((int128)tmp_q[6] * 2865464289689L) + ((-((int128)tmp_q[7] * 2334629458781L) - ((int128)tmp_q[8] * 504042314465L) - ((int128)tmp_q[9] * 3024994324511L) + ((int128)tmp_q[10] * 4325420975166L) - ((int128)tmp_q[11] * 57768311181L)) * 14);
	tmp_zero[7] = -((int128)tmp_q[0] * 57768311181L) - ((int128)tmp_q[1] * 6299679492494L) - ((int128)tmp_q[2] * 2992185645775L) + ((int128)tmp_q[3] * 1950565702210L) + ((int128)tmp_q[4] * 5354231837075L) - ((int128)tmp_q[5] * 2598956588428L) + ((int128)tmp_q[6] * 1451346432163L) + ((int128)tmp_q[7] * 2865464289689L) + ((-((int128)tmp_q[8] * 2334629458781L) - ((int128)tmp_q[9] * 504042314465L) - ((int128)tmp_q[10] * 3024994324511L) + ((int128)tmp_q[11] * 4325420975166L)) * 14);
	tmp_zero[8] = ((int128)tmp_q[0] * 4325420975166L) - ((int128)tmp_q[1] * 57768311181L) - ((int128)tmp_q[2] * 6299679492494L) - ((int128)tmp_q[3] * 2992185645775L) + ((int128)tmp_q[4] * 1950565702210L) + ((int128)tmp_q[5] * 5354231837075L) - ((int128)tmp_q[6] * 2598956588428L) + ((int128)tmp_q[7] * 1451346432163L) + ((int128)tmp_q[8] * 2865464289689L) + ((-((int128)tmp_q[9] * 2334629458781L) - ((int128)tmp_q[10] * 504042314465L) - ((int128)tmp_q[11] * 3024994324511L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 3024994324511L) + ((int128)tmp_q[1] * 4325420975166L) - ((int128)tmp_q[2] * 57768311181L) - ((int128)tmp_q[3] * 6299679492494L) - ((int128)tmp_q[4] * 2992185645775L) + ((int128)tmp_q[5] * 1950565702210L) + ((int128)tmp_q[6] * 5354231837075L) - ((int128)tmp_q[7] * 2598956588428L) + ((int128)tmp_q[8] * 1451346432163L) + ((int128)tmp_q[9] * 2865464289689L) + ((-((int128)tmp_q[10] * 2334629458781L) - ((int128)tmp_q[11] * 504042314465L)) * 14);
	tmp_zero[10] = -((int128)tmp_q[0] * 504042314465L) - ((int128)tmp_q[1] * 3024994324511L) + ((int128)tmp_q[2] * 4325420975166L) - ((int128)tmp_q[3] * 57768311181L) - ((int128)tmp_q[4] * 6299679492494L) - ((int128)tmp_q[5] * 2992185645775L) + ((int128)tmp_q[6] * 1950565702210L) + ((int128)tmp_q[7] * 5354231837075L) - ((int128)tmp_q[8] * 2598956588428L) + ((int128)tmp_q[9] * 1451346432163L) + ((int128)tmp_q[10] * 2865464289689L) - ((int128)tmp_q[11] * 32684812422934L);
	tmp_zero[11] = -((int128)tmp_q[0] * 2334629458781L) - ((int128)tmp_q[1] * 504042314465L) - ((int128)tmp_q[2] * 3024994324511L) + ((int128)tmp_q[3] * 4325420975166L) - ((int128)tmp_q[4] * 57768311181L) - ((int128)tmp_q[5] * 6299679492494L) - ((int128)tmp_q[6] * 2992185645775L) + ((int128)tmp_q[7] * 1950565702210L) + ((int128)tmp_q[8] * 5354231837075L) - ((int128)tmp_q[9] * 2598956588428L) + ((int128)tmp_q[10] * 1451346432163L) + ((int128)tmp_q[11] * 2865464289689L);

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

