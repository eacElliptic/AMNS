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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7968411653177305759UL) + ((((uint64_t)op[1] * 10717968107541554343UL) + ((uint64_t)op[2] * 14540438733435588375UL) + ((uint64_t)op[3] * 14779017751820133786UL) + ((uint64_t)op[4] * 1790764592689957032UL) + ((uint64_t)op[5] * 10658170695309691022UL) + ((uint64_t)op[6] * 4243034441394170293UL) + ((uint64_t)op[7] * 6538614096098137565UL) + ((uint64_t)op[8] * 10958957112638522151UL) + ((uint64_t)op[9] * 7661202633570664434UL) + ((uint64_t)op[10] * 14428301483122420977UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 14428301483122420977UL) + ((uint64_t)op[1] * 7968411653177305759UL) + ((((uint64_t)op[2] * 10717968107541554343UL) + ((uint64_t)op[3] * 14540438733435588375UL) + ((uint64_t)op[4] * 14779017751820133786UL) + ((uint64_t)op[5] * 1790764592689957032UL) + ((uint64_t)op[6] * 10658170695309691022UL) + ((uint64_t)op[7] * 4243034441394170293UL) + ((uint64_t)op[8] * 6538614096098137565UL) + ((uint64_t)op[9] * 10958957112638522151UL) + ((uint64_t)op[10] * 7661202633570664434UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 7661202633570664434UL) + ((uint64_t)op[1] * 14428301483122420977UL) + ((uint64_t)op[2] * 7968411653177305759UL) + ((((uint64_t)op[3] * 10717968107541554343UL) + ((uint64_t)op[4] * 14540438733435588375UL) + ((uint64_t)op[5] * 14779017751820133786UL) + ((uint64_t)op[6] * 1790764592689957032UL) + ((uint64_t)op[7] * 10658170695309691022UL) + ((uint64_t)op[8] * 4243034441394170293UL) + ((uint64_t)op[9] * 6538614096098137565UL) + ((uint64_t)op[10] * 10958957112638522151UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 10958957112638522151UL) + ((uint64_t)op[1] * 7661202633570664434UL) + ((uint64_t)op[2] * 14428301483122420977UL) + ((uint64_t)op[3] * 7968411653177305759UL) + ((((uint64_t)op[4] * 10717968107541554343UL) + ((uint64_t)op[5] * 14540438733435588375UL) + ((uint64_t)op[6] * 14779017751820133786UL) + ((uint64_t)op[7] * 1790764592689957032UL) + ((uint64_t)op[8] * 10658170695309691022UL) + ((uint64_t)op[9] * 4243034441394170293UL) + ((uint64_t)op[10] * 6538614096098137565UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 6538614096098137565UL) + ((uint64_t)op[1] * 10958957112638522151UL) + ((uint64_t)op[2] * 7661202633570664434UL) + ((uint64_t)op[3] * 14428301483122420977UL) + ((uint64_t)op[4] * 7968411653177305759UL) + ((((uint64_t)op[5] * 10717968107541554343UL) + ((uint64_t)op[6] * 14540438733435588375UL) + ((uint64_t)op[7] * 14779017751820133786UL) + ((uint64_t)op[8] * 1790764592689957032UL) + ((uint64_t)op[9] * 10658170695309691022UL) + ((uint64_t)op[10] * 4243034441394170293UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 4243034441394170293UL) + ((uint64_t)op[1] * 6538614096098137565UL) + ((uint64_t)op[2] * 10958957112638522151UL) + ((uint64_t)op[3] * 7661202633570664434UL) + ((uint64_t)op[4] * 14428301483122420977UL) + ((uint64_t)op[5] * 7968411653177305759UL) + ((((uint64_t)op[6] * 10717968107541554343UL) + ((uint64_t)op[7] * 14540438733435588375UL) + ((uint64_t)op[8] * 14779017751820133786UL) + ((uint64_t)op[9] * 1790764592689957032UL) + ((uint64_t)op[10] * 10658170695309691022UL)) * 18446744073709551612);
	tmp_q[6] = ((uint64_t)op[0] * 10658170695309691022UL) + ((uint64_t)op[1] * 4243034441394170293UL) + ((uint64_t)op[2] * 6538614096098137565UL) + ((uint64_t)op[3] * 10958957112638522151UL) + ((uint64_t)op[4] * 7661202633570664434UL) + ((uint64_t)op[5] * 14428301483122420977UL) + ((uint64_t)op[6] * 7968411653177305759UL) + ((((uint64_t)op[7] * 10717968107541554343UL) + ((uint64_t)op[8] * 14540438733435588375UL) + ((uint64_t)op[9] * 14779017751820133786UL) + ((uint64_t)op[10] * 1790764592689957032UL)) * 18446744073709551612);
	tmp_q[7] = ((uint64_t)op[0] * 1790764592689957032UL) + ((uint64_t)op[1] * 10658170695309691022UL) + ((uint64_t)op[2] * 4243034441394170293UL) + ((uint64_t)op[3] * 6538614096098137565UL) + ((uint64_t)op[4] * 10958957112638522151UL) + ((uint64_t)op[5] * 7661202633570664434UL) + ((uint64_t)op[6] * 14428301483122420977UL) + ((uint64_t)op[7] * 7968411653177305759UL) + ((((uint64_t)op[8] * 10717968107541554343UL) + ((uint64_t)op[9] * 14540438733435588375UL) + ((uint64_t)op[10] * 14779017751820133786UL)) * 18446744073709551612);
	tmp_q[8] = ((uint64_t)op[0] * 14779017751820133786UL) + ((uint64_t)op[1] * 1790764592689957032UL) + ((uint64_t)op[2] * 10658170695309691022UL) + ((uint64_t)op[3] * 4243034441394170293UL) + ((uint64_t)op[4] * 6538614096098137565UL) + ((uint64_t)op[5] * 10958957112638522151UL) + ((uint64_t)op[6] * 7661202633570664434UL) + ((uint64_t)op[7] * 14428301483122420977UL) + ((uint64_t)op[8] * 7968411653177305759UL) + ((((uint64_t)op[9] * 10717968107541554343UL) + ((uint64_t)op[10] * 14540438733435588375UL)) * 18446744073709551612);
	tmp_q[9] = ((uint64_t)op[0] * 14540438733435588375UL) + ((uint64_t)op[1] * 14779017751820133786UL) + ((uint64_t)op[2] * 1790764592689957032UL) + ((uint64_t)op[3] * 10658170695309691022UL) + ((uint64_t)op[4] * 4243034441394170293UL) + ((uint64_t)op[5] * 6538614096098137565UL) + ((uint64_t)op[6] * 10958957112638522151UL) + ((uint64_t)op[7] * 7661202633570664434UL) + ((uint64_t)op[8] * 14428301483122420977UL) + ((uint64_t)op[9] * 7968411653177305759UL) + ((uint64_t)op[10] * 12468359790962437476UL);
	tmp_q[10] = ((uint64_t)op[0] * 10717968107541554343UL) + ((uint64_t)op[1] * 14540438733435588375UL) + ((uint64_t)op[2] * 14779017751820133786UL) + ((uint64_t)op[3] * 1790764592689957032UL) + ((uint64_t)op[4] * 10658170695309691022UL) + ((uint64_t)op[5] * 4243034441394170293UL) + ((uint64_t)op[6] * 6538614096098137565UL) + ((uint64_t)op[7] * 10958957112638522151UL) + ((uint64_t)op[8] * 7661202633570664434UL) + ((uint64_t)op[9] * 14428301483122420977UL) + ((uint64_t)op[10] * 7968411653177305759UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 90254424573361L) - ((-((int128)tmp_q[1] * 29052543339198L) + ((int128)tmp_q[2] * 61007061716404L) + ((int128)tmp_q[3] * 39737056303755L) + ((int128)tmp_q[4] * 47564874354292L) - ((int128)tmp_q[5] * 89216675495257L) - ((int128)tmp_q[6] * 25165419269123L) + ((int128)tmp_q[7] * 22688543605302L) - ((int128)tmp_q[8] * 34262172809020L) - ((int128)tmp_q[9] * 25549293365369L) - ((int128)tmp_q[10] * 54778546518871L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 54778546518871L) + ((int128)tmp_q[1] * 90254424573361L) - ((-((int128)tmp_q[2] * 29052543339198L) + ((int128)tmp_q[3] * 61007061716404L) + ((int128)tmp_q[4] * 39737056303755L) + ((int128)tmp_q[5] * 47564874354292L) - ((int128)tmp_q[6] * 89216675495257L) - ((int128)tmp_q[7] * 25165419269123L) + ((int128)tmp_q[8] * 22688543605302L) - ((int128)tmp_q[9] * 34262172809020L) - ((int128)tmp_q[10] * 25549293365369L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 25549293365369L) - ((int128)tmp_q[1] * 54778546518871L) + ((int128)tmp_q[2] * 90254424573361L) - ((-((int128)tmp_q[3] * 29052543339198L) + ((int128)tmp_q[4] * 61007061716404L) + ((int128)tmp_q[5] * 39737056303755L) + ((int128)tmp_q[6] * 47564874354292L) - ((int128)tmp_q[7] * 89216675495257L) - ((int128)tmp_q[8] * 25165419269123L) + ((int128)tmp_q[9] * 22688543605302L) - ((int128)tmp_q[10] * 34262172809020L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 34262172809020L) - ((int128)tmp_q[1] * 25549293365369L) - ((int128)tmp_q[2] * 54778546518871L) + ((int128)tmp_q[3] * 90254424573361L) - ((-((int128)tmp_q[4] * 29052543339198L) + ((int128)tmp_q[5] * 61007061716404L) + ((int128)tmp_q[6] * 39737056303755L) + ((int128)tmp_q[7] * 47564874354292L) - ((int128)tmp_q[8] * 89216675495257L) - ((int128)tmp_q[9] * 25165419269123L) + ((int128)tmp_q[10] * 22688543605302L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 22688543605302L) - ((int128)tmp_q[1] * 34262172809020L) - ((int128)tmp_q[2] * 25549293365369L) - ((int128)tmp_q[3] * 54778546518871L) + ((int128)tmp_q[4] * 90254424573361L) - ((-((int128)tmp_q[5] * 29052543339198L) + ((int128)tmp_q[6] * 61007061716404L) + ((int128)tmp_q[7] * 39737056303755L) + ((int128)tmp_q[8] * 47564874354292L) - ((int128)tmp_q[9] * 89216675495257L) - ((int128)tmp_q[10] * 25165419269123L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 25165419269123L) + ((int128)tmp_q[1] * 22688543605302L) - ((int128)tmp_q[2] * 34262172809020L) - ((int128)tmp_q[3] * 25549293365369L) - ((int128)tmp_q[4] * 54778546518871L) + ((int128)tmp_q[5] * 90254424573361L) - ((-((int128)tmp_q[6] * 29052543339198L) + ((int128)tmp_q[7] * 61007061716404L) + ((int128)tmp_q[8] * 39737056303755L) + ((int128)tmp_q[9] * 47564874354292L) - ((int128)tmp_q[10] * 89216675495257L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 89216675495257L) - ((int128)tmp_q[1] * 25165419269123L) + ((int128)tmp_q[2] * 22688543605302L) - ((int128)tmp_q[3] * 34262172809020L) - ((int128)tmp_q[4] * 25549293365369L) - ((int128)tmp_q[5] * 54778546518871L) + ((int128)tmp_q[6] * 90254424573361L) - ((-((int128)tmp_q[7] * 29052543339198L) + ((int128)tmp_q[8] * 61007061716404L) + ((int128)tmp_q[9] * 39737056303755L) + ((int128)tmp_q[10] * 47564874354292L)) * 4);
	tmp_zero[7] = ((int128)tmp_q[0] * 47564874354292L) - ((int128)tmp_q[1] * 89216675495257L) - ((int128)tmp_q[2] * 25165419269123L) + ((int128)tmp_q[3] * 22688543605302L) - ((int128)tmp_q[4] * 34262172809020L) - ((int128)tmp_q[5] * 25549293365369L) - ((int128)tmp_q[6] * 54778546518871L) + ((int128)tmp_q[7] * 90254424573361L) - ((-((int128)tmp_q[8] * 29052543339198L) + ((int128)tmp_q[9] * 61007061716404L) + ((int128)tmp_q[10] * 39737056303755L)) * 4);
	tmp_zero[8] = ((int128)tmp_q[0] * 39737056303755L) + ((int128)tmp_q[1] * 47564874354292L) - ((int128)tmp_q[2] * 89216675495257L) - ((int128)tmp_q[3] * 25165419269123L) + ((int128)tmp_q[4] * 22688543605302L) - ((int128)tmp_q[5] * 34262172809020L) - ((int128)tmp_q[6] * 25549293365369L) - ((int128)tmp_q[7] * 54778546518871L) + ((int128)tmp_q[8] * 90254424573361L) - ((-((int128)tmp_q[9] * 29052543339198L) + ((int128)tmp_q[10] * 61007061716404L)) * 4);
	tmp_zero[9] = ((int128)tmp_q[0] * 61007061716404L) + ((int128)tmp_q[1] * 39737056303755L) + ((int128)tmp_q[2] * 47564874354292L) - ((int128)tmp_q[3] * 89216675495257L) - ((int128)tmp_q[4] * 25165419269123L) + ((int128)tmp_q[5] * 22688543605302L) - ((int128)tmp_q[6] * 34262172809020L) - ((int128)tmp_q[7] * 25549293365369L) - ((int128)tmp_q[8] * 54778546518871L) + ((int128)tmp_q[9] * 90254424573361L) + ((int128)tmp_q[10] * 116210173356792L);
	tmp_zero[10] = -((int128)tmp_q[0] * 29052543339198L) + ((int128)tmp_q[1] * 61007061716404L) + ((int128)tmp_q[2] * 39737056303755L) + ((int128)tmp_q[3] * 47564874354292L) - ((int128)tmp_q[4] * 89216675495257L) - ((int128)tmp_q[5] * 25165419269123L) + ((int128)tmp_q[6] * 22688543605302L) - ((int128)tmp_q[7] * 34262172809020L) - ((int128)tmp_q[8] * 25549293365369L) - ((int128)tmp_q[9] * 54778546518871L) + ((int128)tmp_q[10] * 90254424573361L);

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

