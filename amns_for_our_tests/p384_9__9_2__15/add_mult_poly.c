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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[8] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8123737794134173605UL) + ((((uint64_t)op[1] * 8301129105893788212UL) + ((uint64_t)op[2] * 8556497955246888990UL) + ((uint64_t)op[3] * 14258841179466662863UL) + ((uint64_t)op[4] * 11783079497825577399UL) + ((uint64_t)op[5] * 15010541613232910131UL) + ((uint64_t)op[6] * 16614490805361746268UL) + ((uint64_t)op[7] * 211092645255721665UL) + ((uint64_t)op[8] * 4228332567585338785UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 4228332567585338785UL) + ((uint64_t)op[1] * 8123737794134173605UL) + ((((uint64_t)op[2] * 8301129105893788212UL) + ((uint64_t)op[3] * 8556497955246888990UL) + ((uint64_t)op[4] * 14258841179466662863UL) + ((uint64_t)op[5] * 11783079497825577399UL) + ((uint64_t)op[6] * 15010541613232910131UL) + ((uint64_t)op[7] * 16614490805361746268UL) + ((uint64_t)op[8] * 211092645255721665UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 211092645255721665UL) + ((uint64_t)op[1] * 4228332567585338785UL) + ((uint64_t)op[2] * 8123737794134173605UL) + ((((uint64_t)op[3] * 8301129105893788212UL) + ((uint64_t)op[4] * 8556497955246888990UL) + ((uint64_t)op[5] * 14258841179466662863UL) + ((uint64_t)op[6] * 11783079497825577399UL) + ((uint64_t)op[7] * 15010541613232910131UL) + ((uint64_t)op[8] * 16614490805361746268UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 16614490805361746268UL) + ((uint64_t)op[1] * 211092645255721665UL) + ((uint64_t)op[2] * 4228332567585338785UL) + ((uint64_t)op[3] * 8123737794134173605UL) + ((((uint64_t)op[4] * 8301129105893788212UL) + ((uint64_t)op[5] * 8556497955246888990UL) + ((uint64_t)op[6] * 14258841179466662863UL) + ((uint64_t)op[7] * 11783079497825577399UL) + ((uint64_t)op[8] * 15010541613232910131UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 15010541613232910131UL) + ((uint64_t)op[1] * 16614490805361746268UL) + ((uint64_t)op[2] * 211092645255721665UL) + ((uint64_t)op[3] * 4228332567585338785UL) + ((uint64_t)op[4] * 8123737794134173605UL) + ((((uint64_t)op[5] * 8301129105893788212UL) + ((uint64_t)op[6] * 8556497955246888990UL) + ((uint64_t)op[7] * 14258841179466662863UL) + ((uint64_t)op[8] * 11783079497825577399UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 11783079497825577399UL) + ((uint64_t)op[1] * 15010541613232910131UL) + ((uint64_t)op[2] * 16614490805361746268UL) + ((uint64_t)op[3] * 211092645255721665UL) + ((uint64_t)op[4] * 4228332567585338785UL) + ((uint64_t)op[5] * 8123737794134173605UL) + ((((uint64_t)op[6] * 8301129105893788212UL) + ((uint64_t)op[7] * 8556497955246888990UL) + ((uint64_t)op[8] * 14258841179466662863UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 14258841179466662863UL) + ((uint64_t)op[1] * 11783079497825577399UL) + ((uint64_t)op[2] * 15010541613232910131UL) + ((uint64_t)op[3] * 16614490805361746268UL) + ((uint64_t)op[4] * 211092645255721665UL) + ((uint64_t)op[5] * 4228332567585338785UL) + ((uint64_t)op[6] * 8123737794134173605UL) + ((((uint64_t)op[7] * 8301129105893788212UL) + ((uint64_t)op[8] * 8556497955246888990UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 8556497955246888990UL) + ((uint64_t)op[1] * 14258841179466662863UL) + ((uint64_t)op[2] * 11783079497825577399UL) + ((uint64_t)op[3] * 15010541613232910131UL) + ((uint64_t)op[4] * 16614490805361746268UL) + ((uint64_t)op[5] * 211092645255721665UL) + ((uint64_t)op[6] * 4228332567585338785UL) + ((uint64_t)op[7] * 8123737794134173605UL) + ((uint64_t)op[8] * 16602258211787576424UL);
	tmp_q[8] = ((uint64_t)op[0] * 8301129105893788212UL) + ((uint64_t)op[1] * 8556497955246888990UL) + ((uint64_t)op[2] * 14258841179466662863UL) + ((uint64_t)op[3] * 11783079497825577399UL) + ((uint64_t)op[4] * 15010541613232910131UL) + ((uint64_t)op[5] * 16614490805361746268UL) + ((uint64_t)op[6] * 211092645255721665UL) + ((uint64_t)op[7] * 4228332567585338785UL) + ((uint64_t)op[8] * 8123737794134173605UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4306067586171L) + ((-((int128)tmp_q[1] * 2110127236786L) + ((int128)tmp_q[2] * 1325518675802L) + ((int128)tmp_q[3] * 456876977857L) + ((int128)tmp_q[4] * 2847416302421L) - ((int128)tmp_q[5] * 2821045359078L) - ((int128)tmp_q[6] * 943012173847L) + ((int128)tmp_q[7] * 886690350596L) + ((int128)tmp_q[8] * 4321116898561L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 4321116898561L) + ((int128)tmp_q[1] * 4306067586171L) + ((-((int128)tmp_q[2] * 2110127236786L) + ((int128)tmp_q[3] * 1325518675802L) + ((int128)tmp_q[4] * 456876977857L) + ((int128)tmp_q[5] * 2847416302421L) - ((int128)tmp_q[6] * 2821045359078L) - ((int128)tmp_q[7] * 943012173847L) + ((int128)tmp_q[8] * 886690350596L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 886690350596L) + ((int128)tmp_q[1] * 4321116898561L) + ((int128)tmp_q[2] * 4306067586171L) + ((-((int128)tmp_q[3] * 2110127236786L) + ((int128)tmp_q[4] * 1325518675802L) + ((int128)tmp_q[5] * 456876977857L) + ((int128)tmp_q[6] * 2847416302421L) - ((int128)tmp_q[7] * 2821045359078L) - ((int128)tmp_q[8] * 943012173847L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 943012173847L) + ((int128)tmp_q[1] * 886690350596L) + ((int128)tmp_q[2] * 4321116898561L) + ((int128)tmp_q[3] * 4306067586171L) + ((-((int128)tmp_q[4] * 2110127236786L) + ((int128)tmp_q[5] * 1325518675802L) + ((int128)tmp_q[6] * 456876977857L) + ((int128)tmp_q[7] * 2847416302421L) - ((int128)tmp_q[8] * 2821045359078L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 2821045359078L) - ((int128)tmp_q[1] * 943012173847L) + ((int128)tmp_q[2] * 886690350596L) + ((int128)tmp_q[3] * 4321116898561L) + ((int128)tmp_q[4] * 4306067586171L) + ((-((int128)tmp_q[5] * 2110127236786L) + ((int128)tmp_q[6] * 1325518675802L) + ((int128)tmp_q[7] * 456876977857L) + ((int128)tmp_q[8] * 2847416302421L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 2847416302421L) - ((int128)tmp_q[1] * 2821045359078L) - ((int128)tmp_q[2] * 943012173847L) + ((int128)tmp_q[3] * 886690350596L) + ((int128)tmp_q[4] * 4321116898561L) + ((int128)tmp_q[5] * 4306067586171L) + ((-((int128)tmp_q[6] * 2110127236786L) + ((int128)tmp_q[7] * 1325518675802L) + ((int128)tmp_q[8] * 456876977857L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 456876977857L) + ((int128)tmp_q[1] * 2847416302421L) - ((int128)tmp_q[2] * 2821045359078L) - ((int128)tmp_q[3] * 943012173847L) + ((int128)tmp_q[4] * 886690350596L) + ((int128)tmp_q[5] * 4321116898561L) + ((int128)tmp_q[6] * 4306067586171L) + ((-((int128)tmp_q[7] * 2110127236786L) + ((int128)tmp_q[8] * 1325518675802L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 1325518675802L) + ((int128)tmp_q[1] * 456876977857L) + ((int128)tmp_q[2] * 2847416302421L) - ((int128)tmp_q[3] * 2821045359078L) - ((int128)tmp_q[4] * 943012173847L) + ((int128)tmp_q[5] * 886690350596L) + ((int128)tmp_q[6] * 4321116898561L) + ((int128)tmp_q[7] * 4306067586171L) - ((int128)tmp_q[8] * 4220254473572L);
	tmp_zero[8] = -((int128)tmp_q[0] * 2110127236786L) + ((int128)tmp_q[1] * 1325518675802L) + ((int128)tmp_q[2] * 456876977857L) + ((int128)tmp_q[3] * 2847416302421L) - ((int128)tmp_q[4] * 2821045359078L) - ((int128)tmp_q[5] * 943012173847L) + ((int128)tmp_q[6] * 886690350596L) + ((int128)tmp_q[7] * 4321116898561L) + ((int128)tmp_q[8] * 4306067586171L);

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
}

