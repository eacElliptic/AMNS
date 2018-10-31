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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14978621653772995207UL) + ((((uint64_t)op[1] * 1876691394638825326UL) + ((uint64_t)op[2] * 17661544290795815318UL) + ((uint64_t)op[3] * 4937479850530143461UL) + ((uint64_t)op[4] * 3559909662137972224UL) + ((uint64_t)op[5] * 7788521772015481236UL) + ((uint64_t)op[6] * 17226400798528396567UL) + ((uint64_t)op[7] * 12455287797648325778UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 12455287797648325778UL) + ((uint64_t)op[1] * 14978621653772995207UL) + ((((uint64_t)op[2] * 1876691394638825326UL) + ((uint64_t)op[3] * 17661544290795815318UL) + ((uint64_t)op[4] * 4937479850530143461UL) + ((uint64_t)op[5] * 3559909662137972224UL) + ((uint64_t)op[6] * 7788521772015481236UL) + ((uint64_t)op[7] * 17226400798528396567UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 17226400798528396567UL) + ((uint64_t)op[1] * 12455287797648325778UL) + ((uint64_t)op[2] * 14978621653772995207UL) + ((((uint64_t)op[3] * 1876691394638825326UL) + ((uint64_t)op[4] * 17661544290795815318UL) + ((uint64_t)op[5] * 4937479850530143461UL) + ((uint64_t)op[6] * 3559909662137972224UL) + ((uint64_t)op[7] * 7788521772015481236UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 7788521772015481236UL) + ((uint64_t)op[1] * 17226400798528396567UL) + ((uint64_t)op[2] * 12455287797648325778UL) + ((uint64_t)op[3] * 14978621653772995207UL) + ((((uint64_t)op[4] * 1876691394638825326UL) + ((uint64_t)op[5] * 17661544290795815318UL) + ((uint64_t)op[6] * 4937479850530143461UL) + ((uint64_t)op[7] * 3559909662137972224UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 3559909662137972224UL) + ((uint64_t)op[1] * 7788521772015481236UL) + ((uint64_t)op[2] * 17226400798528396567UL) + ((uint64_t)op[3] * 12455287797648325778UL) + ((uint64_t)op[4] * 14978621653772995207UL) + ((((uint64_t)op[5] * 1876691394638825326UL) + ((uint64_t)op[6] * 17661544290795815318UL) + ((uint64_t)op[7] * 4937479850530143461UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 4937479850530143461UL) + ((uint64_t)op[1] * 3559909662137972224UL) + ((uint64_t)op[2] * 7788521772015481236UL) + ((uint64_t)op[3] * 17226400798528396567UL) + ((uint64_t)op[4] * 12455287797648325778UL) + ((uint64_t)op[5] * 14978621653772995207UL) + ((((uint64_t)op[6] * 1876691394638825326UL) + ((uint64_t)op[7] * 17661544290795815318UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 17661544290795815318UL) + ((uint64_t)op[1] * 4937479850530143461UL) + ((uint64_t)op[2] * 3559909662137972224UL) + ((uint64_t)op[3] * 7788521772015481236UL) + ((uint64_t)op[4] * 17226400798528396567UL) + ((uint64_t)op[5] * 12455287797648325778UL) + ((uint64_t)op[6] * 14978621653772995207UL) + ((uint64_t)op[7] * 13136839762471777282UL);
	tmp_q[7] = ((uint64_t)op[0] * 1876691394638825326UL) + ((uint64_t)op[1] * 17661544290795815318UL) + ((uint64_t)op[2] * 4937479850530143461UL) + ((uint64_t)op[3] * 3559909662137972224UL) + ((uint64_t)op[4] * 7788521772015481236UL) + ((uint64_t)op[5] * 17226400798528396567UL) + ((uint64_t)op[6] * 12455287797648325778UL) + ((uint64_t)op[7] * 14978621653772995207UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 46377049270410L) + ((-((int128)tmp_q[1] * 46942040565380L) + ((int128)tmp_q[2] * 84279041682612L) + ((int128)tmp_q[3] * 26204854603393L) - ((int128)tmp_q[4] * 133544437612879L) + ((int128)tmp_q[5] * 106975586913695L) - ((int128)tmp_q[6] * 106006152409769L) + ((int128)tmp_q[7] * 98486624093119L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 98486624093119L) + ((int128)tmp_q[1] * 46377049270410L) + ((-((int128)tmp_q[2] * 46942040565380L) + ((int128)tmp_q[3] * 84279041682612L) + ((int128)tmp_q[4] * 26204854603393L) - ((int128)tmp_q[5] * 133544437612879L) + ((int128)tmp_q[6] * 106975586913695L) - ((int128)tmp_q[7] * 106006152409769L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 106006152409769L) + ((int128)tmp_q[1] * 98486624093119L) + ((int128)tmp_q[2] * 46377049270410L) + ((-((int128)tmp_q[3] * 46942040565380L) + ((int128)tmp_q[4] * 84279041682612L) + ((int128)tmp_q[5] * 26204854603393L) - ((int128)tmp_q[6] * 133544437612879L) + ((int128)tmp_q[7] * 106975586913695L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 106975586913695L) - ((int128)tmp_q[1] * 106006152409769L) + ((int128)tmp_q[2] * 98486624093119L) + ((int128)tmp_q[3] * 46377049270410L) + ((-((int128)tmp_q[4] * 46942040565380L) + ((int128)tmp_q[5] * 84279041682612L) + ((int128)tmp_q[6] * 26204854603393L) - ((int128)tmp_q[7] * 133544437612879L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 133544437612879L) + ((int128)tmp_q[1] * 106975586913695L) - ((int128)tmp_q[2] * 106006152409769L) + ((int128)tmp_q[3] * 98486624093119L) + ((int128)tmp_q[4] * 46377049270410L) + ((-((int128)tmp_q[5] * 46942040565380L) + ((int128)tmp_q[6] * 84279041682612L) + ((int128)tmp_q[7] * 26204854603393L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 26204854603393L) - ((int128)tmp_q[1] * 133544437612879L) + ((int128)tmp_q[2] * 106975586913695L) - ((int128)tmp_q[3] * 106006152409769L) + ((int128)tmp_q[4] * 98486624093119L) + ((int128)tmp_q[5] * 46377049270410L) + ((-((int128)tmp_q[6] * 46942040565380L) + ((int128)tmp_q[7] * 84279041682612L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 84279041682612L) + ((int128)tmp_q[1] * 26204854603393L) - ((int128)tmp_q[2] * 133544437612879L) + ((int128)tmp_q[3] * 106975586913695L) - ((int128)tmp_q[4] * 106006152409769L) + ((int128)tmp_q[5] * 98486624093119L) + ((int128)tmp_q[6] * 46377049270410L) - ((int128)tmp_q[7] * 328594283957660L);
	tmp_zero[7] = -((int128)tmp_q[0] * 46942040565380L) + ((int128)tmp_q[1] * 84279041682612L) + ((int128)tmp_q[2] * 26204854603393L) - ((int128)tmp_q[3] * 133544437612879L) + ((int128)tmp_q[4] * 106975586913695L) - ((int128)tmp_q[5] * 106006152409769L) + ((int128)tmp_q[6] * 98486624093119L) + ((int128)tmp_q[7] * 46377049270410L);

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

