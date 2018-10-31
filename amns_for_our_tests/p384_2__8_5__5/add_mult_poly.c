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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12270211996601333961UL) + ((((uint64_t)op[1] * 11381364077378960026UL) + ((uint64_t)op[2] * 17144223289320924379UL) + ((uint64_t)op[3] * 92243857812892498UL) + ((uint64_t)op[4] * 9595673782333880592UL) + ((uint64_t)op[5] * 15225204007191025947UL) + ((uint64_t)op[6] * 60862535160684387UL) + ((uint64_t)op[7] * 12418497341878075359UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 12418497341878075359UL) + ((uint64_t)op[1] * 12270211996601333961UL) + ((((uint64_t)op[2] * 11381364077378960026UL) + ((uint64_t)op[3] * 17144223289320924379UL) + ((uint64_t)op[4] * 92243857812892498UL) + ((uint64_t)op[5] * 9595673782333880592UL) + ((uint64_t)op[6] * 15225204007191025947UL) + ((uint64_t)op[7] * 60862535160684387UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 60862535160684387UL) + ((uint64_t)op[1] * 12418497341878075359UL) + ((uint64_t)op[2] * 12270211996601333961UL) + ((((uint64_t)op[3] * 11381364077378960026UL) + ((uint64_t)op[4] * 17144223289320924379UL) + ((uint64_t)op[5] * 92243857812892498UL) + ((uint64_t)op[6] * 9595673782333880592UL) + ((uint64_t)op[7] * 15225204007191025947UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 15225204007191025947UL) + ((uint64_t)op[1] * 60862535160684387UL) + ((uint64_t)op[2] * 12418497341878075359UL) + ((uint64_t)op[3] * 12270211996601333961UL) + ((((uint64_t)op[4] * 11381364077378960026UL) + ((uint64_t)op[5] * 17144223289320924379UL) + ((uint64_t)op[6] * 92243857812892498UL) + ((uint64_t)op[7] * 9595673782333880592UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 9595673782333880592UL) + ((uint64_t)op[1] * 15225204007191025947UL) + ((uint64_t)op[2] * 60862535160684387UL) + ((uint64_t)op[3] * 12418497341878075359UL) + ((uint64_t)op[4] * 12270211996601333961UL) + ((((uint64_t)op[5] * 11381364077378960026UL) + ((uint64_t)op[6] * 17144223289320924379UL) + ((uint64_t)op[7] * 92243857812892498UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 92243857812892498UL) + ((uint64_t)op[1] * 9595673782333880592UL) + ((uint64_t)op[2] * 15225204007191025947UL) + ((uint64_t)op[3] * 60862535160684387UL) + ((uint64_t)op[4] * 12418497341878075359UL) + ((uint64_t)op[5] * 12270211996601333961UL) + ((((uint64_t)op[6] * 11381364077378960026UL) + ((uint64_t)op[7] * 17144223289320924379UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 17144223289320924379UL) + ((uint64_t)op[1] * 92243857812892498UL) + ((uint64_t)op[2] * 9595673782333880592UL) + ((uint64_t)op[3] * 15225204007191025947UL) + ((uint64_t)op[4] * 60862535160684387UL) + ((uint64_t)op[5] * 12418497341878075359UL) + ((uint64_t)op[6] * 12270211996601333961UL) + ((uint64_t)op[7] * 1566588165766145282UL);
	tmp_q[7] = ((uint64_t)op[0] * 11381364077378960026UL) + ((uint64_t)op[1] * 17144223289320924379UL) + ((uint64_t)op[2] * 92243857812892498UL) + ((uint64_t)op[3] * 9595673782333880592UL) + ((uint64_t)op[4] * 15225204007191025947UL) + ((uint64_t)op[5] * 60862535160684387UL) + ((uint64_t)op[6] * 12418497341878075359UL) + ((uint64_t)op[7] * 12270211996601333961UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 126109919135213L) + ((-((int128)tmp_q[1] * 79101463793897L) - ((int128)tmp_q[2] * 27228293940888L) + ((int128)tmp_q[3] * 125642135104891L) - ((int128)tmp_q[4] * 15130689064028L) + ((int128)tmp_q[5] * 88055635701946L) + ((int128)tmp_q[6] * 10449742006392L) - ((int128)tmp_q[7] * 121706339176274L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 121706339176274L) + ((int128)tmp_q[1] * 126109919135213L) + ((-((int128)tmp_q[2] * 79101463793897L) - ((int128)tmp_q[3] * 27228293940888L) + ((int128)tmp_q[4] * 125642135104891L) - ((int128)tmp_q[5] * 15130689064028L) + ((int128)tmp_q[6] * 88055635701946L) + ((int128)tmp_q[7] * 10449742006392L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 10449742006392L) - ((int128)tmp_q[1] * 121706339176274L) + ((int128)tmp_q[2] * 126109919135213L) + ((-((int128)tmp_q[3] * 79101463793897L) - ((int128)tmp_q[4] * 27228293940888L) + ((int128)tmp_q[5] * 125642135104891L) - ((int128)tmp_q[6] * 15130689064028L) + ((int128)tmp_q[7] * 88055635701946L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 88055635701946L) + ((int128)tmp_q[1] * 10449742006392L) - ((int128)tmp_q[2] * 121706339176274L) + ((int128)tmp_q[3] * 126109919135213L) + ((-((int128)tmp_q[4] * 79101463793897L) - ((int128)tmp_q[5] * 27228293940888L) + ((int128)tmp_q[6] * 125642135104891L) - ((int128)tmp_q[7] * 15130689064028L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 15130689064028L) + ((int128)tmp_q[1] * 88055635701946L) + ((int128)tmp_q[2] * 10449742006392L) - ((int128)tmp_q[3] * 121706339176274L) + ((int128)tmp_q[4] * 126109919135213L) + ((-((int128)tmp_q[5] * 79101463793897L) - ((int128)tmp_q[6] * 27228293940888L) + ((int128)tmp_q[7] * 125642135104891L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 125642135104891L) - ((int128)tmp_q[1] * 15130689064028L) + ((int128)tmp_q[2] * 88055635701946L) + ((int128)tmp_q[3] * 10449742006392L) - ((int128)tmp_q[4] * 121706339176274L) + ((int128)tmp_q[5] * 126109919135213L) + ((-((int128)tmp_q[6] * 79101463793897L) - ((int128)tmp_q[7] * 27228293940888L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 27228293940888L) + ((int128)tmp_q[1] * 125642135104891L) - ((int128)tmp_q[2] * 15130689064028L) + ((int128)tmp_q[3] * 88055635701946L) + ((int128)tmp_q[4] * 10449742006392L) - ((int128)tmp_q[5] * 121706339176274L) + ((int128)tmp_q[6] * 126109919135213L) - ((int128)tmp_q[7] * 395507318969485L);
	tmp_zero[7] = -((int128)tmp_q[0] * 79101463793897L) - ((int128)tmp_q[1] * 27228293940888L) + ((int128)tmp_q[2] * 125642135104891L) - ((int128)tmp_q[3] * 15130689064028L) + ((int128)tmp_q[4] * 88055635701946L) + ((int128)tmp_q[5] * 10449742006392L) - ((int128)tmp_q[6] * 121706339176274L) + ((int128)tmp_q[7] * 126109919135213L);

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

