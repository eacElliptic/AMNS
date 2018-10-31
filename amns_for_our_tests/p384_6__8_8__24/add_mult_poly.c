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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 4);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3852602748539471407UL) + ((((uint64_t)op[1] * 8659473103612201828UL) + ((uint64_t)op[2] * 17719515985682935607UL) + ((uint64_t)op[3] * 3325283878453716505UL) + ((uint64_t)op[4] * 8068708888058194828UL) + ((uint64_t)op[5] * 7764891021128700555UL) + ((uint64_t)op[6] * 4363299405391723147UL) + ((uint64_t)op[7] * 5476074584255113861UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 5476074584255113861UL) + ((uint64_t)op[1] * 3852602748539471407UL) + ((((uint64_t)op[2] * 8659473103612201828UL) + ((uint64_t)op[3] * 17719515985682935607UL) + ((uint64_t)op[4] * 3325283878453716505UL) + ((uint64_t)op[5] * 8068708888058194828UL) + ((uint64_t)op[6] * 7764891021128700555UL) + ((uint64_t)op[7] * 4363299405391723147UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 4363299405391723147UL) + ((uint64_t)op[1] * 5476074584255113861UL) + ((uint64_t)op[2] * 3852602748539471407UL) + ((((uint64_t)op[3] * 8659473103612201828UL) + ((uint64_t)op[4] * 17719515985682935607UL) + ((uint64_t)op[5] * 3325283878453716505UL) + ((uint64_t)op[6] * 8068708888058194828UL) + ((uint64_t)op[7] * 7764891021128700555UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 7764891021128700555UL) + ((uint64_t)op[1] * 4363299405391723147UL) + ((uint64_t)op[2] * 5476074584255113861UL) + ((uint64_t)op[3] * 3852602748539471407UL) + ((((uint64_t)op[4] * 8659473103612201828UL) + ((uint64_t)op[5] * 17719515985682935607UL) + ((uint64_t)op[6] * 3325283878453716505UL) + ((uint64_t)op[7] * 8068708888058194828UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 8068708888058194828UL) + ((uint64_t)op[1] * 7764891021128700555UL) + ((uint64_t)op[2] * 4363299405391723147UL) + ((uint64_t)op[3] * 5476074584255113861UL) + ((uint64_t)op[4] * 3852602748539471407UL) + ((((uint64_t)op[5] * 8659473103612201828UL) + ((uint64_t)op[6] * 17719515985682935607UL) + ((uint64_t)op[7] * 3325283878453716505UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 3325283878453716505UL) + ((uint64_t)op[1] * 8068708888058194828UL) + ((uint64_t)op[2] * 7764891021128700555UL) + ((uint64_t)op[3] * 4363299405391723147UL) + ((uint64_t)op[4] * 5476074584255113861UL) + ((uint64_t)op[5] * 3852602748539471407UL) + ((((uint64_t)op[6] * 8659473103612201828UL) + ((uint64_t)op[7] * 17719515985682935607UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 17719515985682935607UL) + ((uint64_t)op[1] * 3325283878453716505UL) + ((uint64_t)op[2] * 8068708888058194828UL) + ((uint64_t)op[3] * 7764891021128700555UL) + ((uint64_t)op[4] * 4363299405391723147UL) + ((uint64_t)op[5] * 5476074584255113861UL) + ((uint64_t)op[6] * 3852602748539471407UL) + ((uint64_t)op[7] * 13935552607768959776UL);
	tmp_q[7] = ((uint64_t)op[0] * 8659473103612201828UL) + ((uint64_t)op[1] * 17719515985682935607UL) + ((uint64_t)op[2] * 3325283878453716505UL) + ((uint64_t)op[3] * 8068708888058194828UL) + ((uint64_t)op[4] * 7764891021128700555UL) + ((uint64_t)op[5] * 4363299405391723147UL) + ((uint64_t)op[6] * 5476074584255113861UL) + ((uint64_t)op[7] * 3852602748539471407UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2789347986225L) + ((((int128)tmp_q[1] * 98555900239259L) - ((int128)tmp_q[2] * 50548984617823L) - ((int128)tmp_q[3] * 108494868110572L) + ((int128)tmp_q[4] * 74889330999077L) + ((int128)tmp_q[5] * 38538857120798L) - ((int128)tmp_q[6] * 48408179478804L) + ((int128)tmp_q[7] * 36750497881597L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 36750497881597L) + ((int128)tmp_q[1] * 2789347986225L) + ((((int128)tmp_q[2] * 98555900239259L) - ((int128)tmp_q[3] * 50548984617823L) - ((int128)tmp_q[4] * 108494868110572L) + ((int128)tmp_q[5] * 74889330999077L) + ((int128)tmp_q[6] * 38538857120798L) - ((int128)tmp_q[7] * 48408179478804L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 48408179478804L) + ((int128)tmp_q[1] * 36750497881597L) + ((int128)tmp_q[2] * 2789347986225L) + ((((int128)tmp_q[3] * 98555900239259L) - ((int128)tmp_q[4] * 50548984617823L) - ((int128)tmp_q[5] * 108494868110572L) + ((int128)tmp_q[6] * 74889330999077L) + ((int128)tmp_q[7] * 38538857120798L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 38538857120798L) - ((int128)tmp_q[1] * 48408179478804L) + ((int128)tmp_q[2] * 36750497881597L) + ((int128)tmp_q[3] * 2789347986225L) + ((((int128)tmp_q[4] * 98555900239259L) - ((int128)tmp_q[5] * 50548984617823L) - ((int128)tmp_q[6] * 108494868110572L) + ((int128)tmp_q[7] * 74889330999077L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 74889330999077L) + ((int128)tmp_q[1] * 38538857120798L) - ((int128)tmp_q[2] * 48408179478804L) + ((int128)tmp_q[3] * 36750497881597L) + ((int128)tmp_q[4] * 2789347986225L) + ((((int128)tmp_q[5] * 98555900239259L) - ((int128)tmp_q[6] * 50548984617823L) - ((int128)tmp_q[7] * 108494868110572L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 108494868110572L) + ((int128)tmp_q[1] * 74889330999077L) + ((int128)tmp_q[2] * 38538857120798L) - ((int128)tmp_q[3] * 48408179478804L) + ((int128)tmp_q[4] * 36750497881597L) + ((int128)tmp_q[5] * 2789347986225L) + ((((int128)tmp_q[6] * 98555900239259L) - ((int128)tmp_q[7] * 50548984617823L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 50548984617823L) - ((int128)tmp_q[1] * 108494868110572L) + ((int128)tmp_q[2] * 74889330999077L) + ((int128)tmp_q[3] * 38538857120798L) - ((int128)tmp_q[4] * 48408179478804L) + ((int128)tmp_q[5] * 36750497881597L) + ((int128)tmp_q[6] * 2789347986225L) + ((int128)tmp_q[7] * 788447201914072L);
	tmp_zero[7] = ((int128)tmp_q[0] * 98555900239259L) - ((int128)tmp_q[1] * 50548984617823L) - ((int128)tmp_q[2] * 108494868110572L) + ((int128)tmp_q[3] * 74889330999077L) + ((int128)tmp_q[4] * 38538857120798L) - ((int128)tmp_q[5] * 48408179478804L) + ((int128)tmp_q[6] * 36750497881597L) + ((int128)tmp_q[7] * 2789347986225L);

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

