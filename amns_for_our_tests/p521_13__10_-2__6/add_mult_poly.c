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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[9] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14404640412319668243UL) + ((((uint64_t)op[1] * 4999845156299562036UL) + ((uint64_t)op[2] * 5607572832843259142UL) + ((uint64_t)op[3] * 17231123867865659150UL) + ((uint64_t)op[4] * 1842361179046800940UL) + ((uint64_t)op[5] * 7904643730151119742UL) + ((uint64_t)op[6] * 9181300598519214998UL) + ((uint64_t)op[7] * 10870649077256938336UL) + ((uint64_t)op[8] * 8573294372331940552UL) + ((uint64_t)op[9] * 4819168460625324940UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 4819168460625324940UL) + ((uint64_t)op[1] * 14404640412319668243UL) + ((((uint64_t)op[2] * 4999845156299562036UL) + ((uint64_t)op[3] * 5607572832843259142UL) + ((uint64_t)op[4] * 17231123867865659150UL) + ((uint64_t)op[5] * 1842361179046800940UL) + ((uint64_t)op[6] * 7904643730151119742UL) + ((uint64_t)op[7] * 9181300598519214998UL) + ((uint64_t)op[8] * 10870649077256938336UL) + ((uint64_t)op[9] * 8573294372331940552UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 8573294372331940552UL) + ((uint64_t)op[1] * 4819168460625324940UL) + ((uint64_t)op[2] * 14404640412319668243UL) + ((((uint64_t)op[3] * 4999845156299562036UL) + ((uint64_t)op[4] * 5607572832843259142UL) + ((uint64_t)op[5] * 17231123867865659150UL) + ((uint64_t)op[6] * 1842361179046800940UL) + ((uint64_t)op[7] * 7904643730151119742UL) + ((uint64_t)op[8] * 9181300598519214998UL) + ((uint64_t)op[9] * 10870649077256938336UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 10870649077256938336UL) + ((uint64_t)op[1] * 8573294372331940552UL) + ((uint64_t)op[2] * 4819168460625324940UL) + ((uint64_t)op[3] * 14404640412319668243UL) + ((((uint64_t)op[4] * 4999845156299562036UL) + ((uint64_t)op[5] * 5607572832843259142UL) + ((uint64_t)op[6] * 17231123867865659150UL) + ((uint64_t)op[7] * 1842361179046800940UL) + ((uint64_t)op[8] * 7904643730151119742UL) + ((uint64_t)op[9] * 9181300598519214998UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 9181300598519214998UL) + ((uint64_t)op[1] * 10870649077256938336UL) + ((uint64_t)op[2] * 8573294372331940552UL) + ((uint64_t)op[3] * 4819168460625324940UL) + ((uint64_t)op[4] * 14404640412319668243UL) + ((((uint64_t)op[5] * 4999845156299562036UL) + ((uint64_t)op[6] * 5607572832843259142UL) + ((uint64_t)op[7] * 17231123867865659150UL) + ((uint64_t)op[8] * 1842361179046800940UL) + ((uint64_t)op[9] * 7904643730151119742UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 7904643730151119742UL) + ((uint64_t)op[1] * 9181300598519214998UL) + ((uint64_t)op[2] * 10870649077256938336UL) + ((uint64_t)op[3] * 8573294372331940552UL) + ((uint64_t)op[4] * 4819168460625324940UL) + ((uint64_t)op[5] * 14404640412319668243UL) + ((((uint64_t)op[6] * 4999845156299562036UL) + ((uint64_t)op[7] * 5607572832843259142UL) + ((uint64_t)op[8] * 17231123867865659150UL) + ((uint64_t)op[9] * 1842361179046800940UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 1842361179046800940UL) + ((uint64_t)op[1] * 7904643730151119742UL) + ((uint64_t)op[2] * 9181300598519214998UL) + ((uint64_t)op[3] * 10870649077256938336UL) + ((uint64_t)op[4] * 8573294372331940552UL) + ((uint64_t)op[5] * 4819168460625324940UL) + ((uint64_t)op[6] * 14404640412319668243UL) + ((((uint64_t)op[7] * 4999845156299562036UL) + ((uint64_t)op[8] * 5607572832843259142UL) + ((uint64_t)op[9] * 17231123867865659150UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 17231123867865659150UL) + ((uint64_t)op[1] * 1842361179046800940UL) + ((uint64_t)op[2] * 7904643730151119742UL) + ((uint64_t)op[3] * 9181300598519214998UL) + ((uint64_t)op[4] * 10870649077256938336UL) + ((uint64_t)op[5] * 8573294372331940552UL) + ((uint64_t)op[6] * 4819168460625324940UL) + ((uint64_t)op[7] * 14404640412319668243UL) + ((((uint64_t)op[8] * 4999845156299562036UL) + ((uint64_t)op[9] * 5607572832843259142UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 5607572832843259142UL) + ((uint64_t)op[1] * 17231123867865659150UL) + ((uint64_t)op[2] * 1842361179046800940UL) + ((uint64_t)op[3] * 7904643730151119742UL) + ((uint64_t)op[4] * 9181300598519214998UL) + ((uint64_t)op[5] * 10870649077256938336UL) + ((uint64_t)op[6] * 8573294372331940552UL) + ((uint64_t)op[7] * 4819168460625324940UL) + ((uint64_t)op[8] * 14404640412319668243UL) + ((uint64_t)op[9] * 8447053761110427544UL);
	tmp_q[9] = ((uint64_t)op[0] * 4999845156299562036UL) + ((uint64_t)op[1] * 5607572832843259142UL) + ((uint64_t)op[2] * 17231123867865659150UL) + ((uint64_t)op[3] * 1842361179046800940UL) + ((uint64_t)op[4] * 7904643730151119742UL) + ((uint64_t)op[5] * 9181300598519214998UL) + ((uint64_t)op[6] * 10870649077256938336UL) + ((uint64_t)op[7] * 8573294372331940552UL) + ((uint64_t)op[8] * 4819168460625324940UL) + ((uint64_t)op[9] * 14404640412319668243UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 897494113030397L) - ((((int128)tmp_q[1] * 1078823401354540L) - ((int128)tmp_q[2] * 828300062112966L) - ((int128)tmp_q[3] * 1061275849892818L) - ((int128)tmp_q[4] * 1713147644876124L) - ((int128)tmp_q[5] * 2571098790027282L) - ((int128)tmp_q[6] * 1844957467998610L) + ((int128)tmp_q[7] * 2699132360489056L) - ((int128)tmp_q[8] * 877975421296440L) - ((int128)tmp_q[9] * 1098509906391716L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1098509906391716L) + ((int128)tmp_q[1] * 897494113030397L) - ((((int128)tmp_q[2] * 1078823401354540L) - ((int128)tmp_q[3] * 828300062112966L) - ((int128)tmp_q[4] * 1061275849892818L) - ((int128)tmp_q[5] * 1713147644876124L) - ((int128)tmp_q[6] * 2571098790027282L) - ((int128)tmp_q[7] * 1844957467998610L) + ((int128)tmp_q[8] * 2699132360489056L) - ((int128)tmp_q[9] * 877975421296440L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 877975421296440L) - ((int128)tmp_q[1] * 1098509906391716L) + ((int128)tmp_q[2] * 897494113030397L) - ((((int128)tmp_q[3] * 1078823401354540L) - ((int128)tmp_q[4] * 828300062112966L) - ((int128)tmp_q[5] * 1061275849892818L) - ((int128)tmp_q[6] * 1713147644876124L) - ((int128)tmp_q[7] * 2571098790027282L) - ((int128)tmp_q[8] * 1844957467998610L) + ((int128)tmp_q[9] * 2699132360489056L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 2699132360489056L) - ((int128)tmp_q[1] * 877975421296440L) - ((int128)tmp_q[2] * 1098509906391716L) + ((int128)tmp_q[3] * 897494113030397L) - ((((int128)tmp_q[4] * 1078823401354540L) - ((int128)tmp_q[5] * 828300062112966L) - ((int128)tmp_q[6] * 1061275849892818L) - ((int128)tmp_q[7] * 1713147644876124L) - ((int128)tmp_q[8] * 2571098790027282L) - ((int128)tmp_q[9] * 1844957467998610L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1844957467998610L) + ((int128)tmp_q[1] * 2699132360489056L) - ((int128)tmp_q[2] * 877975421296440L) - ((int128)tmp_q[3] * 1098509906391716L) + ((int128)tmp_q[4] * 897494113030397L) - ((((int128)tmp_q[5] * 1078823401354540L) - ((int128)tmp_q[6] * 828300062112966L) - ((int128)tmp_q[7] * 1061275849892818L) - ((int128)tmp_q[8] * 1713147644876124L) - ((int128)tmp_q[9] * 2571098790027282L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 2571098790027282L) - ((int128)tmp_q[1] * 1844957467998610L) + ((int128)tmp_q[2] * 2699132360489056L) - ((int128)tmp_q[3] * 877975421296440L) - ((int128)tmp_q[4] * 1098509906391716L) + ((int128)tmp_q[5] * 897494113030397L) - ((((int128)tmp_q[6] * 1078823401354540L) - ((int128)tmp_q[7] * 828300062112966L) - ((int128)tmp_q[8] * 1061275849892818L) - ((int128)tmp_q[9] * 1713147644876124L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 1713147644876124L) - ((int128)tmp_q[1] * 2571098790027282L) - ((int128)tmp_q[2] * 1844957467998610L) + ((int128)tmp_q[3] * 2699132360489056L) - ((int128)tmp_q[4] * 877975421296440L) - ((int128)tmp_q[5] * 1098509906391716L) + ((int128)tmp_q[6] * 897494113030397L) - ((((int128)tmp_q[7] * 1078823401354540L) - ((int128)tmp_q[8] * 828300062112966L) - ((int128)tmp_q[9] * 1061275849892818L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 1061275849892818L) - ((int128)tmp_q[1] * 1713147644876124L) - ((int128)tmp_q[2] * 2571098790027282L) - ((int128)tmp_q[3] * 1844957467998610L) + ((int128)tmp_q[4] * 2699132360489056L) - ((int128)tmp_q[5] * 877975421296440L) - ((int128)tmp_q[6] * 1098509906391716L) + ((int128)tmp_q[7] * 897494113030397L) - ((((int128)tmp_q[8] * 1078823401354540L) - ((int128)tmp_q[9] * 828300062112966L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 828300062112966L) - ((int128)tmp_q[1] * 1061275849892818L) - ((int128)tmp_q[2] * 1713147644876124L) - ((int128)tmp_q[3] * 2571098790027282L) - ((int128)tmp_q[4] * 1844957467998610L) + ((int128)tmp_q[5] * 2699132360489056L) - ((int128)tmp_q[6] * 877975421296440L) - ((int128)tmp_q[7] * 1098509906391716L) + ((int128)tmp_q[8] * 897494113030397L) - ((int128)tmp_q[9] * 2157646802709080L);
	tmp_zero[9] = ((int128)tmp_q[0] * 1078823401354540L) - ((int128)tmp_q[1] * 828300062112966L) - ((int128)tmp_q[2] * 1061275849892818L) - ((int128)tmp_q[3] * 1713147644876124L) - ((int128)tmp_q[4] * 2571098790027282L) - ((int128)tmp_q[5] * 1844957467998610L) + ((int128)tmp_q[6] * 2699132360489056L) - ((int128)tmp_q[7] * 877975421296440L) - ((int128)tmp_q[8] * 1098509906391716L) + ((int128)tmp_q[9] * 897494113030397L);

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
}

