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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12217779775044312857UL) + ((((uint64_t)op[1] * 1575406346337557549UL) + ((uint64_t)op[2] * 6506724869304209079UL) + ((uint64_t)op[3] * 6548769526438843968UL) + ((uint64_t)op[4] * 6151916648169663154UL) + ((uint64_t)op[5] * 16327681304234039592UL) + ((uint64_t)op[6] * 2082984686238398486UL) + ((uint64_t)op[7] * 18212950745406372340UL) + ((uint64_t)op[8] * 15082431264408174528UL) + ((uint64_t)op[9] * 6784505950645607240UL) + ((uint64_t)op[10] * 12158238110728817904UL) + ((uint64_t)op[11] * 14210398611436207207UL)) * 18446744073709551602);
	tmp_q[1] = ((uint64_t)op[0] * 14210398611436207207UL) + ((uint64_t)op[1] * 12217779775044312857UL) + ((((uint64_t)op[2] * 1575406346337557549UL) + ((uint64_t)op[3] * 6506724869304209079UL) + ((uint64_t)op[4] * 6548769526438843968UL) + ((uint64_t)op[5] * 6151916648169663154UL) + ((uint64_t)op[6] * 16327681304234039592UL) + ((uint64_t)op[7] * 2082984686238398486UL) + ((uint64_t)op[8] * 18212950745406372340UL) + ((uint64_t)op[9] * 15082431264408174528UL) + ((uint64_t)op[10] * 6784505950645607240UL) + ((uint64_t)op[11] * 12158238110728817904UL)) * 18446744073709551602);
	tmp_q[2] = ((uint64_t)op[0] * 12158238110728817904UL) + ((uint64_t)op[1] * 14210398611436207207UL) + ((uint64_t)op[2] * 12217779775044312857UL) + ((((uint64_t)op[3] * 1575406346337557549UL) + ((uint64_t)op[4] * 6506724869304209079UL) + ((uint64_t)op[5] * 6548769526438843968UL) + ((uint64_t)op[6] * 6151916648169663154UL) + ((uint64_t)op[7] * 16327681304234039592UL) + ((uint64_t)op[8] * 2082984686238398486UL) + ((uint64_t)op[9] * 18212950745406372340UL) + ((uint64_t)op[10] * 15082431264408174528UL) + ((uint64_t)op[11] * 6784505950645607240UL)) * 18446744073709551602);
	tmp_q[3] = ((uint64_t)op[0] * 6784505950645607240UL) + ((uint64_t)op[1] * 12158238110728817904UL) + ((uint64_t)op[2] * 14210398611436207207UL) + ((uint64_t)op[3] * 12217779775044312857UL) + ((((uint64_t)op[4] * 1575406346337557549UL) + ((uint64_t)op[5] * 6506724869304209079UL) + ((uint64_t)op[6] * 6548769526438843968UL) + ((uint64_t)op[7] * 6151916648169663154UL) + ((uint64_t)op[8] * 16327681304234039592UL) + ((uint64_t)op[9] * 2082984686238398486UL) + ((uint64_t)op[10] * 18212950745406372340UL) + ((uint64_t)op[11] * 15082431264408174528UL)) * 18446744073709551602);
	tmp_q[4] = ((uint64_t)op[0] * 15082431264408174528UL) + ((uint64_t)op[1] * 6784505950645607240UL) + ((uint64_t)op[2] * 12158238110728817904UL) + ((uint64_t)op[3] * 14210398611436207207UL) + ((uint64_t)op[4] * 12217779775044312857UL) + ((((uint64_t)op[5] * 1575406346337557549UL) + ((uint64_t)op[6] * 6506724869304209079UL) + ((uint64_t)op[7] * 6548769526438843968UL) + ((uint64_t)op[8] * 6151916648169663154UL) + ((uint64_t)op[9] * 16327681304234039592UL) + ((uint64_t)op[10] * 2082984686238398486UL) + ((uint64_t)op[11] * 18212950745406372340UL)) * 18446744073709551602);
	tmp_q[5] = ((uint64_t)op[0] * 18212950745406372340UL) + ((uint64_t)op[1] * 15082431264408174528UL) + ((uint64_t)op[2] * 6784505950645607240UL) + ((uint64_t)op[3] * 12158238110728817904UL) + ((uint64_t)op[4] * 14210398611436207207UL) + ((uint64_t)op[5] * 12217779775044312857UL) + ((((uint64_t)op[6] * 1575406346337557549UL) + ((uint64_t)op[7] * 6506724869304209079UL) + ((uint64_t)op[8] * 6548769526438843968UL) + ((uint64_t)op[9] * 6151916648169663154UL) + ((uint64_t)op[10] * 16327681304234039592UL) + ((uint64_t)op[11] * 2082984686238398486UL)) * 18446744073709551602);
	tmp_q[6] = ((uint64_t)op[0] * 2082984686238398486UL) + ((uint64_t)op[1] * 18212950745406372340UL) + ((uint64_t)op[2] * 15082431264408174528UL) + ((uint64_t)op[3] * 6784505950645607240UL) + ((uint64_t)op[4] * 12158238110728817904UL) + ((uint64_t)op[5] * 14210398611436207207UL) + ((uint64_t)op[6] * 12217779775044312857UL) + ((((uint64_t)op[7] * 1575406346337557549UL) + ((uint64_t)op[8] * 6506724869304209079UL) + ((uint64_t)op[9] * 6548769526438843968UL) + ((uint64_t)op[10] * 6151916648169663154UL) + ((uint64_t)op[11] * 16327681304234039592UL)) * 18446744073709551602);
	tmp_q[7] = ((uint64_t)op[0] * 16327681304234039592UL) + ((uint64_t)op[1] * 2082984686238398486UL) + ((uint64_t)op[2] * 18212950745406372340UL) + ((uint64_t)op[3] * 15082431264408174528UL) + ((uint64_t)op[4] * 6784505950645607240UL) + ((uint64_t)op[5] * 12158238110728817904UL) + ((uint64_t)op[6] * 14210398611436207207UL) + ((uint64_t)op[7] * 12217779775044312857UL) + ((((uint64_t)op[8] * 1575406346337557549UL) + ((uint64_t)op[9] * 6506724869304209079UL) + ((uint64_t)op[10] * 6548769526438843968UL) + ((uint64_t)op[11] * 6151916648169663154UL)) * 18446744073709551602);
	tmp_q[8] = ((uint64_t)op[0] * 6151916648169663154UL) + ((uint64_t)op[1] * 16327681304234039592UL) + ((uint64_t)op[2] * 2082984686238398486UL) + ((uint64_t)op[3] * 18212950745406372340UL) + ((uint64_t)op[4] * 15082431264408174528UL) + ((uint64_t)op[5] * 6784505950645607240UL) + ((uint64_t)op[6] * 12158238110728817904UL) + ((uint64_t)op[7] * 14210398611436207207UL) + ((uint64_t)op[8] * 12217779775044312857UL) + ((((uint64_t)op[9] * 1575406346337557549UL) + ((uint64_t)op[10] * 6506724869304209079UL) + ((uint64_t)op[11] * 6548769526438843968UL)) * 18446744073709551602);
	tmp_q[9] = ((uint64_t)op[0] * 6548769526438843968UL) + ((uint64_t)op[1] * 6151916648169663154UL) + ((uint64_t)op[2] * 16327681304234039592UL) + ((uint64_t)op[3] * 2082984686238398486UL) + ((uint64_t)op[4] * 18212950745406372340UL) + ((uint64_t)op[5] * 15082431264408174528UL) + ((uint64_t)op[6] * 6784505950645607240UL) + ((uint64_t)op[7] * 12158238110728817904UL) + ((uint64_t)op[8] * 14210398611436207207UL) + ((uint64_t)op[9] * 12217779775044312857UL) + ((((uint64_t)op[10] * 1575406346337557549UL) + ((uint64_t)op[11] * 6506724869304209079UL)) * 18446744073709551602);
	tmp_q[10] = ((uint64_t)op[0] * 6506724869304209079UL) + ((uint64_t)op[1] * 6548769526438843968UL) + ((uint64_t)op[2] * 6151916648169663154UL) + ((uint64_t)op[3] * 16327681304234039592UL) + ((uint64_t)op[4] * 2082984686238398486UL) + ((uint64_t)op[5] * 18212950745406372340UL) + ((uint64_t)op[6] * 15082431264408174528UL) + ((uint64_t)op[7] * 6784505950645607240UL) + ((uint64_t)op[8] * 12158238110728817904UL) + ((uint64_t)op[9] * 14210398611436207207UL) + ((uint64_t)op[10] * 12217779775044312857UL) + ((uint64_t)op[11] * 14837799298693297546UL);
	tmp_q[11] = ((uint64_t)op[0] * 1575406346337557549UL) + ((uint64_t)op[1] * 6506724869304209079UL) + ((uint64_t)op[2] * 6548769526438843968UL) + ((uint64_t)op[3] * 6151916648169663154UL) + ((uint64_t)op[4] * 16327681304234039592UL) + ((uint64_t)op[5] * 2082984686238398486UL) + ((uint64_t)op[6] * 18212950745406372340UL) + ((uint64_t)op[7] * 15082431264408174528UL) + ((uint64_t)op[8] * 6784505950645607240UL) + ((uint64_t)op[9] * 12158238110728817904UL) + ((uint64_t)op[10] * 14210398611436207207UL) + ((uint64_t)op[11] * 12217779775044312857UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3736758004087L) - ((((int128)tmp_q[1] * 4079010188924L) - ((int128)tmp_q[2] * 5582193606596L) + ((int128)tmp_q[3] * 3770777945937L) - ((int128)tmp_q[4] * 749391643099L) - ((int128)tmp_q[5] * 6319071913621L) + ((int128)tmp_q[6] * 3368594156833L) - ((int128)tmp_q[7] * 1419624041485L) + ((int128)tmp_q[8] * 3061086247823L) - ((int128)tmp_q[9] * 4029230822617L) + ((int128)tmp_q[10] * 5852120778951L) - ((int128)tmp_q[11] * 4374303583529L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 4374303583529L) + ((int128)tmp_q[1] * 3736758004087L) - ((((int128)tmp_q[2] * 4079010188924L) - ((int128)tmp_q[3] * 5582193606596L) + ((int128)tmp_q[4] * 3770777945937L) - ((int128)tmp_q[5] * 749391643099L) - ((int128)tmp_q[6] * 6319071913621L) + ((int128)tmp_q[7] * 3368594156833L) - ((int128)tmp_q[8] * 1419624041485L) + ((int128)tmp_q[9] * 3061086247823L) - ((int128)tmp_q[10] * 4029230822617L) + ((int128)tmp_q[11] * 5852120778951L)) * 14);
	tmp_zero[2] = ((int128)tmp_q[0] * 5852120778951L) - ((int128)tmp_q[1] * 4374303583529L) + ((int128)tmp_q[2] * 3736758004087L) - ((((int128)tmp_q[3] * 4079010188924L) - ((int128)tmp_q[4] * 5582193606596L) + ((int128)tmp_q[5] * 3770777945937L) - ((int128)tmp_q[6] * 749391643099L) - ((int128)tmp_q[7] * 6319071913621L) + ((int128)tmp_q[8] * 3368594156833L) - ((int128)tmp_q[9] * 1419624041485L) + ((int128)tmp_q[10] * 3061086247823L) - ((int128)tmp_q[11] * 4029230822617L)) * 14);
	tmp_zero[3] = -((int128)tmp_q[0] * 4029230822617L) + ((int128)tmp_q[1] * 5852120778951L) - ((int128)tmp_q[2] * 4374303583529L) + ((int128)tmp_q[3] * 3736758004087L) - ((((int128)tmp_q[4] * 4079010188924L) - ((int128)tmp_q[5] * 5582193606596L) + ((int128)tmp_q[6] * 3770777945937L) - ((int128)tmp_q[7] * 749391643099L) - ((int128)tmp_q[8] * 6319071913621L) + ((int128)tmp_q[9] * 3368594156833L) - ((int128)tmp_q[10] * 1419624041485L) + ((int128)tmp_q[11] * 3061086247823L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 3061086247823L) - ((int128)tmp_q[1] * 4029230822617L) + ((int128)tmp_q[2] * 5852120778951L) - ((int128)tmp_q[3] * 4374303583529L) + ((int128)tmp_q[4] * 3736758004087L) - ((((int128)tmp_q[5] * 4079010188924L) - ((int128)tmp_q[6] * 5582193606596L) + ((int128)tmp_q[7] * 3770777945937L) - ((int128)tmp_q[8] * 749391643099L) - ((int128)tmp_q[9] * 6319071913621L) + ((int128)tmp_q[10] * 3368594156833L) - ((int128)tmp_q[11] * 1419624041485L)) * 14);
	tmp_zero[5] = -((int128)tmp_q[0] * 1419624041485L) + ((int128)tmp_q[1] * 3061086247823L) - ((int128)tmp_q[2] * 4029230822617L) + ((int128)tmp_q[3] * 5852120778951L) - ((int128)tmp_q[4] * 4374303583529L) + ((int128)tmp_q[5] * 3736758004087L) - ((((int128)tmp_q[6] * 4079010188924L) - ((int128)tmp_q[7] * 5582193606596L) + ((int128)tmp_q[8] * 3770777945937L) - ((int128)tmp_q[9] * 749391643099L) - ((int128)tmp_q[10] * 6319071913621L) + ((int128)tmp_q[11] * 3368594156833L)) * 14);
	tmp_zero[6] = ((int128)tmp_q[0] * 3368594156833L) - ((int128)tmp_q[1] * 1419624041485L) + ((int128)tmp_q[2] * 3061086247823L) - ((int128)tmp_q[3] * 4029230822617L) + ((int128)tmp_q[4] * 5852120778951L) - ((int128)tmp_q[5] * 4374303583529L) + ((int128)tmp_q[6] * 3736758004087L) - ((((int128)tmp_q[7] * 4079010188924L) - ((int128)tmp_q[8] * 5582193606596L) + ((int128)tmp_q[9] * 3770777945937L) - ((int128)tmp_q[10] * 749391643099L) - ((int128)tmp_q[11] * 6319071913621L)) * 14);
	tmp_zero[7] = -((int128)tmp_q[0] * 6319071913621L) + ((int128)tmp_q[1] * 3368594156833L) - ((int128)tmp_q[2] * 1419624041485L) + ((int128)tmp_q[3] * 3061086247823L) - ((int128)tmp_q[4] * 4029230822617L) + ((int128)tmp_q[5] * 5852120778951L) - ((int128)tmp_q[6] * 4374303583529L) + ((int128)tmp_q[7] * 3736758004087L) - ((((int128)tmp_q[8] * 4079010188924L) - ((int128)tmp_q[9] * 5582193606596L) + ((int128)tmp_q[10] * 3770777945937L) - ((int128)tmp_q[11] * 749391643099L)) * 14);
	tmp_zero[8] = -((int128)tmp_q[0] * 749391643099L) - ((int128)tmp_q[1] * 6319071913621L) + ((int128)tmp_q[2] * 3368594156833L) - ((int128)tmp_q[3] * 1419624041485L) + ((int128)tmp_q[4] * 3061086247823L) - ((int128)tmp_q[5] * 4029230822617L) + ((int128)tmp_q[6] * 5852120778951L) - ((int128)tmp_q[7] * 4374303583529L) + ((int128)tmp_q[8] * 3736758004087L) - ((((int128)tmp_q[9] * 4079010188924L) - ((int128)tmp_q[10] * 5582193606596L) + ((int128)tmp_q[11] * 3770777945937L)) * 14);
	tmp_zero[9] = ((int128)tmp_q[0] * 3770777945937L) - ((int128)tmp_q[1] * 749391643099L) - ((int128)tmp_q[2] * 6319071913621L) + ((int128)tmp_q[3] * 3368594156833L) - ((int128)tmp_q[4] * 1419624041485L) + ((int128)tmp_q[5] * 3061086247823L) - ((int128)tmp_q[6] * 4029230822617L) + ((int128)tmp_q[7] * 5852120778951L) - ((int128)tmp_q[8] * 4374303583529L) + ((int128)tmp_q[9] * 3736758004087L) - ((((int128)tmp_q[10] * 4079010188924L) - ((int128)tmp_q[11] * 5582193606596L)) * 14);
	tmp_zero[10] = -((int128)tmp_q[0] * 5582193606596L) + ((int128)tmp_q[1] * 3770777945937L) - ((int128)tmp_q[2] * 749391643099L) - ((int128)tmp_q[3] * 6319071913621L) + ((int128)tmp_q[4] * 3368594156833L) - ((int128)tmp_q[5] * 1419624041485L) + ((int128)tmp_q[6] * 3061086247823L) - ((int128)tmp_q[7] * 4029230822617L) + ((int128)tmp_q[8] * 5852120778951L) - ((int128)tmp_q[9] * 4374303583529L) + ((int128)tmp_q[10] * 3736758004087L) - ((int128)tmp_q[11] * 57106142644936L);
	tmp_zero[11] = ((int128)tmp_q[0] * 4079010188924L) - ((int128)tmp_q[1] * 5582193606596L) + ((int128)tmp_q[2] * 3770777945937L) - ((int128)tmp_q[3] * 749391643099L) - ((int128)tmp_q[4] * 6319071913621L) + ((int128)tmp_q[5] * 3368594156833L) - ((int128)tmp_q[6] * 1419624041485L) + ((int128)tmp_q[7] * 3061086247823L) - ((int128)tmp_q[8] * 4029230822617L) + ((int128)tmp_q[9] * 5852120778951L) - ((int128)tmp_q[10] * 4374303583529L) + ((int128)tmp_q[11] * 3736758004087L);

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

