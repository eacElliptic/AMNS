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
	tmp_q[0] = ((uint64_t)op[0] * 3664029019510608413UL) + ((((uint64_t)op[1] * 6956578364029959957UL) + ((uint64_t)op[2] * 13550543776688132187UL) + ((uint64_t)op[3] * 8223372864140780527UL) + ((uint64_t)op[4] * 14951958020510598410UL) + ((uint64_t)op[5] * 10687214165320278058UL) + ((uint64_t)op[6] * 7951085809770910855UL) + ((uint64_t)op[7] * 9494338042660551362UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 9494338042660551362UL) + ((uint64_t)op[1] * 3664029019510608413UL) + ((((uint64_t)op[2] * 6956578364029959957UL) + ((uint64_t)op[3] * 13550543776688132187UL) + ((uint64_t)op[4] * 8223372864140780527UL) + ((uint64_t)op[5] * 14951958020510598410UL) + ((uint64_t)op[6] * 10687214165320278058UL) + ((uint64_t)op[7] * 7951085809770910855UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 7951085809770910855UL) + ((uint64_t)op[1] * 9494338042660551362UL) + ((uint64_t)op[2] * 3664029019510608413UL) + ((((uint64_t)op[3] * 6956578364029959957UL) + ((uint64_t)op[4] * 13550543776688132187UL) + ((uint64_t)op[5] * 8223372864140780527UL) + ((uint64_t)op[6] * 14951958020510598410UL) + ((uint64_t)op[7] * 10687214165320278058UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 10687214165320278058UL) + ((uint64_t)op[1] * 7951085809770910855UL) + ((uint64_t)op[2] * 9494338042660551362UL) + ((uint64_t)op[3] * 3664029019510608413UL) + ((((uint64_t)op[4] * 6956578364029959957UL) + ((uint64_t)op[5] * 13550543776688132187UL) + ((uint64_t)op[6] * 8223372864140780527UL) + ((uint64_t)op[7] * 14951958020510598410UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 14951958020510598410UL) + ((uint64_t)op[1] * 10687214165320278058UL) + ((uint64_t)op[2] * 7951085809770910855UL) + ((uint64_t)op[3] * 9494338042660551362UL) + ((uint64_t)op[4] * 3664029019510608413UL) + ((((uint64_t)op[5] * 6956578364029959957UL) + ((uint64_t)op[6] * 13550543776688132187UL) + ((uint64_t)op[7] * 8223372864140780527UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 8223372864140780527UL) + ((uint64_t)op[1] * 14951958020510598410UL) + ((uint64_t)op[2] * 10687214165320278058UL) + ((uint64_t)op[3] * 7951085809770910855UL) + ((uint64_t)op[4] * 9494338042660551362UL) + ((uint64_t)op[5] * 3664029019510608413UL) + ((((uint64_t)op[6] * 6956578364029959957UL) + ((uint64_t)op[7] * 13550543776688132187UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 13550543776688132187UL) + ((uint64_t)op[1] * 8223372864140780527UL) + ((uint64_t)op[2] * 14951958020510598410UL) + ((uint64_t)op[3] * 10687214165320278058UL) + ((uint64_t)op[4] * 7951085809770910855UL) + ((uint64_t)op[5] * 9494338042660551362UL) + ((uint64_t)op[6] * 3664029019510608413UL) + ((uint64_t)op[7] * 11802560400790616467UL);
	tmp_q[7] = ((uint64_t)op[0] * 6956578364029959957UL) + ((uint64_t)op[1] * 13550543776688132187UL) + ((uint64_t)op[2] * 8223372864140780527UL) + ((uint64_t)op[3] * 14951958020510598410UL) + ((uint64_t)op[4] * 10687214165320278058UL) + ((uint64_t)op[5] * 7951085809770910855UL) + ((uint64_t)op[6] * 9494338042660551362UL) + ((uint64_t)op[7] * 3664029019510608413UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 95322352035825L) + ((((int128)tmp_q[1] * 116618292816492L) - ((int128)tmp_q[2] * 81418596777982L) - ((int128)tmp_q[3] * 51028938291694L) + ((int128)tmp_q[4] * 51385566222504L) - ((int128)tmp_q[5] * 22863315104899L) - ((int128)tmp_q[6] * 127328677645674L) - ((int128)tmp_q[7] * 18725894713773L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 18725894713773L) + ((int128)tmp_q[1] * 95322352035825L) + ((((int128)tmp_q[2] * 116618292816492L) - ((int128)tmp_q[3] * 81418596777982L) - ((int128)tmp_q[4] * 51028938291694L) + ((int128)tmp_q[5] * 51385566222504L) - ((int128)tmp_q[6] * 22863315104899L) - ((int128)tmp_q[7] * 127328677645674L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 127328677645674L) - ((int128)tmp_q[1] * 18725894713773L) + ((int128)tmp_q[2] * 95322352035825L) + ((((int128)tmp_q[3] * 116618292816492L) - ((int128)tmp_q[4] * 81418596777982L) - ((int128)tmp_q[5] * 51028938291694L) + ((int128)tmp_q[6] * 51385566222504L) - ((int128)tmp_q[7] * 22863315104899L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 22863315104899L) - ((int128)tmp_q[1] * 127328677645674L) - ((int128)tmp_q[2] * 18725894713773L) + ((int128)tmp_q[3] * 95322352035825L) + ((((int128)tmp_q[4] * 116618292816492L) - ((int128)tmp_q[5] * 81418596777982L) - ((int128)tmp_q[6] * 51028938291694L) + ((int128)tmp_q[7] * 51385566222504L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 51385566222504L) - ((int128)tmp_q[1] * 22863315104899L) - ((int128)tmp_q[2] * 127328677645674L) - ((int128)tmp_q[3] * 18725894713773L) + ((int128)tmp_q[4] * 95322352035825L) + ((((int128)tmp_q[5] * 116618292816492L) - ((int128)tmp_q[6] * 81418596777982L) - ((int128)tmp_q[7] * 51028938291694L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 51028938291694L) + ((int128)tmp_q[1] * 51385566222504L) - ((int128)tmp_q[2] * 22863315104899L) - ((int128)tmp_q[3] * 127328677645674L) - ((int128)tmp_q[4] * 18725894713773L) + ((int128)tmp_q[5] * 95322352035825L) + ((((int128)tmp_q[6] * 116618292816492L) - ((int128)tmp_q[7] * 81418596777982L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 81418596777982L) - ((int128)tmp_q[1] * 51028938291694L) + ((int128)tmp_q[2] * 51385566222504L) - ((int128)tmp_q[3] * 22863315104899L) - ((int128)tmp_q[4] * 127328677645674L) - ((int128)tmp_q[5] * 18725894713773L) + ((int128)tmp_q[6] * 95322352035825L) + ((int128)tmp_q[7] * 816328049715444L);
	tmp_zero[7] = ((int128)tmp_q[0] * 116618292816492L) - ((int128)tmp_q[1] * 81418596777982L) - ((int128)tmp_q[2] * 51028938291694L) + ((int128)tmp_q[3] * 51385566222504L) - ((int128)tmp_q[4] * 22863315104899L) - ((int128)tmp_q[5] * 127328677645674L) - ((int128)tmp_q[6] * 18725894713773L) + ((int128)tmp_q[7] * 95322352035825L);

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

