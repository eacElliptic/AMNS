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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14858014388295588991UL) + ((((uint64_t)op[1] * 10814680603542316822UL) + ((uint64_t)op[2] * 16912097618779553891UL) + ((uint64_t)op[3] * 6413548915157939964UL) + ((uint64_t)op[4] * 12368982343783832556UL) + ((uint64_t)op[5] * 17191368739010045695UL) + ((uint64_t)op[6] * 11534248871012361107UL) + ((uint64_t)op[7] * 9649315762040861895UL) + ((uint64_t)op[8] * 11573958958276316157UL) + ((uint64_t)op[9] * 1985380725696648162UL) + ((uint64_t)op[10] * 10698115972720211172UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 10698115972720211172UL) + ((uint64_t)op[1] * 14858014388295588991UL) + ((((uint64_t)op[2] * 10814680603542316822UL) + ((uint64_t)op[3] * 16912097618779553891UL) + ((uint64_t)op[4] * 6413548915157939964UL) + ((uint64_t)op[5] * 12368982343783832556UL) + ((uint64_t)op[6] * 17191368739010045695UL) + ((uint64_t)op[7] * 11534248871012361107UL) + ((uint64_t)op[8] * 9649315762040861895UL) + ((uint64_t)op[9] * 11573958958276316157UL) + ((uint64_t)op[10] * 1985380725696648162UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 1985380725696648162UL) + ((uint64_t)op[1] * 10698115972720211172UL) + ((uint64_t)op[2] * 14858014388295588991UL) + ((((uint64_t)op[3] * 10814680603542316822UL) + ((uint64_t)op[4] * 16912097618779553891UL) + ((uint64_t)op[5] * 6413548915157939964UL) + ((uint64_t)op[6] * 12368982343783832556UL) + ((uint64_t)op[7] * 17191368739010045695UL) + ((uint64_t)op[8] * 11534248871012361107UL) + ((uint64_t)op[9] * 9649315762040861895UL) + ((uint64_t)op[10] * 11573958958276316157UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 11573958958276316157UL) + ((uint64_t)op[1] * 1985380725696648162UL) + ((uint64_t)op[2] * 10698115972720211172UL) + ((uint64_t)op[3] * 14858014388295588991UL) + ((((uint64_t)op[4] * 10814680603542316822UL) + ((uint64_t)op[5] * 16912097618779553891UL) + ((uint64_t)op[6] * 6413548915157939964UL) + ((uint64_t)op[7] * 12368982343783832556UL) + ((uint64_t)op[8] * 17191368739010045695UL) + ((uint64_t)op[9] * 11534248871012361107UL) + ((uint64_t)op[10] * 9649315762040861895UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 9649315762040861895UL) + ((uint64_t)op[1] * 11573958958276316157UL) + ((uint64_t)op[2] * 1985380725696648162UL) + ((uint64_t)op[3] * 10698115972720211172UL) + ((uint64_t)op[4] * 14858014388295588991UL) + ((((uint64_t)op[5] * 10814680603542316822UL) + ((uint64_t)op[6] * 16912097618779553891UL) + ((uint64_t)op[7] * 6413548915157939964UL) + ((uint64_t)op[8] * 12368982343783832556UL) + ((uint64_t)op[9] * 17191368739010045695UL) + ((uint64_t)op[10] * 11534248871012361107UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 11534248871012361107UL) + ((uint64_t)op[1] * 9649315762040861895UL) + ((uint64_t)op[2] * 11573958958276316157UL) + ((uint64_t)op[3] * 1985380725696648162UL) + ((uint64_t)op[4] * 10698115972720211172UL) + ((uint64_t)op[5] * 14858014388295588991UL) + ((((uint64_t)op[6] * 10814680603542316822UL) + ((uint64_t)op[7] * 16912097618779553891UL) + ((uint64_t)op[8] * 6413548915157939964UL) + ((uint64_t)op[9] * 12368982343783832556UL) + ((uint64_t)op[10] * 17191368739010045695UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 17191368739010045695UL) + ((uint64_t)op[1] * 11534248871012361107UL) + ((uint64_t)op[2] * 9649315762040861895UL) + ((uint64_t)op[3] * 11573958958276316157UL) + ((uint64_t)op[4] * 1985380725696648162UL) + ((uint64_t)op[5] * 10698115972720211172UL) + ((uint64_t)op[6] * 14858014388295588991UL) + ((((uint64_t)op[7] * 10814680603542316822UL) + ((uint64_t)op[8] * 16912097618779553891UL) + ((uint64_t)op[9] * 6413548915157939964UL) + ((uint64_t)op[10] * 12368982343783832556UL)) * 18446744073709551608);
	tmp_q[7] = ((uint64_t)op[0] * 12368982343783832556UL) + ((uint64_t)op[1] * 17191368739010045695UL) + ((uint64_t)op[2] * 11534248871012361107UL) + ((uint64_t)op[3] * 9649315762040861895UL) + ((uint64_t)op[4] * 11573958958276316157UL) + ((uint64_t)op[5] * 1985380725696648162UL) + ((uint64_t)op[6] * 10698115972720211172UL) + ((uint64_t)op[7] * 14858014388295588991UL) + ((((uint64_t)op[8] * 10814680603542316822UL) + ((uint64_t)op[9] * 16912097618779553891UL) + ((uint64_t)op[10] * 6413548915157939964UL)) * 18446744073709551608);
	tmp_q[8] = ((uint64_t)op[0] * 6413548915157939964UL) + ((uint64_t)op[1] * 12368982343783832556UL) + ((uint64_t)op[2] * 17191368739010045695UL) + ((uint64_t)op[3] * 11534248871012361107UL) + ((uint64_t)op[4] * 9649315762040861895UL) + ((uint64_t)op[5] * 11573958958276316157UL) + ((uint64_t)op[6] * 1985380725696648162UL) + ((uint64_t)op[7] * 10698115972720211172UL) + ((uint64_t)op[8] * 14858014388295588991UL) + ((((uint64_t)op[9] * 10814680603542316822UL) + ((uint64_t)op[10] * 16912097618779553891UL)) * 18446744073709551608);
	tmp_q[9] = ((uint64_t)op[0] * 16912097618779553891UL) + ((uint64_t)op[1] * 6413548915157939964UL) + ((uint64_t)op[2] * 12368982343783832556UL) + ((uint64_t)op[3] * 17191368739010045695UL) + ((uint64_t)op[4] * 11534248871012361107UL) + ((uint64_t)op[5] * 9649315762040861895UL) + ((uint64_t)op[6] * 11573958958276316157UL) + ((uint64_t)op[7] * 1985380725696648162UL) + ((uint64_t)op[8] * 10698115972720211172UL) + ((uint64_t)op[9] * 14858014388295588991UL) + ((uint64_t)op[10] * 5716275540209223504UL);
	tmp_q[10] = ((uint64_t)op[0] * 10814680603542316822UL) + ((uint64_t)op[1] * 16912097618779553891UL) + ((uint64_t)op[2] * 6413548915157939964UL) + ((uint64_t)op[3] * 12368982343783832556UL) + ((uint64_t)op[4] * 17191368739010045695UL) + ((uint64_t)op[5] * 11534248871012361107UL) + ((uint64_t)op[6] * 9649315762040861895UL) + ((uint64_t)op[7] * 11573958958276316157UL) + ((uint64_t)op[8] * 1985380725696648162UL) + ((uint64_t)op[9] * 10698115972720211172UL) + ((uint64_t)op[10] * 14858014388295588991UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 91137891150095L) - ((((int128)tmp_q[1] * 33827399366844L) - ((int128)tmp_q[2] * 45196930029956L) - ((int128)tmp_q[3] * 51668878885527L) - ((int128)tmp_q[4] * 66474678313026L) + ((int128)tmp_q[5] * 50558539773468L) - ((int128)tmp_q[6] * 6120438578337L) + ((int128)tmp_q[7] * 88537843228083L) - ((int128)tmp_q[8] * 49494128098163L) + ((int128)tmp_q[9] * 33854729160722L) + ((int128)tmp_q[10] * 92654145922244L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 92654145922244L) - ((int128)tmp_q[1] * 91137891150095L) - ((((int128)tmp_q[2] * 33827399366844L) - ((int128)tmp_q[3] * 45196930029956L) - ((int128)tmp_q[4] * 51668878885527L) - ((int128)tmp_q[5] * 66474678313026L) + ((int128)tmp_q[6] * 50558539773468L) - ((int128)tmp_q[7] * 6120438578337L) + ((int128)tmp_q[8] * 88537843228083L) - ((int128)tmp_q[9] * 49494128098163L) + ((int128)tmp_q[10] * 33854729160722L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 33854729160722L) + ((int128)tmp_q[1] * 92654145922244L) - ((int128)tmp_q[2] * 91137891150095L) - ((((int128)tmp_q[3] * 33827399366844L) - ((int128)tmp_q[4] * 45196930029956L) - ((int128)tmp_q[5] * 51668878885527L) - ((int128)tmp_q[6] * 66474678313026L) + ((int128)tmp_q[7] * 50558539773468L) - ((int128)tmp_q[8] * 6120438578337L) + ((int128)tmp_q[9] * 88537843228083L) - ((int128)tmp_q[10] * 49494128098163L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 49494128098163L) + ((int128)tmp_q[1] * 33854729160722L) + ((int128)tmp_q[2] * 92654145922244L) - ((int128)tmp_q[3] * 91137891150095L) - ((((int128)tmp_q[4] * 33827399366844L) - ((int128)tmp_q[5] * 45196930029956L) - ((int128)tmp_q[6] * 51668878885527L) - ((int128)tmp_q[7] * 66474678313026L) + ((int128)tmp_q[8] * 50558539773468L) - ((int128)tmp_q[9] * 6120438578337L) + ((int128)tmp_q[10] * 88537843228083L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 88537843228083L) - ((int128)tmp_q[1] * 49494128098163L) + ((int128)tmp_q[2] * 33854729160722L) + ((int128)tmp_q[3] * 92654145922244L) - ((int128)tmp_q[4] * 91137891150095L) - ((((int128)tmp_q[5] * 33827399366844L) - ((int128)tmp_q[6] * 45196930029956L) - ((int128)tmp_q[7] * 51668878885527L) - ((int128)tmp_q[8] * 66474678313026L) + ((int128)tmp_q[9] * 50558539773468L) - ((int128)tmp_q[10] * 6120438578337L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 6120438578337L) + ((int128)tmp_q[1] * 88537843228083L) - ((int128)tmp_q[2] * 49494128098163L) + ((int128)tmp_q[3] * 33854729160722L) + ((int128)tmp_q[4] * 92654145922244L) - ((int128)tmp_q[5] * 91137891150095L) - ((((int128)tmp_q[6] * 33827399366844L) - ((int128)tmp_q[7] * 45196930029956L) - ((int128)tmp_q[8] * 51668878885527L) - ((int128)tmp_q[9] * 66474678313026L) + ((int128)tmp_q[10] * 50558539773468L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 50558539773468L) - ((int128)tmp_q[1] * 6120438578337L) + ((int128)tmp_q[2] * 88537843228083L) - ((int128)tmp_q[3] * 49494128098163L) + ((int128)tmp_q[4] * 33854729160722L) + ((int128)tmp_q[5] * 92654145922244L) - ((int128)tmp_q[6] * 91137891150095L) - ((((int128)tmp_q[7] * 33827399366844L) - ((int128)tmp_q[8] * 45196930029956L) - ((int128)tmp_q[9] * 51668878885527L) - ((int128)tmp_q[10] * 66474678313026L)) * 8);
	tmp_zero[7] = -((int128)tmp_q[0] * 66474678313026L) + ((int128)tmp_q[1] * 50558539773468L) - ((int128)tmp_q[2] * 6120438578337L) + ((int128)tmp_q[3] * 88537843228083L) - ((int128)tmp_q[4] * 49494128098163L) + ((int128)tmp_q[5] * 33854729160722L) + ((int128)tmp_q[6] * 92654145922244L) - ((int128)tmp_q[7] * 91137891150095L) - ((((int128)tmp_q[8] * 33827399366844L) - ((int128)tmp_q[9] * 45196930029956L) - ((int128)tmp_q[10] * 51668878885527L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 51668878885527L) - ((int128)tmp_q[1] * 66474678313026L) + ((int128)tmp_q[2] * 50558539773468L) - ((int128)tmp_q[3] * 6120438578337L) + ((int128)tmp_q[4] * 88537843228083L) - ((int128)tmp_q[5] * 49494128098163L) + ((int128)tmp_q[6] * 33854729160722L) + ((int128)tmp_q[7] * 92654145922244L) - ((int128)tmp_q[8] * 91137891150095L) - ((((int128)tmp_q[9] * 33827399366844L) - ((int128)tmp_q[10] * 45196930029956L)) * 8);
	tmp_zero[9] = -((int128)tmp_q[0] * 45196930029956L) - ((int128)tmp_q[1] * 51668878885527L) - ((int128)tmp_q[2] * 66474678313026L) + ((int128)tmp_q[3] * 50558539773468L) - ((int128)tmp_q[4] * 6120438578337L) + ((int128)tmp_q[5] * 88537843228083L) - ((int128)tmp_q[6] * 49494128098163L) + ((int128)tmp_q[7] * 33854729160722L) + ((int128)tmp_q[8] * 92654145922244L) - ((int128)tmp_q[9] * 91137891150095L) - ((int128)tmp_q[10] * 270619194934752L);
	tmp_zero[10] = ((int128)tmp_q[0] * 33827399366844L) - ((int128)tmp_q[1] * 45196930029956L) - ((int128)tmp_q[2] * 51668878885527L) - ((int128)tmp_q[3] * 66474678313026L) + ((int128)tmp_q[4] * 50558539773468L) - ((int128)tmp_q[5] * 6120438578337L) + ((int128)tmp_q[6] * 88537843228083L) - ((int128)tmp_q[7] * 49494128098163L) + ((int128)tmp_q[8] * 33854729160722L) + ((int128)tmp_q[9] * 92654145922244L) - ((int128)tmp_q[10] * 91137891150095L);

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

