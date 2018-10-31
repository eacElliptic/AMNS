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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3032100978319218926UL) + ((((uint64_t)op[1] * 2239388948852210279UL) + ((uint64_t)op[2] * 15273176679450849295UL) + ((uint64_t)op[3] * 12273388804262012310UL) + ((uint64_t)op[4] * 11221737545519111142UL) + ((uint64_t)op[5] * 12201587493475465996UL) + ((uint64_t)op[6] * 13678276380802697737UL) + ((uint64_t)op[7] * 7051930127756445995UL) + ((uint64_t)op[8] * 5742670642601390004UL) + ((uint64_t)op[9] * 6888864087897378859UL) + ((uint64_t)op[10] * 13162482694068408160UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 13162482694068408160UL) + ((uint64_t)op[1] * 3032100978319218926UL) + ((((uint64_t)op[2] * 2239388948852210279UL) + ((uint64_t)op[3] * 15273176679450849295UL) + ((uint64_t)op[4] * 12273388804262012310UL) + ((uint64_t)op[5] * 11221737545519111142UL) + ((uint64_t)op[6] * 12201587493475465996UL) + ((uint64_t)op[7] * 13678276380802697737UL) + ((uint64_t)op[8] * 7051930127756445995UL) + ((uint64_t)op[9] * 5742670642601390004UL) + ((uint64_t)op[10] * 6888864087897378859UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 6888864087897378859UL) + ((uint64_t)op[1] * 13162482694068408160UL) + ((uint64_t)op[2] * 3032100978319218926UL) + ((((uint64_t)op[3] * 2239388948852210279UL) + ((uint64_t)op[4] * 15273176679450849295UL) + ((uint64_t)op[5] * 12273388804262012310UL) + ((uint64_t)op[6] * 11221737545519111142UL) + ((uint64_t)op[7] * 12201587493475465996UL) + ((uint64_t)op[8] * 13678276380802697737UL) + ((uint64_t)op[9] * 7051930127756445995UL) + ((uint64_t)op[10] * 5742670642601390004UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 5742670642601390004UL) + ((uint64_t)op[1] * 6888864087897378859UL) + ((uint64_t)op[2] * 13162482694068408160UL) + ((uint64_t)op[3] * 3032100978319218926UL) + ((((uint64_t)op[4] * 2239388948852210279UL) + ((uint64_t)op[5] * 15273176679450849295UL) + ((uint64_t)op[6] * 12273388804262012310UL) + ((uint64_t)op[7] * 11221737545519111142UL) + ((uint64_t)op[8] * 12201587493475465996UL) + ((uint64_t)op[9] * 13678276380802697737UL) + ((uint64_t)op[10] * 7051930127756445995UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 7051930127756445995UL) + ((uint64_t)op[1] * 5742670642601390004UL) + ((uint64_t)op[2] * 6888864087897378859UL) + ((uint64_t)op[3] * 13162482694068408160UL) + ((uint64_t)op[4] * 3032100978319218926UL) + ((((uint64_t)op[5] * 2239388948852210279UL) + ((uint64_t)op[6] * 15273176679450849295UL) + ((uint64_t)op[7] * 12273388804262012310UL) + ((uint64_t)op[8] * 11221737545519111142UL) + ((uint64_t)op[9] * 12201587493475465996UL) + ((uint64_t)op[10] * 13678276380802697737UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 13678276380802697737UL) + ((uint64_t)op[1] * 7051930127756445995UL) + ((uint64_t)op[2] * 5742670642601390004UL) + ((uint64_t)op[3] * 6888864087897378859UL) + ((uint64_t)op[4] * 13162482694068408160UL) + ((uint64_t)op[5] * 3032100978319218926UL) + ((((uint64_t)op[6] * 2239388948852210279UL) + ((uint64_t)op[7] * 15273176679450849295UL) + ((uint64_t)op[8] * 12273388804262012310UL) + ((uint64_t)op[9] * 11221737545519111142UL) + ((uint64_t)op[10] * 12201587493475465996UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 12201587493475465996UL) + ((uint64_t)op[1] * 13678276380802697737UL) + ((uint64_t)op[2] * 7051930127756445995UL) + ((uint64_t)op[3] * 5742670642601390004UL) + ((uint64_t)op[4] * 6888864087897378859UL) + ((uint64_t)op[5] * 13162482694068408160UL) + ((uint64_t)op[6] * 3032100978319218926UL) + ((((uint64_t)op[7] * 2239388948852210279UL) + ((uint64_t)op[8] * 15273176679450849295UL) + ((uint64_t)op[9] * 12273388804262012310UL) + ((uint64_t)op[10] * 11221737545519111142UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 11221737545519111142UL) + ((uint64_t)op[1] * 12201587493475465996UL) + ((uint64_t)op[2] * 13678276380802697737UL) + ((uint64_t)op[3] * 7051930127756445995UL) + ((uint64_t)op[4] * 5742670642601390004UL) + ((uint64_t)op[5] * 6888864087897378859UL) + ((uint64_t)op[6] * 13162482694068408160UL) + ((uint64_t)op[7] * 3032100978319218926UL) + ((((uint64_t)op[8] * 2239388948852210279UL) + ((uint64_t)op[9] * 15273176679450849295UL) + ((uint64_t)op[10] * 12273388804262012310UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 12273388804262012310UL) + ((uint64_t)op[1] * 11221737545519111142UL) + ((uint64_t)op[2] * 12201587493475465996UL) + ((uint64_t)op[3] * 13678276380802697737UL) + ((uint64_t)op[4] * 7051930127756445995UL) + ((uint64_t)op[5] * 5742670642601390004UL) + ((uint64_t)op[6] * 6888864087897378859UL) + ((uint64_t)op[7] * 13162482694068408160UL) + ((uint64_t)op[8] * 3032100978319218926UL) + ((((uint64_t)op[9] * 2239388948852210279UL) + ((uint64_t)op[10] * 15273176679450849295UL)) * 18446744073709551613);
	tmp_q[9] = ((uint64_t)op[0] * 15273176679450849295UL) + ((uint64_t)op[1] * 12273388804262012310UL) + ((uint64_t)op[2] * 11221737545519111142UL) + ((uint64_t)op[3] * 12201587493475465996UL) + ((uint64_t)op[4] * 13678276380802697737UL) + ((uint64_t)op[5] * 7051930127756445995UL) + ((uint64_t)op[6] * 5742670642601390004UL) + ((uint64_t)op[7] * 6888864087897378859UL) + ((uint64_t)op[8] * 13162482694068408160UL) + ((uint64_t)op[9] * 3032100978319218926UL) + ((uint64_t)op[10] * 11728577227152920779UL);
	tmp_q[10] = ((uint64_t)op[0] * 2239388948852210279UL) + ((uint64_t)op[1] * 15273176679450849295UL) + ((uint64_t)op[2] * 12273388804262012310UL) + ((uint64_t)op[3] * 11221737545519111142UL) + ((uint64_t)op[4] * 12201587493475465996UL) + ((uint64_t)op[5] * 13678276380802697737UL) + ((uint64_t)op[6] * 7051930127756445995UL) + ((uint64_t)op[7] * 5742670642601390004UL) + ((uint64_t)op[8] * 6888864087897378859UL) + ((uint64_t)op[9] * 13162482694068408160UL) + ((uint64_t)op[10] * 3032100978319218926UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 15383650762052L) - ((((int128)tmp_q[1] * 85774198029201L) + ((int128)tmp_q[2] * 5283654177453L) - ((int128)tmp_q[3] * 41540001577567L) - ((int128)tmp_q[4] * 41786313475016L) + ((int128)tmp_q[5] * 111114764157874L) + ((int128)tmp_q[6] * 59814511131620L) - ((int128)tmp_q[7] * 13329011321713L) + ((int128)tmp_q[8] * 66937362876471L) + ((int128)tmp_q[9] * 50357183814631L) - ((int128)tmp_q[10] * 1520523367665L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 1520523367665L) - ((int128)tmp_q[1] * 15383650762052L) - ((((int128)tmp_q[2] * 85774198029201L) + ((int128)tmp_q[3] * 5283654177453L) - ((int128)tmp_q[4] * 41540001577567L) - ((int128)tmp_q[5] * 41786313475016L) + ((int128)tmp_q[6] * 111114764157874L) + ((int128)tmp_q[7] * 59814511131620L) - ((int128)tmp_q[8] * 13329011321713L) + ((int128)tmp_q[9] * 66937362876471L) + ((int128)tmp_q[10] * 50357183814631L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 50357183814631L) - ((int128)tmp_q[1] * 1520523367665L) - ((int128)tmp_q[2] * 15383650762052L) - ((((int128)tmp_q[3] * 85774198029201L) + ((int128)tmp_q[4] * 5283654177453L) - ((int128)tmp_q[5] * 41540001577567L) - ((int128)tmp_q[6] * 41786313475016L) + ((int128)tmp_q[7] * 111114764157874L) + ((int128)tmp_q[8] * 59814511131620L) - ((int128)tmp_q[9] * 13329011321713L) + ((int128)tmp_q[10] * 66937362876471L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 66937362876471L) + ((int128)tmp_q[1] * 50357183814631L) - ((int128)tmp_q[2] * 1520523367665L) - ((int128)tmp_q[3] * 15383650762052L) - ((((int128)tmp_q[4] * 85774198029201L) + ((int128)tmp_q[5] * 5283654177453L) - ((int128)tmp_q[6] * 41540001577567L) - ((int128)tmp_q[7] * 41786313475016L) + ((int128)tmp_q[8] * 111114764157874L) + ((int128)tmp_q[9] * 59814511131620L) - ((int128)tmp_q[10] * 13329011321713L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 13329011321713L) + ((int128)tmp_q[1] * 66937362876471L) + ((int128)tmp_q[2] * 50357183814631L) - ((int128)tmp_q[3] * 1520523367665L) - ((int128)tmp_q[4] * 15383650762052L) - ((((int128)tmp_q[5] * 85774198029201L) + ((int128)tmp_q[6] * 5283654177453L) - ((int128)tmp_q[7] * 41540001577567L) - ((int128)tmp_q[8] * 41786313475016L) + ((int128)tmp_q[9] * 111114764157874L) + ((int128)tmp_q[10] * 59814511131620L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 59814511131620L) - ((int128)tmp_q[1] * 13329011321713L) + ((int128)tmp_q[2] * 66937362876471L) + ((int128)tmp_q[3] * 50357183814631L) - ((int128)tmp_q[4] * 1520523367665L) - ((int128)tmp_q[5] * 15383650762052L) - ((((int128)tmp_q[6] * 85774198029201L) + ((int128)tmp_q[7] * 5283654177453L) - ((int128)tmp_q[8] * 41540001577567L) - ((int128)tmp_q[9] * 41786313475016L) + ((int128)tmp_q[10] * 111114764157874L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 111114764157874L) + ((int128)tmp_q[1] * 59814511131620L) - ((int128)tmp_q[2] * 13329011321713L) + ((int128)tmp_q[3] * 66937362876471L) + ((int128)tmp_q[4] * 50357183814631L) - ((int128)tmp_q[5] * 1520523367665L) - ((int128)tmp_q[6] * 15383650762052L) - ((((int128)tmp_q[7] * 85774198029201L) + ((int128)tmp_q[8] * 5283654177453L) - ((int128)tmp_q[9] * 41540001577567L) - ((int128)tmp_q[10] * 41786313475016L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 41786313475016L) + ((int128)tmp_q[1] * 111114764157874L) + ((int128)tmp_q[2] * 59814511131620L) - ((int128)tmp_q[3] * 13329011321713L) + ((int128)tmp_q[4] * 66937362876471L) + ((int128)tmp_q[5] * 50357183814631L) - ((int128)tmp_q[6] * 1520523367665L) - ((int128)tmp_q[7] * 15383650762052L) - ((((int128)tmp_q[8] * 85774198029201L) + ((int128)tmp_q[9] * 5283654177453L) - ((int128)tmp_q[10] * 41540001577567L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 41540001577567L) - ((int128)tmp_q[1] * 41786313475016L) + ((int128)tmp_q[2] * 111114764157874L) + ((int128)tmp_q[3] * 59814511131620L) - ((int128)tmp_q[4] * 13329011321713L) + ((int128)tmp_q[5] * 66937362876471L) + ((int128)tmp_q[6] * 50357183814631L) - ((int128)tmp_q[7] * 1520523367665L) - ((int128)tmp_q[8] * 15383650762052L) - ((((int128)tmp_q[9] * 85774198029201L) + ((int128)tmp_q[10] * 5283654177453L)) * 3);
	tmp_zero[9] = ((int128)tmp_q[0] * 5283654177453L) - ((int128)tmp_q[1] * 41540001577567L) - ((int128)tmp_q[2] * 41786313475016L) + ((int128)tmp_q[3] * 111114764157874L) + ((int128)tmp_q[4] * 59814511131620L) - ((int128)tmp_q[5] * 13329011321713L) + ((int128)tmp_q[6] * 66937362876471L) + ((int128)tmp_q[7] * 50357183814631L) - ((int128)tmp_q[8] * 1520523367665L) - ((int128)tmp_q[9] * 15383650762052L) - ((int128)tmp_q[10] * 257322594087603L);
	tmp_zero[10] = ((int128)tmp_q[0] * 85774198029201L) + ((int128)tmp_q[1] * 5283654177453L) - ((int128)tmp_q[2] * 41540001577567L) - ((int128)tmp_q[3] * 41786313475016L) + ((int128)tmp_q[4] * 111114764157874L) + ((int128)tmp_q[5] * 59814511131620L) - ((int128)tmp_q[6] * 13329011321713L) + ((int128)tmp_q[7] * 66937362876471L) + ((int128)tmp_q[8] * 50357183814631L) - ((int128)tmp_q[9] * 1520523367665L) - ((int128)tmp_q[10] * 15383650762052L);

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

