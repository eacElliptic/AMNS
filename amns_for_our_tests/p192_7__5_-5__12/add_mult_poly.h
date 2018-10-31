#ifndef POLY_MULT_ADD
#define POLY_MULT_ADD


void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void add_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void neg_poly(int64_t *rop, int64_t *op);
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar);

void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb);
void square_mod_poly(int64_t *rop, int64_t *pa);

void internal_reduction(int64_t *rop, int128 *op);

#endif

