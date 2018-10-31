#ifndef POLY_MULT_ADD
#define POLY_MULT_ADD


void sub_poly(int *rop, int *pa, int *pb);
void add_poly(int *rop, int *pa, int *pb);
void neg_poly(int *rop, int *op);
void scalar_mult_poly(int *rop, int *op, int scalar);

void mult_mod_poly(int *rop, int *pa, int *pb);
void square_mod_poly(int *rop, int *pa);

void internal_reduction(int *rop, llong *op);

#endif

