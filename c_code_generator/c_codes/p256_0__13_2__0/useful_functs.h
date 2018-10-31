#ifndef USEFUL_FUNCTS
#define USEFUL_FUNCTS


void from_int_to_amns(int *rop, mpz_t op);

void from_amns_to_int(mpz_t rop, int *op);


int cmp_polys(int *pa, int *pb);

void copy_poly(int *rop, int *op);


void add_lpoly(llong *rop, llong *pa, llong *pb);

void scalar_mult_lpoly(llong *rop, int *op, uint scalar);


void from_mont_domain(int *rop, int *op);


void print_element(int *poly);

#endif

