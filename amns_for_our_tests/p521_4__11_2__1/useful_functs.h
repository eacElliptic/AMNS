#ifndef USEFUL_FUNCTS
#define USEFUL_FUNCTS


void from_int_to_amns(int64_t *rop, mpz_t op);

void from_amns_to_int(mpz_t rop, int64_t *op);


int cmp_polys(int64_t *pa, int64_t *pb);

void copy_poly(int64_t *rop, int64_t *op);


void add_lpoly(int128 *rop, int128 *pa, int128 *pb);

void scalar_mult_lpoly(int128 *rop, int64_t *op, uint64_t scalar);


void from_mont_domain(int64_t *rop, int64_t *op);


void print_element(int64_t *poly);

#endif

