load("GENERATOR_amns.sage")


word_size = 64

p_size = 256

n_min = (p_size//word_size) + 1 
n_max = n_min + 2

abs_lamb_max = 1 << 3

nth_root_max_duration_checks = 60  # seconds

find_all_nthroot = True


p = random_prime(2**p_size, lbound=2**(p_size-1))
print (p)
print (p.is_prime())
print (p.nbits())

	
flt_amns = generate_amns_candidates_with_n_min_max(word_size, p, n_min, n_max, abs_lamb_max, nth_root_max_duration_checks, find_all_nthroot)

print (len(flt_amns))

