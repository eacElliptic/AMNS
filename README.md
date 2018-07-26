# AMNS
This repository contains codes to generate Adapted Modular Number Systems (AMNS) for a given prime.
It also contains a C code generator for arithmetic (and conversion) operations using AMNS.  

The subdirectory 'amns_generator' contains codes to generate AMNS given a prime and some parameters; see an example in file 'gen_example.py' of this subdirectory.

The subdirectory 'c_code_generator' contains codes to generate a C code from an AMNS generated using codes in the subdirectory 'amns_generator'. This C code allows to perform arithmetic (and conversion) operations efficiently. The entry point of this C code generator is the file 'the_GENERATOR.sage' in the subdirectory 'c_code_generator'; this file contains an example of usage. The generated code will be in the subdirectory 'codes' of this subdirectory.

The subdirectory 'amns_c_codes_for_our_tests' contains the C codes of the AMNS we used to build the table of performances in the article. The subdirectory 'ns_p521' contains the code for the recommended NIST prime 'p = 2^521 - 1' for elliptic curve cryptography. The others subdirectories 'pX' contain C codes for the X-bits primes we used in the article.
Notice that in each subdirectory, there is a file 'test.txt' which contains the results we used to build the table of performances.
