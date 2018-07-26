# AMNS
This repository contains codes to generate Adapted Modular Number Systems (AMNS) for a given prime.
It also contains a C code generator for arithmetic (and conversion) operations using AMNS.  

The first subdirectory 'amns_generator' contains codes to generate AMNS given a prime and some parameters; see an example in file 'gen_example.py' of this subdirectory.

The second subdirectory 'c_code_generator' contains codes to generate a C code from an AMNS generated using codes in the subdirectory 'amns_generator'. This C code allows to perform arithmetic (and conversion) operations efficiently. The entry point of this C code generator is the file 'the_GENERATOR.sage' in the subdirectory 'c_code_generator'; this file contains an example of usage. The generated code will be in the subdirectory 'codes' of this subdirectory.
