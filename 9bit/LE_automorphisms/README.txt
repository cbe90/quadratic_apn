For each class of canonical LE-automorphisms, quad_apn.c performs a search for quadratic APN functions admitting that particular LE-automorphism. For lots of classes, this search terminates in short time (e.g., a few seconds up to minutes). However, for some classes, the deterministic search is infeasible. The program quad_apn_rand.c performs a randomized search for quadratic APN functions admitting that particular LE-automorphism.

After about 36 CPUh of computation for each of the classes, we found quadratic APN functions in Classes 3, 6, 8, 35, 49, and 82. 

Note that Classes 32, 34, 37-40, 50, 58, 59, 72, 74, 75, 78-81, 87-89, 91, 97-99, 102, 103, 107, 112-154, 156, 157 and all Classes >= 165 cannot yield APN functions, so they are not included in our search.

We did not perform any search for Classes 111 and 164, since B has >= 2^8 fixpoints and there are too many quadratic 8-bit APNs.
