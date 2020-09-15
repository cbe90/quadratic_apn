For each class of canonical LE-automorphisms, quad_apn.c performs a search for quadratic APN functions admitting that particular LE-automorphism. For lots of classes, this search terminates in short time (e.g., a few seconds up to minutes). However, for some classes, the deterministic search is infeasible. The program quad_apn_rand.c performs a randomized search for quadratic APN functions admitting that particular LE-automorphism.

We found quadratic APN functions in Classes 5, 47, 49, 75, 78, 163, and 167. Note that we know that solutions must exist in Classes 140 and 218, since those classes correspond to polynomials with coefficients in a subfield of GF(2^10). However, after running the randomized search for roughly 16 CPU days for both of those classes, we still did not find solutions. 

Note that Classes 10-12, 39-44, 66-71, 101-103, 106, 107, 112, 141, 157-162, 178, 179, 182-184, 187, 188, 200, 201, 203, 205, 206, 211, 219-226, 230-232, 236, 237, 242, 248-314, and all Classes >= 325 cannot yield APN functions, so they are not included in our search.

We did not perform any search for Classes 216, 241, 246, 247, 318, 323, and 324, since B has >= 2^8 fixpoints and there are too many quadratic 8-bit APNs.
