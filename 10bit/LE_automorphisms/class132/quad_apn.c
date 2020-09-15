/* Compile with: g++ -O2 -march=native -funroll-loops quad_apn.c -o quad_apn
Authors: C. Beierle, G. Leander  -- 2020

This program searches for all N-bit quadratic APN functions F s.t. AF = FB with A, B corresponding to the tuple defined in the header file "matrix_tuple.h", up to EA-equivalence
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "matrix_tuple.h"

// the following three lines are needed for the degree bound
#define DEGREE 2 // we only search for functions with algebraic degree less than or equal to DEGREE
int cube_ctr[(1<<N)];
int cube_sum[(1<<N)];

// to track running time
time_t start;
time_t runtime;

int solutions;
long long iterations;
int n_skipped;
int skipped[128];

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

void printArray(int A[1<<N]) {
	for (int i=0;i<(1<<N);i++)
		printf("%2x ",i);
	printf("\n");
	for (int i=0;i<(1<<N);i++)
		if (A[i]!=-1) printf("%2x ",A[i]);
			else printf("-- ");
	printf("\n\n");
}

void printSkippedArray(int A[128]) {
	for (int i=0;i<32;i++)
		printf("%d ",A[i]);
	printf("\n\n");
}

void printArray2(int A[1<<N]) {
	printf("{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		printf(",0x%02x",A[i]);
	printf("};\n");
}

void fprintArray2(FILE *fp, int A[1<<N]) {
	fprintf(fp, "{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		fprintf(fp, ",0x%02x",A[i]);
	fprintf(fp, "};\n");
}

// checks if $c \preceq x$
int is_prec(int c, int x) {
  int mask = 1;
  for (int i=0; i<N; i++) {
    if ((c&mask) > (x&mask)) return 0;
    mask = (mask << 1);
  }
  return 1;
}

// multiplicative order of x for matrix
int order_mat(int A[1<<N], int x) {
	int y=A[x];
	int i=1;
	while(y!=x) {
		y=A[y];
		i++;
	}
	return i;
}

// returns 1 if A is the smallest function up to conjugation with matrices in commuting_matrices. We use the order of the lookup tables defined in lex.
int is_smallest_in_class(int A[1<<N]){
    int E;
    for (int i=0;i<N_COMMUTING_MATRICES_B;i++) {
        for (int j=0; j<N_COMMUTING_MATRICES_A;j++) {
            for (int x=0;x<(1<<N);x++){
                if (A[lex[x]]==-1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                if (A[commuting_matrices_B[i][lex[x]]] == -1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                else {
                    E = commuting_matrices_A[j][A[commuting_matrices_B[i][lex[x]]]];
                    if (E < A[lex[x]]) return 0; // E is smaller, so A cannot be the smallest in class
                    if (E > A[lex[x]]) goto NEXT_MAT; // A is smaller, so consider next matrix
                } // consider next position
            }
            NEXT_MAT:;
        }
    }
    // A is the smallest of all, so return 1
    return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<(1<<(N-1))-1;i++) {
		a = evens[i];	// we only need to check those input differences corresponding to vectors with even Hamming Weight
                // if sbox[c + a] is set, increase DDT entry by 2
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]+=2;
                    // if set above 2, return that it cannot be APN anymore
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]>2) return 0;
                }

	}
	return 1;
}

// reducing DDT entries for new removed point sbox[c]. It reduces all entries that were increased during the APN check of addDDTInformation(c).
int removeDDTInformation(int c) {
	int a;
	for (int i=0;i<(1<<(N-1))-1;i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) return 0; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
  return 1;
}

// checks whether the ANF of sbox does not contain monomials of degree larger than DEGREE
int addDegreeInformation(int c) {
  for (int u=0; u<(1<<N); u++) {
    if (hw[u]>=(DEGREE+1) && is_prec(c,u)) {
      cube_ctr[u]++;
      cube_sum[u] = cube_sum[u]^sbox[c];
      if (cube_ctr[u] == (1<<hw[u])) {
        // if the sum over the cube is 0, the degree must be higher
        if (cube_sum[u]!=0) return 0;
      }
    }
  }
  return 1;
}

int removeDegreeInformation(int c) {
  for (int u=0; u<(1<<N); u++) {
    if (hw[u]>=(DEGREE+1) && is_prec(c,u)) {
      cube_ctr[u]--;
      cube_sum[u] = cube_sum[u]^sbox[c];
      if (cube_ctr[u] == (1<<hw[u])-1) {
        if (cube_sum[u]!=sbox[c]) return 0;
      }
    }
  }
  return 1;
}

int addPoint(int c) {
  // first check if the sbox can be apn
  if (addDDTInformation(c)) {
    // if there is no contradiction, check the algebraic degree
    return addDegreeInformation(c);
  }
  return 0;
}

void removePoint(int c) {
  if (removeDDTInformation(c)) {
    removeDegreeInformation(c);
  }
}

int maxDepth=0;

// recursive construction of the sbox.
void nextValue(int depth, char* filename) {
	int xS,yS;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=60) { // display the status every 60 seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli  skipped: %d \nDepth skipped: \n",iterations,n_skipped);
		printSkippedArray(skipped);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
	}

	//complete
	if (depth==N_POS) {
			solutions++;
			printf("found a new apn function: #%d\n",solutions);
			printArray2(sbox);
			fflush(stdout);
			FILE *fp = fopen(filename, "a");
			if (fp == NULL)
            {
                printf("Error opening file!\n");
		fflush(stdout);
            }
            else
            {
                fprintArray2(fp, sbox);
            }
            fclose(fp);
		return;
	}

	//not complete

	int x=pos[depth];	// get next free position in the look-up table
	for (int y=0;y<(1<<N);y++) {
		if (order_mat(B,x)>=order_mat(A,y)) {
			xS = x;
			yS = y;
			for (int i=0; i<order_mat(B,x)-1; i++) {
				sbox[xS]=yS;
				if (!addPoint(xS)) goto UNDO; // undo if it can't be APN or quadratic anymore
				xS=B[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addPoint(xS)) goto UNDO; // undo if it can't be APN or quadratic anymore

							if (depth <= 7) {
								if (is_smallest_in_class(sbox)) nextValue(depth+1, filename);
								else {n_skipped++; skipped[depth]++;}
							}
							else {	// omit check for smallest representative
										nextValue(depth+1, filename);
							}

			UNDO:
			// we have to remove the points in reverse order
			do {
				removePoint(xS);
				sbox[xS] = -1;
				xS = B_inv[xS];
			} while (xS != B_inv[x]);
		}
	}
}

void test_for_matrix(char* filename) {
    solutions=0;
    iterations=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrix: \n");
    printArray2(A);
    printArray2(B);
    printf("\n\n\n");
    fflush(stdout);

    for (int a=0;a<N_APN_FP;a++) {
	memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
  memset(cube_ctr,0,(1<<N)*sizeof(int));
  memset(cube_sum,0,(1<<N)*sizeof(int));
	for (int x=0;x<(1<<N);x++) {
        	sbox[x]=-1;
    	}

      // fix a quadratic apn on the fixpoints of B
      printf("Set new apn on the fixpoints.. \n");
      fflush(stdout);
      for (int x=0; x<N_FIXPOINTS_B; x++) {
        sbox[pos_fix_B[x]] = pos_fix_A[APN[a][x]];
        if (!addPoint(pos_fix_B[x])) exit(0);
      }

    	// start the recursive search
    	nextValue(0,filename);
  }
}


int main(int argc, char* argv[])
{
	if (argc != 2) {
        	printf("Usage: %s <outfile> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);

	test_for_matrix(argv[1]);
	printf("\nFinished. Total solutions: %d. Total running time: %li sec\n", solutions, time(NULL)-runtime);
	exit(0);
}
