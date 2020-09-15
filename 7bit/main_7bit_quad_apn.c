/* Compile with: g++ -O2 -march=native -funroll-loops main_7bit_quad_apn.c -o 7bit_quad_apn
Authors: C. Beierle, G. Leander  -- 2020

This program randomly searches for quadratic 7-bit APN functions.
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <random>

// dimension of the function to find
#define  N 7
#define RE_SHUFFLE 10

// the following three lines are needed for the degree bound
#define DEGREE 2 // we only search for permutations with algebraic degree less than or equal to DEGREE
int cube_ctr[(1<<N)];
int cube_sum[(1<<N)];

// to track the running time
time_t start;
time_t runtime;

// numbers of even Hamming weight up to 127 (needed for more efficient APN check)
int evens[((1<<(N-1))-1)] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126};

// the Hamming weights of range(1<<N). The lookup is slighly faster than __builtin_popcount
int hw[(1<<N)] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6};

int P[(1<<N)][(1<<N)];	// contains the random initial state

int solutions;
long long iterations;
int proceed;

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

// random number generators to generate integers uniformly in [0,(1<<N))
std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_int_distribution<int> dis(0,(1<<N)-1);

// applys the Fisher-Yates shuffle to the array arr to get a uniform permutation
void shuffle_array(int* arr) {
    uint32_t j = 0;
    uint32_t tmp = 0;
    for (int i=(1<<N)-1; i!=0; i--) {
        decltype(dis.param()) nrange (0, i);
        dis.param(nrange);
        j = dis(gen);
        // swap
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

void printArray(int A[1<<N]) {
	for (int i=0;i<(1<<N);i++)
		printf("%2x ",i);
	printf("\n");
	for (int i=0;i<(1<<N);i++)
		if (A[i]!=-1) printf("%2x ",A[i]);
			else printf("-- ");
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

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<((1<<(N-1))-1);i++) {
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
	for (int i=0;i<((1<<(N-1))-1);i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) return 0; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
  return 1;
}

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
	int y;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=RE_SHUFFLE) { // display the status every RE_SHUFFLE seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli \n",iterations);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
		proceed = 0;
	}

	//complete
	if (depth==((1<<N)-1)) {
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
            proceed = 0;
		return;
	}

	//not complete
	int x=depth+1;	// get next free position in the look-up table
	for (int z=0;z<(1<<N);z++) {
		y = P[depth][z];
			sbox[x]=y;
			if (!addPoint(x)) goto UNDO;

      if (proceed) nextValue(depth+1, filename);

			//undo the changes
			UNDO:
			removePoint(x);
			sbox[x]=-1;
	}
}

void search(char* filename, int max_runtime) {
    solutions=0;
    iterations=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);


		for (int t=0;t<(1<<N);t++) {
		 for (int x=0;x<(1<<N);x++) {
	 P[t][x] = x;
		 }
	 }

		 while(time(NULL)-runtime < max_runtime) {
	 memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
	 memset(cube_ctr,0,(1<<N)*sizeof(int));
	 memset(cube_sum,0,(1<<N)*sizeof(int));
	 for (int x=0;x<(1<<N);x++) {
					 sbox[x]=-1;
			 }
	 sbox[0] = 0;
	 if (!addPoint(0)) exit(0);

			 // shuffle the state
			 for (int t=0;t<(1<<N);t++) {
					shuffle_array(P[t]);
		 }
			 printf("Re-shuffled intitial state ..\n");
			 fflush(stdout);

			 // start the recursive search
			 proceed = 1;
			 nextValue(0,filename);
		 }
	}


	int main(int argc, char* argv[])
	{
		if (argc != 3) {
	        	printf("Usage: %s <outfile> <max_runtime> ...\n", argv[0]);
	        	return -1;
	    	}
		srand (time(NULL));
		start = time(NULL);
		runtime = time(NULL);

		search(argv[1],atoi(argv[2]));
		printf("\nFinished. Total solutions: %d. Total running time: %li sec\n", solutions, time(NULL)-runtime);
		exit(0);
	}
