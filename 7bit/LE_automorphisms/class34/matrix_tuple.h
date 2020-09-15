#define N 7
#define N_FIXPOINTS_A 16
#define N_FIXPOINTS_B 2
#define N_POS 18
#define N_COMMUTING_MATRICES_A 7
#define N_COMMUTING_MATRICES_B 49
#define N_APN_FP 1

int A[(1<<N)] = {0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 64, 65, 66, 67, 68, 69, 70, 71, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 112, 113, 114, 115, 116, 117, 118, 119, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 72, 73, 74, 75, 76, 77, 78, 79, 88, 89, 90, 91, 92, 93, 94, 95, 40, 41, 42, 43, 44, 45, 46, 47, 56, 57, 58, 59, 60, 61, 62, 63, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31};

int B[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127};

int B_inv[(1<<N)] = {0, 64, 1, 65, 2, 66, 3, 67, 4, 68, 5, 69, 6, 70, 7, 71, 8, 72, 9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79, 16, 80, 17, 81, 18, 82, 19, 83, 20, 84, 21, 85, 22, 86, 23, 87, 24, 88, 25, 89, 26, 90, 27, 91, 28, 92, 29, 93, 30, 94, 31, 95, 32, 96, 33, 97, 34, 98, 35, 99, 36, 100, 37, 101, 38, 102, 39, 103, 40, 104, 41, 105, 42, 106, 43, 107, 44, 108, 45, 109, 46, 110, 47, 111, 48, 112, 49, 113, 50, 114, 51, 115, 52, 116, 53, 117, 54, 118, 55, 119, 56, 120, 57, 121, 58, 122, 59, 123, 60, 124, 61, 125, 62, 126, 63, 127};

int pos[N_POS] = {1, 3, 5, 7, 9, 11, 13, 15, 19, 21, 23, 27, 29, 31, 43, 47, 55, 63};

int pos_fix_A[N_FIXPOINTS_A] = {0, 1, 2, 3, 4, 5, 6, 7, 88, 89, 90, 91, 92, 93, 94, 95};

int pos_fix_B[N_FIXPOINTS_B] = {0, 127};

int lex[(1<<N)] = {0, 127, 1, 2, 4, 8, 16, 32, 64, 3, 6, 12, 24, 48, 96, 65, 5, 10, 20, 40, 80, 33, 66, 7, 14, 28, 56, 112, 97, 67, 9, 18, 36, 72, 17, 34, 68, 11, 22, 44, 88, 49, 98, 69, 13, 26, 52, 104, 81, 35, 70, 15, 30, 60, 120, 113, 99, 71, 19, 38, 76, 25, 50, 100, 73, 21, 42, 84, 41, 82, 37, 74, 23, 46, 92, 57, 114, 101, 75, 27, 54, 108, 89, 51, 102, 77, 29, 58, 116, 105, 83, 39, 78, 31, 62, 124, 121, 115, 103, 79, 43, 86, 45, 90, 53, 106, 85, 47, 94, 61, 122, 117, 107, 87, 55, 110, 93, 59, 118, 109, 91, 63, 126, 125, 123, 119, 111, 95};

int evens[(1<<(N-1))-1] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126};

int hw[(1<<N)] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7};

int APN[N_APN_FP][N_FIXPOINTS_B] = {{0, 0}};

int commuting_matrices_A[N_COMMUTING_MATRICES_A][(1<<N)] = {
{0, 1, 2, 3, 4, 5, 6, 7, 112, 113, 114, 115, 116, 117, 118, 119, 8, 9, 10, 11, 12, 13, 14, 15, 120, 121, 122, 123, 124, 125, 126, 127, 16, 17, 18, 19, 20, 21, 22, 23, 96, 97, 98, 99, 100, 101, 102, 103, 24, 25, 26, 27, 28, 29, 30, 31, 104, 105, 106, 107, 108, 109, 110, 111, 32, 33, 34, 35, 36, 37, 38, 39, 80, 81, 82, 83, 84, 85, 86, 87, 40, 41, 42, 43, 44, 45, 46, 47, 88, 89, 90, 91, 92, 93, 94, 95, 48, 49, 50, 51, 52, 53, 54, 55, 64, 65, 66, 67, 68, 69, 70, 71, 56, 57, 58, 59, 60, 61, 62, 63, 72, 73, 74, 75, 76, 77, 78, 79},
{0, 1, 2, 3, 4, 5, 6, 7, 104, 105, 106, 107, 108, 109, 110, 111, 56, 57, 58, 59, 60, 61, 62, 63, 80, 81, 82, 83, 84, 85, 86, 87, 112, 113, 114, 115, 116, 117, 118, 119, 24, 25, 26, 27, 28, 29, 30, 31, 72, 73, 74, 75, 76, 77, 78, 79, 32, 33, 34, 35, 36, 37, 38, 39, 8, 9, 10, 11, 12, 13, 14, 15, 96, 97, 98, 99, 100, 101, 102, 103, 48, 49, 50, 51, 52, 53, 54, 55, 88, 89, 90, 91, 92, 93, 94, 95, 120, 121, 122, 123, 124, 125, 126, 127, 16, 17, 18, 19, 20, 21, 22, 23, 64, 65, 66, 67, 68, 69, 70, 71, 40, 41, 42, 43, 44, 45, 46, 47},
{0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 64, 65, 66, 67, 68, 69, 70, 71, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 112, 113, 114, 115, 116, 117, 118, 119, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 72, 73, 74, 75, 76, 77, 78, 79, 88, 89, 90, 91, 92, 93, 94, 95, 40, 41, 42, 43, 44, 45, 46, 47, 56, 57, 58, 59, 60, 61, 62, 63, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31},
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127},
{0, 1, 2, 3, 4, 5, 6, 7, 64, 65, 66, 67, 68, 69, 70, 71, 104, 105, 106, 107, 108, 109, 110, 111, 40, 41, 42, 43, 44, 45, 46, 47, 56, 57, 58, 59, 60, 61, 62, 63, 120, 121, 122, 123, 124, 125, 126, 127, 80, 81, 82, 83, 84, 85, 86, 87, 16, 17, 18, 19, 20, 21, 22, 23, 112, 113, 114, 115, 116, 117, 118, 119, 48, 49, 50, 51, 52, 53, 54, 55, 24, 25, 26, 27, 28, 29, 30, 31, 88, 89, 90, 91, 92, 93, 94, 95, 72, 73, 74, 75, 76, 77, 78, 79, 8, 9, 10, 11, 12, 13, 14, 15, 32, 33, 34, 35, 36, 37, 38, 39, 96, 97, 98, 99, 100, 101, 102, 103},
{0, 1, 2, 3, 4, 5, 6, 7, 32, 33, 34, 35, 36, 37, 38, 39, 64, 65, 66, 67, 68, 69, 70, 71, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 72, 73, 74, 75, 76, 77, 78, 79, 40, 41, 42, 43, 44, 45, 46, 47, 8, 9, 10, 11, 12, 13, 14, 15, 56, 57, 58, 59, 60, 61, 62, 63, 24, 25, 26, 27, 28, 29, 30, 31, 120, 121, 122, 123, 124, 125, 126, 127, 88, 89, 90, 91, 92, 93, 94, 95, 80, 81, 82, 83, 84, 85, 86, 87, 112, 113, 114, 115, 116, 117, 118, 119, 16, 17, 18, 19, 20, 21, 22, 23, 48, 49, 50, 51, 52, 53, 54, 55},
{0, 1, 2, 3, 4, 5, 6, 7, 56, 57, 58, 59, 60, 61, 62, 63, 112, 113, 114, 115, 116, 117, 118, 119, 72, 73, 74, 75, 76, 77, 78, 79, 8, 9, 10, 11, 12, 13, 14, 15, 48, 49, 50, 51, 52, 53, 54, 55, 120, 121, 122, 123, 124, 125, 126, 127, 64, 65, 66, 67, 68, 69, 70, 71, 16, 17, 18, 19, 20, 21, 22, 23, 40, 41, 42, 43, 44, 45, 46, 47, 96, 97, 98, 99, 100, 101, 102, 103, 88, 89, 90, 91, 92, 93, 94, 95, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 104, 105, 106, 107, 108, 109, 110, 111, 80, 81, 82, 83, 84, 85, 86, 87}
};

int commuting_matrices_B[N_COMMUTING_MATRICES_B][(1<<N)] = {
{0, 64, 1, 65, 2, 66, 3, 67, 4, 68, 5, 69, 6, 70, 7, 71, 8, 72, 9, 73, 10, 74, 11, 75, 12, 76, 13, 77, 14, 78, 15, 79, 16, 80, 17, 81, 18, 82, 19, 83, 20, 84, 21, 85, 22, 86, 23, 87, 24, 88, 25, 89, 26, 90, 27, 91, 28, 92, 29, 93, 30, 94, 31, 95, 32, 96, 33, 97, 34, 98, 35, 99, 36, 100, 37, 101, 38, 102, 39, 103, 40, 104, 41, 105, 42, 106, 43, 107, 44, 108, 45, 109, 46, 110, 47, 111, 48, 112, 49, 113, 50, 114, 51, 115, 52, 116, 53, 117, 54, 118, 55, 119, 56, 120, 57, 121, 58, 122, 59, 123, 60, 124, 61, 125, 62, 126, 63, 127},
{0, 67, 7, 68, 14, 77, 9, 74, 28, 95, 27, 88, 18, 81, 21, 86, 56, 123, 63, 124, 54, 117, 49, 114, 36, 103, 35, 96, 42, 105, 45, 110, 112, 51, 119, 52, 126, 61, 121, 58, 108, 47, 107, 40, 98, 33, 101, 38, 72, 11, 79, 12, 70, 5, 65, 2, 84, 23, 83, 16, 90, 25, 93, 30, 97, 34, 102, 37, 111, 44, 104, 43, 125, 62, 122, 57, 115, 48, 116, 55, 89, 26, 94, 29, 87, 20, 80, 19, 69, 6, 66, 1, 75, 8, 76, 15, 17, 82, 22, 85, 31, 92, 24, 91, 13, 78, 10, 73, 3, 64, 4, 71, 41, 106, 46, 109, 39, 100, 32, 99, 53, 118, 50, 113, 59, 120, 60, 127},
{0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127},
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127},
{0, 97, 67, 34, 7, 102, 68, 37, 14, 111, 77, 44, 9, 104, 74, 43, 28, 125, 95, 62, 27, 122, 88, 57, 18, 115, 81, 48, 21, 116, 86, 55, 56, 89, 123, 26, 63, 94, 124, 29, 54, 87, 117, 20, 49, 80, 114, 19, 36, 69, 103, 6, 35, 66, 96, 1, 42, 75, 105, 8, 45, 76, 110, 15, 112, 17, 51, 82, 119, 22, 52, 85, 126, 31, 61, 92, 121, 24, 58, 91, 108, 13, 47, 78, 107, 10, 40, 73, 98, 3, 33, 64, 101, 4, 38, 71, 72, 41, 11, 106, 79, 46, 12, 109, 70, 39, 5, 100, 65, 32, 2, 99, 84, 53, 23, 118, 83, 50, 16, 113, 90, 59, 25, 120, 93, 60, 30, 127},
{0, 32, 64, 96, 1, 33, 65, 97, 2, 34, 66, 98, 3, 35, 67, 99, 4, 36, 68, 100, 5, 37, 69, 101, 6, 38, 70, 102, 7, 39, 71, 103, 8, 40, 72, 104, 9, 41, 73, 105, 10, 42, 74, 106, 11, 43, 75, 107, 12, 44, 76, 108, 13, 45, 77, 109, 14, 46, 78, 110, 15, 47, 79, 111, 16, 48, 80, 112, 17, 49, 81, 113, 18, 50, 82, 114, 19, 51, 83, 115, 20, 52, 84, 116, 21, 53, 85, 117, 22, 54, 86, 118, 23, 55, 87, 119, 24, 56, 88, 120, 25, 57, 89, 121, 26, 58, 90, 122, 27, 59, 91, 123, 28, 60, 92, 124, 29, 61, 93, 125, 30, 62, 94, 126, 31, 63, 95, 127},
{0, 82, 37, 119, 74, 24, 111, 61, 21, 71, 48, 98, 95, 13, 122, 40, 42, 120, 15, 93, 96, 50, 69, 23, 63, 109, 26, 72, 117, 39, 80, 2, 84, 6, 113, 35, 30, 76, 59, 105, 65, 19, 100, 54, 11, 89, 46, 124, 126, 44, 91, 9, 52, 102, 17, 67, 107, 57, 78, 28, 33, 115, 4, 86, 41, 123, 12, 94, 99, 49, 70, 20, 60, 110, 25, 75, 118, 36, 83, 1, 3, 81, 38, 116, 73, 27, 108, 62, 22, 68, 51, 97, 92, 14, 121, 43, 125, 47, 88, 10, 55, 101, 18, 64, 104, 58, 77, 31, 34, 112, 7, 85, 87, 5, 114, 32, 29, 79, 56, 106, 66, 16, 103, 53, 8, 90, 45, 127},
{0, 16, 32, 48, 64, 80, 96, 112, 1, 17, 33, 49, 65, 81, 97, 113, 2, 18, 34, 50, 66, 82, 98, 114, 3, 19, 35, 51, 67, 83, 99, 115, 4, 20, 36, 52, 68, 84, 100, 116, 5, 21, 37, 53, 69, 85, 101, 117, 6, 22, 38, 54, 70, 86, 102, 118, 7, 23, 39, 55, 71, 87, 103, 119, 8, 24, 40, 56, 72, 88, 104, 120, 9, 25, 41, 57, 73, 89, 105, 121, 10, 26, 42, 58, 74, 90, 106, 122, 11, 27, 43, 59, 75, 91, 107, 123, 12, 28, 44, 60, 76, 92, 108, 124, 13, 29, 45, 61, 77, 93, 109, 125, 14, 30, 46, 62, 78, 94, 110, 126, 15, 31, 47, 63, 79, 95, 111, 127},
{0, 19, 38, 53, 76, 95, 106, 121, 25, 10, 63, 44, 85, 70, 115, 96, 50, 33, 20, 7, 126, 109, 88, 75, 43, 56, 13, 30, 103, 116, 65, 82, 100, 119, 66, 81, 40, 59, 14, 29, 125, 110, 91, 72, 49, 34, 23, 4, 86, 69, 112, 99, 26, 9, 60, 47, 79, 92, 105, 122, 3, 16, 37, 54, 73, 90, 111, 124, 5, 22, 35, 48, 80, 67, 118, 101, 28, 15, 58, 41, 123, 104, 93, 78, 55, 36, 17, 2, 98, 113, 68, 87, 46, 61, 8, 27, 45, 62, 11, 24, 97, 114, 71, 84, 52, 39, 18, 1, 120, 107, 94, 77, 31, 12, 57, 42, 83, 64, 117, 102, 6, 21, 32, 51, 74, 89, 108, 127},
{0, 112, 97, 17, 67, 51, 34, 82, 7, 119, 102, 22, 68, 52, 37, 85, 14, 126, 111, 31, 77, 61, 44, 92, 9, 121, 104, 24, 74, 58, 43, 91, 28, 108, 125, 13, 95, 47, 62, 78, 27, 107, 122, 10, 88, 40, 57, 73, 18, 98, 115, 3, 81, 33, 48, 64, 21, 101, 116, 4, 86, 38, 55, 71, 56, 72, 89, 41, 123, 11, 26, 106, 63, 79, 94, 46, 124, 12, 29, 109, 54, 70, 87, 39, 117, 5, 20, 100, 49, 65, 80, 32, 114, 2, 19, 99, 36, 84, 69, 53, 103, 23, 6, 118, 35, 83, 66, 50, 96, 16, 1, 113, 42, 90, 75, 59, 105, 25, 8, 120, 45, 93, 76, 60, 110, 30, 15, 127},
{0, 115, 103, 20, 79, 60, 40, 91, 31, 108, 120, 11, 80, 35, 55, 68, 62, 77, 89, 42, 113, 2, 22, 101, 33, 82, 70, 53, 110, 29, 9, 122, 124, 15, 27, 104, 51, 64, 84, 39, 99, 16, 4, 119, 44, 95, 75, 56, 66, 49, 37, 86, 13, 126, 106, 25, 93, 46, 58, 73, 18, 97, 117, 6, 121, 10, 30, 109, 54, 69, 81, 34, 102, 21, 1, 114, 41, 90, 78, 61, 71, 52, 32, 83, 8, 123, 111, 28, 88, 43, 63, 76, 23, 100, 112, 3, 5, 118, 98, 17, 74, 57, 45, 94, 26, 105, 125, 14, 85, 38, 50, 65, 59, 72, 92, 47, 116, 7, 19, 96, 36, 87, 67, 48, 107, 24, 12, 127},
{0, 50, 100, 86, 73, 123, 45, 31, 19, 33, 119, 69, 90, 104, 62, 12, 38, 20, 66, 112, 111, 93, 11, 57, 53, 7, 81, 99, 124, 78, 24, 42, 76, 126, 40, 26, 5, 55, 97, 83, 95, 109, 59, 9, 22, 36, 114, 64, 106, 88, 14, 60, 35, 17, 71, 117, 121, 75, 29, 47, 48, 2, 84, 102, 25, 43, 125, 79, 80, 98, 52, 6, 10, 56, 110, 92, 67, 113, 39, 21, 63, 13, 91, 105, 118, 68, 18, 32, 44, 30, 72, 122, 101, 87, 1, 51, 85, 103, 49, 3, 28, 46, 120, 74, 70, 116, 34, 16, 15, 61, 107, 89, 115, 65, 23, 37, 58, 8, 94, 108, 96, 82, 4, 54, 41, 27, 77, 127},
{0, 74, 21, 95, 42, 96, 63, 117, 84, 30, 65, 11, 126, 52, 107, 33, 41, 99, 60, 118, 3, 73, 22, 92, 125, 55, 104, 34, 87, 29, 66, 8, 82, 24, 71, 13, 120, 50, 109, 39, 6, 76, 19, 89, 44, 102, 57, 115, 123, 49, 110, 36, 81, 27, 68, 14, 47, 101, 58, 112, 5, 79, 16, 90, 37, 111, 48, 122, 15, 69, 26, 80, 113, 59, 100, 46, 91, 17, 78, 4, 12, 70, 25, 83, 38, 108, 51, 121, 88, 18, 77, 7, 114, 56, 103, 45, 119, 61, 98, 40, 93, 23, 72, 2, 35, 105, 54, 124, 9, 67, 28, 86, 94, 20, 75, 1, 116, 62, 97, 43, 10, 64, 31, 85, 32, 106, 53, 127},
{0, 73, 19, 90, 38, 111, 53, 124, 76, 5, 95, 22, 106, 35, 121, 48, 25, 80, 10, 67, 63, 118, 44, 101, 85, 28, 70, 15, 115, 58, 96, 41, 50, 123, 33, 104, 20, 93, 7, 78, 126, 55, 109, 36, 88, 17, 75, 2, 43, 98, 56, 113, 13, 68, 30, 87, 103, 46, 116, 61, 65, 8, 82, 27, 100, 45, 119, 62, 66, 11, 81, 24, 40, 97, 59, 114, 14, 71, 29, 84, 125, 52, 110, 39, 91, 18, 72, 1, 49, 120, 34, 107, 23, 94, 4, 77, 86, 31, 69, 12, 112, 57, 99, 42, 26, 83, 9, 64, 60, 117, 47, 102, 79, 6, 92, 21, 105, 32, 122, 51, 3, 74, 16, 89, 37, 108, 54, 127},
{0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 1, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121, 2, 10, 18, 26, 34, 42, 50, 58, 66, 74, 82, 90, 98, 106, 114, 122, 3, 11, 19, 27, 35, 43, 51, 59, 67, 75, 83, 91, 99, 107, 115, 123, 4, 12, 20, 28, 36, 44, 52, 60, 68, 76, 84, 92, 100, 108, 116, 124, 5, 13, 21, 29, 37, 45, 53, 61, 69, 77, 85, 93, 101, 109, 117, 125, 6, 14, 22, 30, 38, 46, 54, 62, 70, 78, 86, 94, 102, 110, 118, 126, 7, 15, 23, 31, 39, 47, 55, 63, 71, 79, 87, 95, 103, 111, 119, 127},
{0, 107, 87, 60, 47, 68, 120, 19, 94, 53, 9, 98, 113, 26, 38, 77, 61, 86, 106, 1, 18, 121, 69, 46, 99, 8, 52, 95, 76, 39, 27, 112, 122, 17, 45, 70, 85, 62, 2, 105, 36, 79, 115, 24, 11, 96, 92, 55, 71, 44, 16, 123, 104, 3, 63, 84, 25, 114, 78, 37, 54, 93, 97, 10, 117, 30, 34, 73, 90, 49, 13, 102, 43, 64, 124, 23, 4, 111, 83, 56, 72, 35, 31, 116, 103, 12, 48, 91, 22, 125, 65, 42, 57, 82, 110, 5, 15, 100, 88, 51, 32, 75, 119, 28, 81, 58, 6, 109, 126, 21, 41, 66, 50, 89, 101, 14, 29, 118, 74, 33, 108, 7, 59, 80, 67, 40, 20, 127},
{0, 42, 84, 126, 41, 3, 125, 87, 82, 120, 6, 44, 123, 81, 47, 5, 37, 15, 113, 91, 12, 38, 88, 114, 119, 93, 35, 9, 94, 116, 10, 32, 74, 96, 30, 52, 99, 73, 55, 29, 24, 50, 76, 102, 49, 27, 101, 79, 111, 69, 59, 17, 70, 108, 18, 56, 61, 23, 105, 67, 20, 62, 64, 106, 21, 63, 65, 107, 60, 22, 104, 66, 71, 109, 19, 57, 110, 68, 58, 16, 48, 26, 100, 78, 25, 51, 77, 103, 98, 72, 54, 28, 75, 97, 31, 53, 95, 117, 11, 33, 118, 92, 34, 8, 13, 39, 89, 115, 36, 14, 112, 90, 122, 80, 46, 4, 83, 121, 7, 45, 40, 2, 124, 86, 1, 43, 85, 127},
{0, 41, 82, 123, 37, 12, 119, 94, 74, 99, 24, 49, 111, 70, 61, 20, 21, 60, 71, 110, 48, 25, 98, 75, 95, 118, 13, 36, 122, 83, 40, 1, 42, 3, 120, 81, 15, 38, 93, 116, 96, 73, 50, 27, 69, 108, 23, 62, 63, 22, 109, 68, 26, 51, 72, 97, 117, 92, 39, 14, 80, 121, 2, 43, 84, 125, 6, 47, 113, 88, 35, 10, 30, 55, 76, 101, 59, 18, 105, 64, 65, 104, 19, 58, 100, 77, 54, 31, 11, 34, 89, 112, 46, 7, 124, 85, 126, 87, 44, 5, 91, 114, 9, 32, 52, 29, 102, 79, 17, 56, 67, 106, 107, 66, 57, 16, 78, 103, 28, 53, 33, 8, 115, 90, 4, 45, 86, 127},
{0, 91, 55, 108, 110, 53, 89, 2, 93, 6, 106, 49, 51, 104, 4, 95, 59, 96, 12, 87, 85, 14, 98, 57, 102, 61, 81, 10, 8, 83, 63, 100, 118, 45, 65, 26, 24, 67, 47, 116, 43, 112, 28, 71, 69, 30, 114, 41, 77, 22, 122, 33, 35, 120, 20, 79, 16, 75, 39, 124, 126, 37, 73, 18, 109, 54, 90, 1, 3, 88, 52, 111, 48, 107, 7, 92, 94, 5, 105, 50, 86, 13, 97, 58, 56, 99, 15, 84, 11, 80, 60, 103, 101, 62, 82, 9, 27, 64, 44, 119, 117, 46, 66, 25, 70, 29, 113, 42, 40, 115, 31, 68, 32, 123, 23, 76, 78, 21, 121, 34, 125, 38, 74, 17, 19, 72, 36, 127},
{0, 25, 50, 43, 100, 125, 86, 79, 73, 80, 123, 98, 45, 52, 31, 6, 19, 10, 33, 56, 119, 110, 69, 92, 90, 67, 104, 113, 62, 39, 12, 21, 38, 63, 20, 13, 66, 91, 112, 105, 111, 118, 93, 68, 11, 18, 57, 32, 53, 44, 7, 30, 81, 72, 99, 122, 124, 101, 78, 87, 24, 1, 42, 51, 76, 85, 126, 103, 40, 49, 26, 3, 5, 28, 55, 46, 97, 120, 83, 74, 95, 70, 109, 116, 59, 34, 9, 16, 22, 15, 36, 61, 114, 107, 64, 89, 106, 115, 88, 65, 14, 23, 60, 37, 35, 58, 17, 8, 71, 94, 117, 108, 121, 96, 75, 82, 29, 4, 47, 54, 48, 41, 2, 27, 84, 77, 102, 127},
{0, 122, 117, 15, 107, 17, 30, 100, 87, 45, 34, 88, 60, 70, 73, 51, 47, 85, 90, 32, 68, 62, 49, 75, 120, 2, 13, 119, 19, 105, 102, 28, 94, 36, 43, 81, 53, 79, 64, 58, 9, 115, 124, 6, 98, 24, 23, 109, 113, 11, 4, 126, 26, 96, 111, 21, 38, 92, 83, 41, 77, 55, 56, 66, 61, 71, 72, 50, 86, 44, 35, 89, 106, 16, 31, 101, 1, 123, 116, 14, 18, 104, 103, 29, 121, 3, 12, 118, 69, 63, 48, 74, 46, 84, 91, 33, 99, 25, 22, 108, 8, 114, 125, 7, 52, 78, 65, 59, 95, 37, 42, 80, 76, 54, 57, 67, 39, 93, 82, 40, 27, 97, 110, 20, 112, 10, 5, 127},
{0, 121, 115, 10, 103, 30, 20, 109, 79, 54, 60, 69, 40, 81, 91, 34, 31, 102, 108, 21, 120, 1, 11, 114, 80, 41, 35, 90, 55, 78, 68, 61, 62, 71, 77, 52, 89, 32, 42, 83, 113, 8, 2, 123, 22, 111, 101, 28, 33, 88, 82, 43, 70, 63, 53, 76, 110, 23, 29, 100, 9, 112, 122, 3, 124, 5, 15, 118, 27, 98, 104, 17, 51, 74, 64, 57, 84, 45, 39, 94, 99, 26, 16, 105, 4, 125, 119, 14, 44, 85, 95, 38, 75, 50, 56, 65, 66, 59, 49, 72, 37, 92, 86, 47, 13, 116, 126, 7, 106, 19, 25, 96, 93, 36, 46, 87, 58, 67, 73, 48, 18, 107, 97, 24, 117, 12, 6, 127},
{0, 56, 112, 72, 97, 89, 17, 41, 67, 123, 51, 11, 34, 26, 82, 106, 7, 63, 119, 79, 102, 94, 22, 46, 68, 124, 52, 12, 37, 29, 85, 109, 14, 54, 126, 70, 111, 87, 31, 39, 77, 117, 61, 5, 44, 20, 92, 100, 9, 49, 121, 65, 104, 80, 24, 32, 74, 114, 58, 2, 43, 19, 91, 99, 28, 36, 108, 84, 125, 69, 13, 53, 95, 103, 47, 23, 62, 6, 78, 118, 27, 35, 107, 83, 122, 66, 10, 50, 88, 96, 40, 16, 57, 1, 73, 113, 18, 42, 98, 90, 115, 75, 3, 59, 81, 105, 33, 25, 48, 8, 64, 120, 21, 45, 101, 93, 116, 76, 4, 60, 86, 110, 38, 30, 55, 15, 71, 127},
{0, 59, 118, 77, 109, 86, 27, 32, 91, 96, 45, 22, 54, 13, 64, 123, 55, 12, 65, 122, 90, 97, 44, 23, 108, 87, 26, 33, 1, 58, 119, 76, 110, 85, 24, 35, 3, 56, 117, 78, 53, 14, 67, 120, 88, 99, 46, 21, 89, 98, 47, 20, 52, 15, 66, 121, 2, 57, 116, 79, 111, 84, 25, 34, 93, 102, 43, 16, 48, 11, 70, 125, 6, 61, 112, 75, 107, 80, 29, 38, 106, 81, 28, 39, 7, 60, 113, 74, 49, 10, 71, 124, 92, 103, 42, 17, 51, 8, 69, 126, 94, 101, 40, 19, 104, 83, 30, 37, 5, 62, 115, 72, 4, 63, 114, 73, 105, 82, 31, 36, 95, 100, 41, 18, 50, 9, 68, 127},
{0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113, 117, 121, 125, 2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74, 78, 82, 86, 90, 94, 98, 102, 106, 110, 114, 118, 122, 126, 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127},
{0, 7, 14, 9, 28, 27, 18, 21, 56, 63, 54, 49, 36, 35, 42, 45, 112, 119, 126, 121, 108, 107, 98, 101, 72, 79, 70, 65, 84, 83, 90, 93, 97, 102, 111, 104, 125, 122, 115, 116, 89, 94, 87, 80, 69, 66, 75, 76, 17, 22, 31, 24, 13, 10, 3, 4, 41, 46, 39, 32, 53, 50, 59, 60, 67, 68, 77, 74, 95, 88, 81, 86, 123, 124, 117, 114, 103, 96, 105, 110, 51, 52, 61, 58, 47, 40, 33, 38, 11, 12, 5, 2, 23, 16, 25, 30, 34, 37, 44, 43, 62, 57, 48, 55, 26, 29, 20, 19, 6, 1, 8, 15, 82, 85, 92, 91, 78, 73, 64, 71, 106, 109, 100, 99, 118, 113, 120, 127},
{0, 100, 73, 45, 19, 119, 90, 62, 38, 66, 111, 11, 53, 81, 124, 24, 76, 40, 5, 97, 95, 59, 22, 114, 106, 14, 35, 71, 121, 29, 48, 84, 25, 125, 80, 52, 10, 110, 67, 39, 63, 91, 118, 18, 44, 72, 101, 1, 85, 49, 28, 120, 70, 34, 15, 107, 115, 23, 58, 94, 96, 4, 41, 77, 50, 86, 123, 31, 33, 69, 104, 12, 20, 112, 93, 57, 7, 99, 78, 42, 126, 26, 55, 83, 109, 9, 36, 64, 88, 60, 17, 117, 75, 47, 2, 102, 43, 79, 98, 6, 56, 92, 113, 21, 13, 105, 68, 32, 30, 122, 87, 51, 103, 3, 46, 74, 116, 16, 61, 89, 65, 37, 8, 108, 82, 54, 27, 127},
{0, 103, 79, 40, 31, 120, 80, 55, 62, 89, 113, 22, 33, 70, 110, 9, 124, 27, 51, 84, 99, 4, 44, 75, 66, 37, 13, 106, 93, 58, 18, 117, 121, 30, 54, 81, 102, 1, 41, 78, 71, 32, 8, 111, 88, 63, 23, 112, 5, 98, 74, 45, 26, 125, 85, 50, 59, 92, 116, 19, 36, 67, 107, 12, 115, 20, 60, 91, 108, 11, 35, 68, 77, 42, 2, 101, 82, 53, 29, 122, 15, 104, 64, 39, 16, 119, 95, 56, 49, 86, 126, 25, 46, 73, 97, 6, 10, 109, 69, 34, 21, 114, 90, 61, 52, 83, 123, 28, 43, 76, 100, 3, 118, 17, 57, 94, 105, 14, 38, 65, 72, 47, 7, 96, 87, 48, 24, 127},
{0, 38, 76, 106, 25, 63, 85, 115, 50, 20, 126, 88, 43, 13, 103, 65, 100, 66, 40, 14, 125, 91, 49, 23, 86, 112, 26, 60, 79, 105, 3, 37, 73, 111, 5, 35, 80, 118, 28, 58, 123, 93, 55, 17, 98, 68, 46, 8, 45, 11, 97, 71, 52, 18, 120, 94, 31, 57, 83, 117, 6, 32, 74, 108, 19, 53, 95, 121, 10, 44, 70, 96, 33, 7, 109, 75, 56, 30, 116, 82, 119, 81, 59, 29, 110, 72, 34, 4, 69, 99, 9, 47, 92, 122, 16, 54, 90, 124, 22, 48, 67, 101, 15, 41, 104, 78, 36, 2, 113, 87, 61, 27, 62, 24, 114, 84, 39, 1, 107, 77, 12, 42, 64, 102, 21, 51, 89, 127},
{0, 37, 74, 111, 21, 48, 95, 122, 42, 15, 96, 69, 63, 26, 117, 80, 84, 113, 30, 59, 65, 100, 11, 46, 126, 91, 52, 17, 107, 78, 33, 4, 41, 12, 99, 70, 60, 25, 118, 83, 3, 38, 73, 108, 22, 51, 92, 121, 125, 88, 55, 18, 104, 77, 34, 7, 87, 114, 29, 56, 66, 103, 8, 45, 82, 119, 24, 61, 71, 98, 13, 40, 120, 93, 50, 23, 109, 72, 39, 2, 6, 35, 76, 105, 19, 54, 89, 124, 44, 9, 102, 67, 57, 28, 115, 86, 123, 94, 49, 20, 110, 75, 36, 1, 81, 116, 27, 62, 68, 97, 14, 43, 47, 10, 101, 64, 58, 31, 112, 85, 5, 32, 79, 106, 16, 53, 90, 127},
{0, 84, 41, 125, 82, 6, 123, 47, 37, 113, 12, 88, 119, 35, 94, 10, 74, 30, 99, 55, 24, 76, 49, 101, 111, 59, 70, 18, 61, 105, 20, 64, 21, 65, 60, 104, 71, 19, 110, 58, 48, 100, 25, 77, 98, 54, 75, 31, 95, 11, 118, 34, 13, 89, 36, 112, 122, 46, 83, 7, 40, 124, 1, 85, 42, 126, 3, 87, 120, 44, 81, 5, 15, 91, 38, 114, 93, 9, 116, 32, 96, 52, 73, 29, 50, 102, 27, 79, 69, 17, 108, 56, 23, 67, 62, 106, 63, 107, 22, 66, 109, 57, 68, 16, 26, 78, 51, 103, 72, 28, 97, 53, 117, 33, 92, 8, 39, 115, 14, 90, 80, 4, 121, 45, 2, 86, 43, 127},
{0, 87, 47, 120, 94, 9, 113, 38, 61, 106, 18, 69, 99, 52, 76, 27, 122, 45, 85, 2, 36, 115, 11, 92, 71, 16, 104, 63, 25, 78, 54, 97, 117, 34, 90, 13, 43, 124, 4, 83, 72, 31, 103, 48, 22, 65, 57, 110, 15, 88, 32, 119, 81, 6, 126, 41, 50, 101, 29, 74, 108, 59, 67, 20, 107, 60, 68, 19, 53, 98, 26, 77, 86, 1, 121, 46, 8, 95, 39, 112, 17, 70, 62, 105, 79, 24, 96, 55, 44, 123, 3, 84, 114, 37, 93, 10, 30, 73, 49, 102, 64, 23, 111, 56, 35, 116, 12, 91, 125, 42, 82, 5, 100, 51, 75, 28, 58, 109, 21, 66, 89, 14, 118, 33, 7, 80, 40, 127},
{0, 21, 42, 63, 84, 65, 126, 107, 41, 60, 3, 22, 125, 104, 87, 66, 82, 71, 120, 109, 6, 19, 44, 57, 123, 110, 81, 68, 47, 58, 5, 16, 37, 48, 15, 26, 113, 100, 91, 78, 12, 25, 38, 51, 88, 77, 114, 103, 119, 98, 93, 72, 35, 54, 9, 28, 94, 75, 116, 97, 10, 31, 32, 53, 74, 95, 96, 117, 30, 11, 52, 33, 99, 118, 73, 92, 55, 34, 29, 8, 24, 13, 50, 39, 76, 89, 102, 115, 49, 36, 27, 14, 101, 112, 79, 90, 111, 122, 69, 80, 59, 46, 17, 4, 70, 83, 108, 121, 18, 7, 56, 45, 61, 40, 23, 2, 105, 124, 67, 86, 20, 1, 62, 43, 64, 85, 106, 127},
{0, 118, 109, 27, 91, 45, 54, 64, 55, 65, 90, 44, 108, 26, 1, 119, 110, 24, 3, 117, 53, 67, 88, 46, 89, 47, 52, 66, 2, 116, 111, 25, 93, 43, 48, 70, 6, 112, 107, 29, 106, 28, 7, 113, 49, 71, 92, 42, 51, 69, 94, 40, 104, 30, 5, 115, 4, 114, 105, 31, 95, 41, 50, 68, 59, 77, 86, 32, 96, 22, 13, 123, 12, 122, 97, 23, 87, 33, 58, 76, 85, 35, 56, 78, 14, 120, 99, 21, 98, 20, 15, 121, 57, 79, 84, 34, 102, 16, 11, 125, 61, 75, 80, 38, 81, 39, 60, 74, 10, 124, 103, 17, 8, 126, 101, 19, 83, 37, 62, 72, 63, 73, 82, 36, 100, 18, 9, 127},
{0, 117, 107, 30, 87, 34, 60, 73, 47, 90, 68, 49, 120, 13, 19, 102, 94, 43, 53, 64, 9, 124, 98, 23, 113, 4, 26, 111, 38, 83, 77, 56, 61, 72, 86, 35, 106, 31, 1, 116, 18, 103, 121, 12, 69, 48, 46, 91, 99, 22, 8, 125, 52, 65, 95, 42, 76, 57, 39, 82, 27, 110, 112, 5, 122, 15, 17, 100, 45, 88, 70, 51, 85, 32, 62, 75, 2, 119, 105, 28, 36, 81, 79, 58, 115, 6, 24, 109, 11, 126, 96, 21, 92, 41, 55, 66, 71, 50, 44, 89, 16, 101, 123, 14, 104, 29, 3, 118, 63, 74, 84, 33, 25, 108, 114, 7, 78, 59, 37, 80, 54, 67, 93, 40, 97, 20, 10, 127},
{0, 55, 110, 89, 93, 106, 51, 4, 59, 12, 85, 98, 102, 81, 8, 63, 118, 65, 24, 47, 43, 28, 69, 114, 77, 122, 35, 20, 16, 39, 126, 73, 109, 90, 3, 52, 48, 7, 94, 105, 86, 97, 56, 15, 11, 60, 101, 82, 27, 44, 117, 66, 70, 113, 40, 31, 32, 23, 78, 121, 125, 74, 19, 36, 91, 108, 53, 2, 6, 49, 104, 95, 96, 87, 14, 57, 61, 10, 83, 100, 45, 26, 67, 116, 112, 71, 30, 41, 22, 33, 120, 79, 75, 124, 37, 18, 54, 1, 88, 111, 107, 92, 5, 50, 13, 58, 99, 84, 80, 103, 62, 9, 64, 119, 46, 25, 29, 42, 115, 68, 123, 76, 21, 34, 38, 17, 72, 127},
{0, 76, 25, 85, 50, 126, 43, 103, 100, 40, 125, 49, 86, 26, 79, 3, 73, 5, 80, 28, 123, 55, 98, 46, 45, 97, 52, 120, 31, 83, 6, 74, 19, 95, 10, 70, 33, 109, 56, 116, 119, 59, 110, 34, 69, 9, 92, 16, 90, 22, 67, 15, 104, 36, 113, 61, 62, 114, 39, 107, 12, 64, 21, 89, 38, 106, 63, 115, 20, 88, 13, 65, 66, 14, 91, 23, 112, 60, 105, 37, 111, 35, 118, 58, 93, 17, 68, 8, 11, 71, 18, 94, 57, 117, 32, 108, 53, 121, 44, 96, 7, 75, 30, 82, 81, 29, 72, 4, 99, 47, 122, 54, 124, 48, 101, 41, 78, 2, 87, 27, 24, 84, 1, 77, 42, 102, 51, 127},
{0, 79, 31, 80, 62, 113, 33, 110, 124, 51, 99, 44, 66, 13, 93, 18, 121, 54, 102, 41, 71, 8, 88, 23, 5, 74, 26, 85, 59, 116, 36, 107, 115, 60, 108, 35, 77, 2, 82, 29, 15, 64, 16, 95, 49, 126, 46, 97, 10, 69, 21, 90, 52, 123, 43, 100, 118, 57, 105, 38, 72, 7, 87, 24, 103, 40, 120, 55, 89, 22, 70, 9, 27, 84, 4, 75, 37, 106, 58, 117, 30, 81, 1, 78, 32, 111, 63, 112, 98, 45, 125, 50, 92, 19, 67, 12, 20, 91, 11, 68, 42, 101, 53, 122, 104, 39, 119, 56, 86, 25, 73, 6, 109, 34, 114, 61, 83, 28, 76, 3, 17, 94, 14, 65, 47, 96, 48, 127},
{0, 14, 28, 18, 56, 54, 36, 42, 112, 126, 108, 98, 72, 70, 84, 90, 97, 111, 125, 115, 89, 87, 69, 75, 17, 31, 13, 3, 41, 39, 53, 59, 67, 77, 95, 81, 123, 117, 103, 105, 51, 61, 47, 33, 11, 5, 23, 25, 34, 44, 62, 48, 26, 20, 6, 8, 82, 92, 78, 64, 106, 100, 118, 120, 7, 9, 27, 21, 63, 49, 35, 45, 119, 121, 107, 101, 79, 65, 83, 93, 102, 104, 122, 116, 94, 80, 66, 76, 22, 24, 10, 4, 46, 32, 50, 60, 68, 74, 88, 86, 124, 114, 96, 110, 52, 58, 40, 38, 12, 2, 16, 30, 37, 43, 57, 55, 29, 19, 1, 15, 85, 91, 73, 71, 109, 99, 113, 127},
{0, 110, 93, 51, 59, 85, 102, 8, 118, 24, 43, 69, 77, 35, 16, 126, 109, 3, 48, 94, 86, 56, 11, 101, 27, 117, 70, 40, 32, 78, 125, 19, 91, 53, 6, 104, 96, 14, 61, 83, 45, 67, 112, 30, 22, 120, 75, 37, 54, 88, 107, 5, 13, 99, 80, 62, 64, 46, 29, 115, 123, 21, 38, 72, 55, 89, 106, 4, 12, 98, 81, 63, 65, 47, 28, 114, 122, 20, 39, 73, 90, 52, 7, 105, 97, 15, 60, 82, 44, 66, 113, 31, 23, 121, 74, 36, 108, 2, 49, 95, 87, 57, 10, 100, 26, 116, 71, 41, 33, 79, 124, 18, 1, 111, 92, 50, 58, 84, 103, 9, 119, 25, 42, 68, 76, 34, 17, 127},
{0, 109, 91, 54, 55, 90, 108, 1, 110, 3, 53, 88, 89, 52, 2, 111, 93, 48, 6, 107, 106, 7, 49, 92, 51, 94, 104, 5, 4, 105, 95, 50, 59, 86, 96, 13, 12, 97, 87, 58, 85, 56, 14, 99, 98, 15, 57, 84, 102, 11, 61, 80, 81, 60, 10, 103, 8, 101, 83, 62, 63, 82, 100, 9, 118, 27, 45, 64, 65, 44, 26, 119, 24, 117, 67, 46, 47, 66, 116, 25, 43, 70, 112, 29, 28, 113, 71, 42, 69, 40, 30, 115, 114, 31, 41, 68, 77, 32, 22, 123, 122, 23, 33, 76, 35, 78, 120, 21, 20, 121, 79, 34, 16, 125, 75, 38, 39, 74, 124, 17, 126, 19, 37, 72, 73, 36, 18, 127},
{0, 47, 94, 113, 61, 18, 99, 76, 122, 85, 36, 11, 71, 104, 25, 54, 117, 90, 43, 4, 72, 103, 22, 57, 15, 32, 81, 126, 50, 29, 108, 67, 107, 68, 53, 26, 86, 121, 8, 39, 17, 62, 79, 96, 44, 3, 114, 93, 30, 49, 64, 111, 35, 12, 125, 82, 100, 75, 58, 21, 89, 118, 7, 40, 87, 120, 9, 38, 106, 69, 52, 27, 45, 2, 115, 92, 16, 63, 78, 97, 34, 13, 124, 83, 31, 48, 65, 110, 88, 119, 6, 41, 101, 74, 59, 20, 60, 19, 98, 77, 1, 46, 95, 112, 70, 105, 24, 55, 123, 84, 37, 10, 73, 102, 23, 56, 116, 91, 42, 5, 51, 28, 109, 66, 14, 33, 80, 127},
{0, 94, 61, 99, 122, 36, 71, 25, 117, 43, 72, 22, 15, 81, 50, 108, 107, 53, 86, 8, 17, 79, 44, 114, 30, 64, 35, 125, 100, 58, 89, 7, 87, 9, 106, 52, 45, 115, 16, 78, 34, 124, 31, 65, 88, 6, 101, 59, 60, 98, 1, 95, 70, 24, 123, 37, 73, 23, 116, 42, 51, 109, 14, 80, 47, 113, 18, 76, 85, 11, 104, 54, 90, 4, 103, 57, 32, 126, 29, 67, 68, 26, 121, 39, 62, 96, 3, 93, 49, 111, 12, 82, 75, 21, 118, 40, 120, 38, 69, 27, 2, 92, 63, 97, 13, 83, 48, 110, 119, 41, 74, 20, 19, 77, 46, 112, 105, 55, 84, 10, 102, 56, 91, 5, 28, 66, 33, 127},
{0, 93, 59, 102, 118, 43, 77, 16, 109, 48, 86, 11, 27, 70, 32, 125, 91, 6, 96, 61, 45, 112, 22, 75, 54, 107, 13, 80, 64, 29, 123, 38, 55, 106, 12, 81, 65, 28, 122, 39, 90, 7, 97, 60, 44, 113, 23, 74, 108, 49, 87, 10, 26, 71, 33, 124, 1, 92, 58, 103, 119, 42, 76, 17, 110, 51, 85, 8, 24, 69, 35, 126, 3, 94, 56, 101, 117, 40, 78, 19, 53, 104, 14, 83, 67, 30, 120, 37, 88, 5, 99, 62, 46, 115, 21, 72, 89, 4, 98, 63, 47, 114, 20, 73, 52, 105, 15, 82, 66, 31, 121, 36, 2, 95, 57, 100, 116, 41, 79, 18, 111, 50, 84, 9, 25, 68, 34, 127},
{0, 28, 56, 36, 112, 108, 72, 84, 97, 125, 89, 69, 17, 13, 41, 53, 67, 95, 123, 103, 51, 47, 11, 23, 34, 62, 26, 6, 82, 78, 106, 118, 7, 27, 63, 35, 119, 107, 79, 83, 102, 122, 94, 66, 22, 10, 46, 50, 68, 88, 124, 96, 52, 40, 12, 16, 37, 57, 29, 1, 85, 73, 109, 113, 14, 18, 54, 42, 126, 98, 70, 90, 111, 115, 87, 75, 31, 3, 39, 59, 77, 81, 117, 105, 61, 33, 5, 25, 44, 48, 20, 8, 92, 64, 100, 120, 9, 21, 49, 45, 121, 101, 65, 93, 104, 116, 80, 76, 24, 4, 32, 60, 74, 86, 114, 110, 58, 38, 2, 30, 43, 55, 19, 15, 91, 71, 99, 127},
{0, 31, 62, 33, 124, 99, 66, 93, 121, 102, 71, 88, 5, 26, 59, 36, 115, 108, 77, 82, 15, 16, 49, 46, 10, 21, 52, 43, 118, 105, 72, 87, 103, 120, 89, 70, 27, 4, 37, 58, 30, 1, 32, 63, 98, 125, 92, 67, 20, 11, 42, 53, 104, 119, 86, 73, 109, 114, 83, 76, 17, 14, 47, 48, 79, 80, 113, 110, 51, 44, 13, 18, 54, 41, 8, 23, 74, 85, 116, 107, 60, 35, 2, 29, 64, 95, 126, 97, 69, 90, 123, 100, 57, 38, 7, 24, 40, 55, 22, 9, 84, 75, 106, 117, 81, 78, 111, 112, 45, 50, 19, 12, 91, 68, 101, 122, 39, 56, 25, 6, 34, 61, 28, 3, 94, 65, 96, 127},
{0, 124, 121, 5, 115, 15, 10, 118, 103, 27, 30, 98, 20, 104, 109, 17, 79, 51, 54, 74, 60, 64, 69, 57, 40, 84, 81, 45, 91, 39, 34, 94, 31, 99, 102, 26, 108, 16, 21, 105, 120, 4, 1, 125, 11, 119, 114, 14, 80, 44, 41, 85, 35, 95, 90, 38, 55, 75, 78, 50, 68, 56, 61, 65, 62, 66, 71, 59, 77, 49, 52, 72, 89, 37, 32, 92, 42, 86, 83, 47, 113, 13, 8, 116, 2, 126, 123, 7, 22, 106, 111, 19, 101, 25, 28, 96, 33, 93, 88, 36, 82, 46, 43, 87, 70, 58, 63, 67, 53, 73, 76, 48, 110, 18, 23, 107, 29, 97, 100, 24, 9, 117, 112, 12, 122, 6, 3, 127},
{0, 62, 124, 66, 121, 71, 5, 59, 115, 77, 15, 49, 10, 52, 118, 72, 103, 89, 27, 37, 30, 32, 98, 92, 20, 42, 104, 86, 109, 83, 17, 47, 79, 113, 51, 13, 54, 8, 74, 116, 60, 2, 64, 126, 69, 123, 57, 7, 40, 22, 84, 106, 81, 111, 45, 19, 91, 101, 39, 25, 34, 28, 94, 96, 31, 33, 99, 93, 102, 88, 26, 36, 108, 82, 16, 46, 21, 43, 105, 87, 120, 70, 4, 58, 1, 63, 125, 67, 11, 53, 119, 73, 114, 76, 14, 48, 80, 110, 44, 18, 41, 23, 85, 107, 35, 29, 95, 97, 90, 100, 38, 24, 55, 9, 75, 117, 78, 112, 50, 12, 68, 122, 56, 6, 61, 3, 65, 127},
{0, 61, 122, 71, 117, 72, 15, 50, 107, 86, 17, 44, 30, 35, 100, 89, 87, 106, 45, 16, 34, 31, 88, 101, 60, 1, 70, 123, 73, 116, 51, 14, 47, 18, 85, 104, 90, 103, 32, 29, 68, 121, 62, 3, 49, 12, 75, 118, 120, 69, 2, 63, 13, 48, 119, 74, 19, 46, 105, 84, 102, 91, 28, 33, 94, 99, 36, 25, 43, 22, 81, 108, 53, 8, 79, 114, 64, 125, 58, 7, 9, 52, 115, 78, 124, 65, 6, 59, 98, 95, 24, 37, 23, 42, 109, 80, 113, 76, 11, 54, 4, 57, 126, 67, 26, 39, 96, 93, 111, 82, 21, 40, 38, 27, 92, 97, 83, 110, 41, 20, 77, 112, 55, 10, 56, 5, 66, 127}
};
