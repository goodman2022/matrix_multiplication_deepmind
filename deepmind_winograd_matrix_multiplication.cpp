
// Implemented by Sayantan Ghorai
#include <bits/stdc++.h>
using namespace std;
static void naive_multiplication(float *A, float *B, float *C, const int M, const int K, const int N)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            float tot = 0.0f;
            for (int k = 0; k < K; k++)
            {
                tot += A[i * K + k] * B[k * N + j];
            }
            C[i * N + j] = tot;
        }
    }
}
static void Brutemultiplication(float *matA, float *matB, float *matC, const int M, const int N, const int K,
                                const int strideA, const int strideB, const int strideC)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            float tot = 0.0f;
            for (int k = 0; k < K; k++)
            {
                tot += matA[i * strideA + k] * matB[k * strideB + j];
            }
            matC[i * strideC + j] = tot;
        }
    }
}

static void Winograd(float *matA, float *matB, float *matC, const int M, const int N, const int K,
                     const int blockA, const int blockB, const int blockC)
{
    if ((M <= 2) || (M % 2 != 0 || N % 2 != 0 || K % 2 != 0))
    {
        return Brutemultiplication(matA, matB, matC, M, N, K, blockA, blockB, blockC);
    }

    float *S1 = (float *)malloc((M / 2) * (K / 2) * sizeof(float));
    float *S2 = (float *)malloc((M / 2) * (K / 2) * sizeof(float));
    float *S3 = (float *)malloc((M / 2) * (K / 2) * sizeof(float));
    float *S4 = (float *)malloc((M / 2) * (K / 2) * sizeof(float));
    for (int i = 0; i < M / 2; i++)
    {
        for (int j = 0; j < K / 2; j++)
        {
            int idxA, to_be_advanced, idxS = i * (K / 2) + j;

            // S1     = A21 + A22
            idxA = (i + (M / 2)) * blockA + j;
            to_be_advanced = K / 2;
            S1[idxS] = matA[idxA] + matA[idxA + to_be_advanced];

            // S2     = S1 - A11
            idxA = i * blockA + j;
            S2[idxS] = S1[idxS] - matA[idxA];

            // S3     = A11 - A21
            to_be_advanced = (M / 2) * blockA;
            S3[idxS] = matA[idxA] - matA[idxA + to_be_advanced];

            // S4     = A12 - S2
            idxA = i * blockA + (K / 2) + j;
            S4[idxS] = matA[idxA] - S2[idxS];
        }
    }

    float *T1 = (float *)malloc((K / 2) * (N / 2) * sizeof(float));
    float *T2 = (float *)malloc((K / 2) * (N / 2) * sizeof(float));
    float *T3 = (float *)malloc((K / 2) * (N / 2) * sizeof(float));
    float *T4 = (float *)malloc((K / 2) * (N / 2) * sizeof(float));

    for (int i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < N / 2; j++)
        {
            int idxB, to_be_advanced, idxT = i * (N / 2) + j;

            // T1     = B12 - B11
            idxB = i * blockB + j;
            to_be_advanced = (N / 2);
            T1[idxT] = matB[idxB + to_be_advanced] - matB[idxB];

            // T2     = B22 - T1
            idxB = (i + (K / 2)) * blockB + (N / 2) + j;
            T2[idxT] = matB[idxB] - T1[idxT];

            // T3     = B22 - B12
            idxB = i * blockB + (N / 2) + j;
            to_be_advanced = ((K / 2)) * blockB;
            T3[idxT] = matB[idxB + to_be_advanced] - matB[idxB];

            // T4     = T2 - B21
            idxB = (i + (K / 2)) * blockB + j;
            T4[idxT] = T2[idxT] - matB[idxB];
        }
    }

    // M1 = A11 * B11
    float *M1 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(matA, matB, &M1[0], M / 2, N / 2, K / 2, blockA, blockB, N / 2);

    // M2 = A12 * B21
    float *M2 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(&matA[K / 2], &matB[(K / 2) * blockB], &M2[0], M / 2, N / 2, K / 2, blockA, blockB, N / 2);

    // M3 = S4 * B22
    float *M3 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(&S4[0], &matB[(K / 2) * blockB + (N / 2)], &M3[0], M / 2, N / 2, K / 2, K / 2, blockB, N / 2);

    // M4 = A22 * T4
    float *M4 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(&matA[(M / 2) * blockA + (K / 2)], &T4[0], &M4[0], M / 2, N / 2, K / 2, blockA, N / 2, N / 2);

    // M5 = S1 * T1
    float *M5 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(&S1[0], &T1[0], &M5[0], M / 2, N / 2, K / 2, K / 2, N / 2, N / 2);

    // M6 = S2 * T2
    float *M6 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(&S2[0], &T2[0], &M6[0], M / 2, N / 2, K / 2, K / 2, N / 2, N / 2);

    // M7 = S3 * T3
    float *M7 = (float *)malloc((M / 2) * (N / 2) * sizeof(float));
    Winograd(&S3[0], &T3[0], &M7[0], M / 2, N / 2, K / 2, K / 2, N / 2, N / 2);

    // C11 = U1 = M1 + M2
    // C12 = U5 = U4 + M3 = U2 + M5 + M3 = M1 + M6 + M5 + M3
    // C21 = U6 = U3 - M4 = U2 + M7 - M4 = M1 + M6 + M7 - M4
    // C22 = U7 = U3 + M5 = U2 + M7 + M5 = M1 + M6 + M7 + M5
    for (int i = 0; i < M / 2; i++)
    {
        for (int j = 0; j < N / 2; j++)
        {
            int idx = i * (N / 2) + j;
            matC[i * blockC + j] = M1[idx] + M2[idx];
            matC[i * blockC + j + (N / 2)] = M1[idx] + M6[idx] + M5[idx] + M3[idx];
            matC[(i + (M / 2)) * blockC + j] = M1[idx] + M6[idx] + M7[idx] - M4[idx];
            matC[(i + (M / 2)) * blockC + j + (N / 2)] = M1[idx] + M6[idx] + M7[idx] + M5[idx];
        }
    }
}

static void Winograd(float *matA, float *matB, float *matC, const int M, const int N, const int K)
{
    Winograd(matA, matB, matC, M, N, K, K, N, N);
}
void deepmind(float *matA, float *matB, float *matC, const int M, const int K, const int N)
{
    if ((M <= 4) || (K <= 5) || (N <= 5) || M % 4 != 0 || N % 5 != 0 || K % 5 != 0)
    {
        return Winograd(matA, matB, matC, M, N, K);
    }
    float *A11 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A12 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A13 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A14 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A15 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A21 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A22 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A23 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A24 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A25 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A31 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A32 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A33 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A34 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A35 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A41 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A42 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A43 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A44 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *A45 = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    int dx = 0, dy = 0;
    for (int i = 0; i < M / 4; i++)
    {

        for (int j = 0; j < K / 5; j++)
        {
            dx = 0, dy = 0;

            A11[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A12[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A13[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A14[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A15[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];

            dx += (M / 4);
            dy = 0;

            A21[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A22[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A23[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A24[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A25[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];

            dx += (M / 4);
            dy = 0;

            A31[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A32[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A33[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A34[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A35[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];

            dx += (M / 4);
            dy = 0;

            A41[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A42[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A43[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A44[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
            dy += (K / 5);
            A45[i * (K / 5) + j] = matA[(i + dx) * K + j + dy];
        }
    }
    float *B11 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B12 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B13 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B14 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B15 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B21 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B22 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B23 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B24 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B25 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B31 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B32 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B33 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B34 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B35 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B41 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B42 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B43 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B44 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B45 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B51 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B52 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B53 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B54 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    float *B55 = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    dx = 0, dy = 0;
    for (int i = 0; i < K / 5; i++)
    {
        for (int j = 0; j < N / 5; j++)
        {
            dx = 0, dy = 0;

            B11[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B12[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B13[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B14[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B15[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];

            dx += (K / 5);
            dy = 0;

            B21[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B22[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B23[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B24[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B25[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];

            dx += (K / 5);
            dy = 0;

            B31[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B32[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B33[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B34[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B35[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];

            dx += (K / 5);
            dy = 0;

            B41[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B42[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B43[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B44[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B45[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];

            dx += (K / 5);
            dy = 0;

            B51[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B52[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B53[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B54[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
            dy += (N / 5);
            B55[i * (N / 5) + j] = matB[(i + dx) * N + j + dy];
        }
    }
    float *tempA = (float *)malloc((M / 4) * (K / 5) * sizeof(float));
    float *tempB = (float *)malloc((K / 5) * (N / 5) * sizeof(float));
    // h1=A32*(-B21-B25-B31)
    float *h1 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A32[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B21[i * (N / 5) + j] - B25[i * (N / 5) + j] - B31[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h1[0], M / 4, K / 5, N / 5);
    }
    // h2=(A22+A25-A35)*(-B25-B51)
    float *h2 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A22[i * (K / 5) + j] + A25[i * (K / 5) + j] - A35[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B25[i * (N / 5) + j] - B51[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h2[0], M / 4, K / 5, N / 5);
    }
    // h3=(-A31-A41+A42)*(-B11+B25)
    float *h3 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A31[i * (K / 5) + j] - A41[i * (K / 5) + j] + A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B11[i * (N / 5) + j] + B25[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h3[0], M / 4, K / 5, N / 5);
    }
    // h4=(A12+A14+A34)*(-B25-B41)
    float *h4 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A12[i * (K / 5) + j] + A14[i * (K / 5) + j] + A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B25[i * (N / 5) + j] - B41[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h4[0], M / 4, K / 5, N / 5);
    }
    // h5=(A15+A22+A25)*(-B24+B51)
    float *h5 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A15[i * (K / 5) + j] + A22[i * (K / 5) + j] + A25[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B24[i * (N / 5) + j] + B51[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h5[0], M / 4, K / 5, N / 5);
    }
    // h6=(-A22-A25-A45)*(B23+B51)
    float *h6 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A22[i * (K / 5) + j] - A25[i * (K / 5) + j] - A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B23[i * (N / 5) + j] + B51[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h6[0], M / 4, K / 5, N / 5);
    }
    // h7=(-A11+A41-A42)*(B11+B24)
    float *h7 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A11[i * (K / 5) + j] + A41[i * (K / 5) + j] - A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B24[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h7[0], M / 4, K / 5, N / 5);
    }
    // h8=(A32-A33-A43)*(-B23+B31)
    float *h8 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A32[i * (K / 5) + j] - A33[i * (K / 5) + j] - A43[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B23[i * (N / 5) + j] + B31[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h8[0], M / 4, K / 5, N / 5);
    }
    // h9=(-A12-A14+A44)*(B23+B41)
    float *h9 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A12[i * (K / 5) + j] - A14[i * (K / 5) + j] + A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B23[i * (N / 5) + j] + B41[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h9[0], M / 4, K / 5, N / 5);
    }
    // h10=(A22+A25)*(B51)
    float *h10 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A22[i * (K / 5) + j] + A25[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B51[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h10[0], M / 4, K / 5, N / 5);
    }
    // h11=(-A21-A41+A42)*(-B11+B22)
    float *h11 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A21[i * (K / 5) + j] - A41[i * (K / 5) + j] + A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B11[i * (N / 5) + j] + B22[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h11[0], M / 4, K / 5, N / 5);
    }
    // h12=(A41-A42)*(B11)
    float *h12 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A41[i * (K / 5) + j] - A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h12[0], M / 4, K / 5, N / 5);
    }
    // h13=(A12+A14+A24)*(B22+B41)
    float *h13 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A12[i * (K / 5) + j] + A14[i * (K / 5) + j] + A24[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B22[i * (N / 5) + j] + B41[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h13[0], M / 4, K / 5, N / 5);
    }
    // h14=(A13-A32+A33)*(B24+B31)
    float *h14 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A13[i * (K / 5) + j] - A32[i * (K / 5) + j] + A33[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B24[i * (N / 5) + j] + B31[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h14[0], M / 4, K / 5, N / 5);
    }
    // h15=(-A12-A14)*(B41)
    float *h15 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A12[i * (K / 5) + j] - A14[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B41[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h15[0], M / 4, K / 5, N / 5);
    }
    // h16=(-A32+A33)*(B31)
    float *h16 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A32[i * (K / 5) + j] + A33[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B31[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h16[0], M / 4, K / 5, N / 5);
    }
    // h17=(A12+A14-A21+A22-A23+A24-A32+A33-A41+A42)*(B22)
    float *h17 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A12[i * (K / 5) + j] + A14[i * (K / 5) + j] - A21[i * (K / 5) + j] + A22[i * (K / 5) + j] - A23[i * (K / 5) + j] + A24[i * (K / 5) + j] - A32[i * (K / 5) + j] + A33[i * (K / 5) + j] - A41[i * (K / 5) + j] + A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B22[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h17[0], M / 4, K / 5, N / 5);
    }
    // h18=(A21)*(B11+B12+B52)
    float *h18 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A21[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B12[i * (N / 5) + j] + B52[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h18[0], M / 4, K / 5, N / 5);
    }
    // h19=(-A23)*(B31+B32+B52)
    float *h19 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A23[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B31[i * (N / 5) + j] + B32[i * (N / 5) + j] + B52[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h19[0], M / 4, K / 5, N / 5);
    }
    // h20=(-A15+A21+A23-A25)*(-B11-B12+B14-B52)
    float *h20 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A15[i * (K / 5) + j] + A21[i * (K / 5) + j] + A23[i * (K / 5) + j] - A25[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B11[i * (N / 5) + j] - B12[i * (N / 5) + j] + B14[i * (N / 5) + j] - B52[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h20[0], M / 4, K / 5, N / 5);
    }
    // h21=(A21+A23-A25)*(B52)
    float *h21 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A21[i * (K / 5) + j] + A23[i * (K / 5) + j] - A25[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B52[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h21[0], M / 4, K / 5, N / 5);
    }
    // h22=(A13-A14-A24)*(B11+B12-B14-B31-B32+B34+B44)
    float *h22 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A13[i * (K / 5) + j] - A14[i * (K / 5) + j] - A24[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B12[i * (N / 5) + j] - B14[i * (N / 5) + j] - B31[i * (N / 5) + j] - B32[i * (N / 5) + j] + B34[i * (N / 5) + j] + B44[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h22[0], M / 4, K / 5, N / 5);
    }
    // h23=(A13)*(-B31+B34+B44)
    float *h23 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A13[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B31[i * (N / 5) + j] + B34[i * (N / 5) + j] + B44[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h23[0], M / 4, K / 5, N / 5);
    }
    // h24=(A15)*(-B44-B51+B54)
    float *h24 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A15[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B44[i * (N / 5) + j] - B51[i * (N / 5) + j] + B54[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h24[0], M / 4, K / 5, N / 5);
    }
    // h25=(-A11)*(B11-B14)
    float *h25 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A11[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] - B14[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h25[0], M / 4, K / 5, N / 5);
    }
    // h26=(-A13+A14+A15)*(B44)
    float *h26 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A13[i * (K / 5) + j] + A14[i * (K / 5) + j] + A15[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B44[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h26[0], M / 4, K / 5, N / 5);
    }
    // h27=(A13-A31+A33)*(B11-B14+B15+B35)
    float *h27 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A13[i * (K / 5) + j] - A31[i * (K / 5) + j] + A33[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] - B14[i * (N / 5) + j] + B15[i * (N / 5) + j] + B35[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h27[0], M / 4, K / 5, N / 5);
    }
    // h28=(-A34)*(-B35-B41-B45)
    float *h28 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B35[i * (N / 5) + j] - B41[i * (N / 5) + j] - B45[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h28[0], M / 4, K / 5, N / 5);
    }
    // h29=(A31)*(B11+B15+B35)
    float *h29 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A31[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B15[i * (N / 5) + j] + B35[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h29[0], M / 4, K / 5, N / 5);
    }
    // h30=(A31-A33+A34)*(B35)
    float *h30 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A31[i * (K / 5) + j] - A33[i * (K / 5) + j] + A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B35[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h30[0], M / 4, K / 5, N / 5);
    }
    // h31=(-A14-A15-A34)*(-B44-B51+B54-B55)
    float *h31 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A14[i * (K / 5) + j] - A15[i * (K / 5) + j] - A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B44[i * (N / 5) + j] - B51[i * (N / 5) + j] + B54[i * (N / 5) + j] - B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h31[0], M / 4, K / 5, N / 5);
    }
    // h32=(A21+A41+A44)*(B13-B41-B42-B43)
    float *h32 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A21[i * (K / 5) + j] + A41[i * (K / 5) + j] + A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B13[i * (N / 5) + j] - B41[i * (N / 5) + j] - B42[i * (N / 5) + j] - B43[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h32[0], M / 4, K / 5, N / 5);
    }
    // h33=(A43)*(-B31-B33)
    float *h33 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A43[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B31[i * (N / 5) + j] - B33[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h33[0], M / 4, K / 5, N / 5);
    }
    // h34=(A44)*(-B13+B41+B43)
    float *h34 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B13[i * (N / 5) + j] + B41[i * (N / 5) + j] + B43[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h34[0], M / 4, K / 5, N / 5);
    }
    // h35=(-A45)*(B13+B51+B53)
    float *h35 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B13[i * (N / 5) + j] + B51[i * (N / 5) + j] + B53[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h35[0], M / 4, K / 5, N / 5);
    }
    // h36=(A23-A25-A45)*(B31+B32+B33+B52)
    float *h36 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A23[i * (K / 5) + j] - A25[i * (K / 5) + j] - A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B31[i * (N / 5) + j] + B32[i * (N / 5) + j] + B33[i * (N / 5) + j] + B52[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h36[0], M / 4, K / 5, N / 5);
    }
    // h37=(-A41-A44+A45)*(B13)
    float *h37 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A41[i * (K / 5) + j] - A44[i * (K / 5) + j] + A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B13[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h37[0], M / 4, K / 5, N / 5);
    }
    // h38=(-A23-A31+A33-A34)*(B35+B41+B42+B45)
    float *h38 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A23[i * (K / 5) + j] - A31[i * (K / 5) + j] + A33[i * (K / 5) + j] - A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B35[i * (N / 5) + j] + B41[i * (N / 5) + j] + B42[i * (N / 5) + j] + B45[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h38[0], M / 4, K / 5, N / 5);
    }
    // h39=(-A31-A41-A44+A45)*(B13+B51+B53+B55)
    float *h39 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A31[i * (K / 5) + j] - A41[i * (K / 5) + j] - A44[i * (K / 5) + j] + A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B13[i * (N / 5) + j] + B51[i * (N / 5) + j] + B53[i * (N / 5) + j] + B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h39[0], M / 4, K / 5, N / 5);
    }
    // h40=(-A13+A14+A15-A44)*(-B31-B33+B34+B44)
    float *h40 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A13[i * (K / 5) + j] + A14[i * (K / 5) + j] + A15[i * (K / 5) + j] - A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B31[i * (N / 5) + j] - B33[i * (N / 5) + j] + B34[i * (N / 5) + j] + B44[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h40[0], M / 4, K / 5, N / 5);
    }
    // h41=(-A11+A41-A45)*(B13+B31+B33-B34+B51+B53-B54)
    float *h41 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A11[i * (K / 5) + j] + A41[i * (K / 5) + j] - A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B13[i * (N / 5) + j] + B31[i * (N / 5) + j] + B33[i * (N / 5) + j] - B34[i * (N / 5) + j] + B51[i * (N / 5) + j] + B53[i * (N / 5) + j] - B54[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h41[0], M / 4, K / 5, N / 5);
    }
    // h42=(-A21+A25-A35)*(-B11-B12-B15+B41+B42+B45-B52)
    float *h42 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A21[i * (K / 5) + j] + A25[i * (K / 5) + j] - A35[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B11[i * (N / 5) + j] - B12[i * (N / 5) + j] - B15[i * (N / 5) + j] + B41[i * (N / 5) + j] + B42[i * (N / 5) + j] + B45[i * (N / 5) + j] - B52[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h42[0], M / 4, K / 5, N / 5);
    }
    // h43=(A24)*(B41+B42)
    float *h43 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A24[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B41[i * (N / 5) + j] + B42[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h43[0], M / 4, K / 5, N / 5);
    }
    // h44=(A23+A32-A33)*(B22-B31)
    float *h44 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A23[i * (K / 5) + j] + A32[i * (K / 5) + j] - A33[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B22[i * (N / 5) + j] - B31[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h44[0], M / 4, K / 5, N / 5);
    }
    // h45=(-A33+A34-A43)*(B35+B41+B43+B45+B51+B53+B55)
    float *h45 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A33[i * (K / 5) + j] + A34[i * (K / 5) + j] - A43[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B35[i * (N / 5) + j] + B41[i * (N / 5) + j] + B43[i * (N / 5) + j] + B45[i * (N / 5) + j] + B51[i * (N / 5) + j] + B53[i * (N / 5) + j] + B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h45[0], M / 4, K / 5, N / 5);
    }
    // h46=(-A35)*(-B51-B55)
    float *h46 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A35[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B51[i * (N / 5) + j] - B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h46[0], M / 4, K / 5, N / 5);
    }
    // h47=(A21-A25-A31+A35)*(B11+B12+B15-B41-B42-B45)
    float *h47 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A21[i * (K / 5) + j] - A25[i * (K / 5) + j] - A31[i * (K / 5) + j] + A35[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B12[i * (N / 5) + j] + B15[i * (N / 5) + j] - B41[i * (N / 5) + j] - B42[i * (N / 5) + j] - B45[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h47[0], M / 4, K / 5, N / 5);
    }
    // h48=(-A23+A33)*(B22+B32+B35+B41+B42+B45)
    float *h48 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A23[i * (K / 5) + j] + A33[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B22[i * (N / 5) + j] + B32[i * (N / 5) + j] + B35[i * (N / 5) + j] + B41[i * (N / 5) + j] + B42[i * (N / 5) + j] + B45[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h48[0], M / 4, K / 5, N / 5);
    }
    // h49=(-A11-A13+A14+A15-A21-A23+A24+A25)*(-B11-B12+B14)
    float *h49 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A11[i * (K / 5) + j] - A13[i * (K / 5) + j] + A14[i * (K / 5) + j] + A15[i * (K / 5) + j] - A21[i * (K / 5) + j] - A23[i * (K / 5) + j] + A24[i * (K / 5) + j] + A25[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B11[i * (N / 5) + j] - B12[i * (N / 5) + j] + B14[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h49[0], M / 4, K / 5, N / 5);
    }
    // h50=(-A14-A24)*(B22-B31-B32+B34-B42+B44)
    float *h50 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A14[i * (K / 5) + j] - A24[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B22[i * (N / 5) + j] - B31[i * (N / 5) + j] - B32[i * (N / 5) + j] + B34[i * (N / 5) + j] - B42[i * (N / 5) + j] + B44[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h50[0], M / 4, K / 5, N / 5);
    }
    // h51=(A22)*(B21+B22-B51)
    float *h51 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A22[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B21[i * (N / 5) + j] + B22[i * (N / 5) + j] - B51[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h51[0], M / 4, K / 5, N / 5);
    }
    // h52=(A42)*(B11+B21+B23)
    float *h52 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B21[i * (N / 5) + j] + B23[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h52[0], M / 4, K / 5, N / 5);
    }
    // h53=(-A12)*(-B21+B24+B41)
    float *h53 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A12[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B21[i * (N / 5) + j] + B24[i * (N / 5) + j] + B41[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h53[0], M / 4, K / 5, N / 5);
    }
    // h54=(A12+A14-A22-A25-A32+A33-A42+A43-A44-A45)*(B23)
    float *h54 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A12[i * (K / 5) + j] + A14[i * (K / 5) + j] - A22[i * (K / 5) + j] - A25[i * (K / 5) + j] - A32[i * (K / 5) + j] + A33[i * (K / 5) + j] - A42[i * (K / 5) + j] + A43[i * (K / 5) + j] - A44[i * (K / 5) + j] - A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B23[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h54[0], M / 4, K / 5, N / 5);
    }
    // h55=(A14-A44)*(-B23+B31+B33-B34+B43-B44)
    float *h55 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A14[i * (K / 5) + j] - A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B23[i * (N / 5) + j] + B31[i * (N / 5) + j] + B33[i * (N / 5) + j] - B34[i * (N / 5) + j] + B43[i * (N / 5) + j] - B44[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h55[0], M / 4, K / 5, N / 5);
    }
    // h56=(A11-A15-A41+A45)*(B31+B33-B34+B51+B53-B54)
    float *h56 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A11[i * (K / 5) + j] - A15[i * (K / 5) + j] - A41[i * (K / 5) + j] + A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B31[i * (N / 5) + j] + B33[i * (N / 5) + j] - B34[i * (N / 5) + j] + B51[i * (N / 5) + j] + B53[i * (N / 5) + j] - B54[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h56[0], M / 4, K / 5, N / 5);
    }
    // h57=(-A31-A41)*(-B13-B15-B25-B51-B53-B55)
    float *h57 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A31[i * (K / 5) + j] - A41[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B13[i * (N / 5) + j] - B15[i * (N / 5) + j] - B25[i * (N / 5) + j] - B51[i * (N / 5) + j] - B53[i * (N / 5) + j] - B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h57[0], M / 4, K / 5, N / 5);
    }
    // h58=(-A14-A15-A34-A35)*(-B51+B54-B55)
    float *h58 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A14[i * (K / 5) + j] - A15[i * (K / 5) + j] - A34[i * (K / 5) + j] - A35[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B51[i * (N / 5) + j] + B54[i * (N / 5) + j] - B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h58[0], M / 4, K / 5, N / 5);
    }
    // h59=(-A33+A34-A43+A44)*(B41+B43+B45+B51+B53+B55)
    float *h59 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A33[i * (K / 5) + j] + A34[i * (K / 5) + j] - A43[i * (K / 5) + j] + A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B41[i * (N / 5) + j] + B43[i * (N / 5) + j] + B45[i * (N / 5) + j] + B51[i * (N / 5) + j] + B53[i * (N / 5) + j] + B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h59[0], M / 4, K / 5, N / 5);
    }
    // h60=(A25+A45)*(B23-B31-B32-B33-B52-B53)
    float *h60 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A25[i * (K / 5) + j] + A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B23[i * (N / 5) + j] - B31[i * (N / 5) + j] - B32[i * (N / 5) + j] - B33[i * (N / 5) + j] - B52[i * (N / 5) + j] - B53[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h60[0], M / 4, K / 5, N / 5);
    }
    // h61=(A14+A34)*(B11-B14+B15-B25-B44+B45-B51+B54-B55)
    float *h61 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A14[i * (K / 5) + j] + A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] - B14[i * (N / 5) + j] + B15[i * (N / 5) + j] - B25[i * (N / 5) + j] - B44[i * (N / 5) + j] + B45[i * (N / 5) + j] - B51[i * (N / 5) + j] + B54[i * (N / 5) + j] - B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h61[0], M / 4, K / 5, N / 5);
    }
    // h62=(A21+A41)*(B12+B13+B22-B41-B42-B43)
    float *h62 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A21[i * (K / 5) + j] + A41[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B12[i * (N / 5) + j] + B13[i * (N / 5) + j] + B22[i * (N / 5) + j] - B41[i * (N / 5) + j] - B42[i * (N / 5) + j] - B43[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h62[0], M / 4, K / 5, N / 5);
    }
    // h63=(-A33-A43)*(-B23-B33-B35-B41-B43-B45)
    float *h63 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A33[i * (K / 5) + j] - A43[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B23[i * (N / 5) + j] - B33[i * (N / 5) + j] - B35[i * (N / 5) + j] - B41[i * (N / 5) + j] - B43[i * (N / 5) + j] - B45[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h63[0], M / 4, K / 5, N / 5);
    }
    // h64=(A11-A13-A14+A31-A33-A34)*(B11-B14+B15)
    float *h64 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A11[i * (K / 5) + j] - A13[i * (K / 5) + j] - A14[i * (K / 5) + j] + A31[i * (K / 5) + j] - A33[i * (K / 5) + j] - A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] - B14[i * (N / 5) + j] + B15[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h64[0], M / 4, K / 5, N / 5);
    }
    // h65=(-A11+A41)*(-B13+B14+B24-B51-B53+B54)
    float *h65 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A11[i * (K / 5) + j] + A41[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B13[i * (N / 5) + j] + B14[i * (N / 5) + j] + B24[i * (N / 5) + j] - B51[i * (N / 5) + j] - B53[i * (N / 5) + j] + B54[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h65[0], M / 4, K / 5, N / 5);
    }
    // h66=(A11-A12+A13-A15-A22-A25-A32+A33-A41+A42)*(B24)
    float *h66 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A11[i * (K / 5) + j] - A12[i * (K / 5) + j] + A13[i * (K / 5) + j] - A15[i * (K / 5) + j] - A22[i * (K / 5) + j] - A25[i * (K / 5) + j] - A32[i * (K / 5) + j] + A33[i * (K / 5) + j] - A41[i * (K / 5) + j] + A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B24[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h66[0], M / 4, K / 5, N / 5);
    }
    // h67=(A25-A35)*(B11+B12+B15-B25-B41-B42-B45+B52+B55)
    float *h67 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A25[i * (K / 5) + j] - A35[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B12[i * (N / 5) + j] + B15[i * (N / 5) + j] - B25[i * (N / 5) + j] - B41[i * (N / 5) + j] - B42[i * (N / 5) + j] - B45[i * (N / 5) + j] + B52[i * (N / 5) + j] + B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h67[0], M / 4, K / 5, N / 5);
    }
    // h68=(A11+A13-A14-A15-A41-A43+A44+A45)*(-B31-B33+B34)
    float *h68 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A11[i * (K / 5) + j] + A13[i * (K / 5) + j] - A14[i * (K / 5) + j] - A15[i * (K / 5) + j] - A41[i * (K / 5) + j] - A43[i * (K / 5) + j] + A44[i * (K / 5) + j] + A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B31[i * (N / 5) + j] - B33[i * (N / 5) + j] + B34[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h68[0], M / 4, K / 5, N / 5);
    }
    // h69=(-A13+A14-A23+A24)*(-B24-B31-B32+B34-B52+B54)
    float *h69 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A13[i * (K / 5) + j] + A14[i * (K / 5) + j] - A23[i * (K / 5) + j] + A24[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B24[i * (N / 5) + j] - B31[i * (N / 5) + j] - B32[i * (N / 5) + j] + B34[i * (N / 5) + j] - B52[i * (N / 5) + j] + B54[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h69[0], M / 4, K / 5, N / 5);
    }
    // h70=(A23-A25+A43-A45)*(-B31-B32-B33)
    float *h70 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A23[i * (K / 5) + j] - A25[i * (K / 5) + j] + A43[i * (K / 5) + j] - A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B31[i * (N / 5) + j] - B32[i * (N / 5) + j] - B33[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h70[0], M / 4, K / 5, N / 5);
    }
    // h71=(-A31+A33-A34+A35-A41+A43-A44+A45)*(-B51-B53-B55)
    float *h71 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A31[i * (K / 5) + j] + A33[i * (K / 5) + j] - A34[i * (K / 5) + j] + A35[i * (K / 5) + j] - A41[i * (K / 5) + j] + A43[i * (K / 5) + j] - A44[i * (K / 5) + j] + A45[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B51[i * (N / 5) + j] - B53[i * (N / 5) + j] - B55[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h71[0], M / 4, K / 5, N / 5);
    }
    // h72=(-A21-A24-A41-A44)*(B41+B42+B43)
    float *h72 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = -A21[i * (K / 5) + j] - A24[i * (K / 5) + j] - A41[i * (K / 5) + j] - A44[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B41[i * (N / 5) + j] + B42[i * (N / 5) + j] + B43[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h72[0], M / 4, K / 5, N / 5);
    }
    // h73=(A13-A14-A15+A23-A24-A25)*(B11+B12-B14+B24+B52-B54)
    float *h73 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A13[i * (K / 5) + j] - A14[i * (K / 5) + j] - A15[i * (K / 5) + j] + A23[i * (K / 5) + j] - A24[i * (K / 5) + j] - A25[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B11[i * (N / 5) + j] + B12[i * (N / 5) + j] - B14[i * (N / 5) + j] + B24[i * (N / 5) + j] + B52[i * (N / 5) + j] - B54[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h73[0], M / 4, K / 5, N / 5);
    }
    // h74=(A21-A23+A24-A31+A33-A34)*(B41+B42+B45)
    float *h74 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A21[i * (K / 5) + j] - A23[i * (K / 5) + j] + A24[i * (K / 5) + j] - A31[i * (K / 5) + j] + A33[i * (K / 5) + j] - A34[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = B41[i * (N / 5) + j] + B42[i * (N / 5) + j] + B45[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h74[0], M / 4, K / 5, N / 5);
    }
    // h75=-(A12+A14-A22-A25-A31+A32+A34+A35-A41+A42)*(B25)
    float *h75 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A12[i * (K / 5) + j] + A14[i * (K / 5) + j] - A22[i * (K / 5) + j] - A25[i * (K / 5) + j] - A31[i * (K / 5) + j] + A32[i * (K / 5) + j] + A34[i * (K / 5) + j] + A35[i * (K / 5) + j] - A41[i * (K / 5) + j] + A42[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B25[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h75[0], M / 4, K / 5, N / 5);
    }
    // h76=(A13+A33)*(-B11+B14-B15+B24+B34-B35)
    float *h76 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    {
        for (int i = 0; i < M / 4; i++)
        {
            for (int j = 0; j < K / 5; j++)
            {
                tempA[i * (K / 5) + j] = A13[i * (K / 5) + j] + A33[i * (K / 5) + j];
            }
        }
        for (int i = 0; i < K / 5; i++)
        {
            for (int j = 0; j < N / 5; j++)
            {
                tempB[i * (N / 5) + j] = -B11[i * (N / 5) + j] + B14[i * (N / 5) + j] - B15[i * (N / 5) + j] + B24[i * (N / 5) + j] + B34[i * (N / 5) + j] - B35[i * (N / 5) + j];
            }
        }
        deepmind(&tempA[0], &tempB[0], &h76[0], M / 4, K / 5, N / 5);
    }
    float *C11 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C12 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C13 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C14 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C15 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C21 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C22 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C23 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C24 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C25 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C31 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C32 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C33 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C34 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C35 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C41 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C42 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C43 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C44 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));
    float *C45 = (float *)malloc((M / 4) * (N / 5) * sizeof(float));

    for (int i = 0; i < M / 4; i++)
    {
        for (int j = 0; j < N / 5; j++)
        {
            // C11=-h10+h12+h14-h15-h16+h53+h5-h66-h7
            C11[i * (N / 5) + j] = -h10[i * (N / 5) + j] + h12[i * (N / 5) + j] + h14[i * (N / 5) + j] - h15[i * (N / 5) + j] - h16[i * (N / 5) + j] + h53[i * (N / 5) + j] + h5[i * (N / 5) + j] - h66[i * (N / 5) + j] - h7[i * (N / 5) + j];
            // C21=h10+h11-h12+h13+h15+h16-h17-h44+h51

            C21[i * (N / 5) + j] = h10[i * (N / 5) + j] + h11[i * (N / 5) + j] - h12[i * (N / 5) + j] + h13[i * (N / 5) + j] + h15[i * (N / 5) + j] + h16[i * (N / 5) + j] - h17[i * (N / 5) + j] - h44[i * (N / 5) + j] + h51[i * (N / 5) + j];
            // C31=h10-h12+h15+h16-h1+h2+h3-h4+h75

            C31[i * (N / 5) + j] = h10[i * (N / 5) + j] - h12[i * (N / 5) + j] + h15[i * (N / 5) + j] + h16[i * (N / 5) + j] - h1[i * (N / 5) + j] + h2[i * (N / 5) + j] + h3[i * (N / 5) + j] - h4[i * (N / 5) + j] + h75[i * (N / 5) + j];

            // C41=-h10+h12-h15-h16+h52+h54-h6-h8+h9

            C41[i * (N / 5) + j] = -h10[i * (N / 5) + j] + h12[i * (N / 5) + j] - h15[i * (N / 5) + j] - h16[i * (N / 5) + j] + h52[i * (N / 5) + j] + h54[i * (N / 5) + j] - h6[i * (N / 5) + j] - h8[i * (N / 5) + j] + h9[i * (N / 5) + j];

            // C12=h13+h15+h20+h21-h22+h23+h25-h43+h49+h50

            C12[i * (N / 5) + j] = h13[i * (N / 5) + j] + h15[i * (N / 5) + j] + h20[i * (N / 5) + j] + h21[i * (N / 5) + j] - h22[i * (N / 5) + j] + h23[i * (N / 5) + j] + h25[i * (N / 5) + j] - h43[i * (N / 5) + j] + h49[i * (N / 5) + j] + h50[i * (N / 5) + j];

            // C22=-h11+h12-h13-h15-h16+h17+h18-h19-h21+h43+h44

            C22[i * (N / 5) + j] = -h11[i * (N / 5) + j] + h12[i * (N / 5) + j] - h13[i * (N / 5) + j] - h15[i * (N / 5) + j] - h16[i * (N / 5) + j] + h17[i * (N / 5) + j] + h18[i * (N / 5) + j] - h19[i * (N / 5) + j] - h21[i * (N / 5) + j] + h43[i * (N / 5) + j] + h44[i * (N / 5) + j];

            // C32=-h16-h19-h21-h28-h29-h38+h42+h44-h47+h48

            C32[i * (N / 5) + j] = -h16[i * (N / 5) + j] - h19[i * (N / 5) + j] - h21[i * (N / 5) + j] - h28[i * (N / 5) + j] - h29[i * (N / 5) + j] - h38[i * (N / 5) + j] + h42[i * (N / 5) + j] + h44[i * (N / 5) + j] - h47[i * (N / 5) + j] + h48[i * (N / 5) + j];

            // C42=h11-h12-h18+h21-h32+h33-h34-h36+h62-h70

            C42[i * (N / 5) + j] = h11[i * (N / 5) + j] - h12[i * (N / 5) + j] - h18[i * (N / 5) + j] + h21[i * (N / 5) + j] - h32[i * (N / 5) + j] + h33[i * (N / 5) + j] - h34[i * (N / 5) + j] - h36[i * (N / 5) + j] + h62[i * (N / 5) + j] - h70[i * (N / 5) + j];

            // C13=h15+h23+h24+h34-h37+h40-h41+h55-h56-h9

            C13[i * (N / 5) + j] = h15[i * (N / 5) + j] + h23[i * (N / 5) + j] + h24[i * (N / 5) + j] + h34[i * (N / 5) + j] - h37[i * (N / 5) + j] + h40[i * (N / 5) + j] - h41[i * (N / 5) + j] + h55[i * (N / 5) + j] - h56[i * (N / 5) + j] - h9[i * (N / 5) + j];

            // C23=-h10+h19+h32+h35+h36+h37-h43-h60-h6-h72

            C23[i * (N / 5) + j] = -h10[i * (N / 5) + j] + h19[i * (N / 5) + j] + h32[i * (N / 5) + j] + h35[i * (N / 5) + j] + h36[i * (N / 5) + j] + h37[i * (N / 5) + j] - h43[i * (N / 5) + j] - h60[i * (N / 5) + j] - h6[i * (N / 5) + j] - h72[i * (N / 5) + j];

            // C33=-h16-h28+h33+h37-h39+h45-h46+h63-h71-h8

            C33[i * (N / 5) + j] = -h16[i * (N / 5) + j] - h28[i * (N / 5) + j] + h33[i * (N / 5) + j] + h37[i * (N / 5) + j] - h39[i * (N / 5) + j] + h45[i * (N / 5) + j] - h46[i * (N / 5) + j] + h63[i * (N / 5) + j] - h71[i * (N / 5) + j] - h8[i * (N / 5) + j];

            // C43=h10+h15+h16-h33+h34-h35-h37-h54+h6+h8-h9

            C43[i * (N / 5) + j] = h10[i * (N / 5) + j] + h15[i * (N / 5) + j] + h16[i * (N / 5) + j] - h33[i * (N / 5) + j] + h34[i * (N / 5) + j] - h35[i * (N / 5) + j] - h37[i * (N / 5) + j] - h54[i * (N / 5) + j] + h6[i * (N / 5) + j] + h8[i * (N / 5) + j] - h9[i * (N / 5) + j];

            // C14=-h10+h12+h14-h16+h23+h24+h25+h26+h5-h66-h7

            C14[i * (N / 5) + j] = -h10[i * (N / 5) + j] + h12[i * (N / 5) + j] + h14[i * (N / 5) + j] - h16[i * (N / 5) + j] + h23[i * (N / 5) + j] + h24[i * (N / 5) + j] + h25[i * (N / 5) + j] + h26[i * (N / 5) + j] + h5[i * (N / 5) + j] - h66[i * (N / 5) + j] - h7[i * (N / 5) + j];

            // C24=h10+h18-h19+h20-h22-h24-h26-h5-h69+h73

            C24[i * (N / 5) + j] = h10[i * (N / 5) + j] + h18[i * (N / 5) + j] - h19[i * (N / 5) + j] + h20[i * (N / 5) + j] - h22[i * (N / 5) + j] - h24[i * (N / 5) + j] - h26[i * (N / 5) + j] - h5[i * (N / 5) + j] - h69[i * (N / 5) + j] + h73[i * (N / 5) + j];

            // C34=-h14+h16-h23-h26+h27+h29+h31+h46-h58+h76

            C34[i * (N / 5) + j] = -h14[i * (N / 5) + j] + h16[i * (N / 5) + j] - h23[i * (N / 5) + j] - h26[i * (N / 5) + j] + h27[i * (N / 5) + j] + h29[i * (N / 5) + j] + h31[i * (N / 5) + j] + h46[i * (N / 5) + j] - h58[i * (N / 5) + j] + h76[i * (N / 5) + j];

            // C44=h12+h25+h26-h33-h35-h40+h41+h65-h68-h7

            C44[i * (N / 5) + j] = h12[i * (N / 5) + j] + h25[i * (N / 5) + j] + h26[i * (N / 5) + j] - h33[i * (N / 5) + j] - h35[i * (N / 5) + j] - h40[i * (N / 5) + j] + h41[i * (N / 5) + j] + h65[i * (N / 5) + j] - h68[i * (N / 5) + j] - h7[i * (N / 5) + j];

            // C15=h15+h24+h25+h27-h28+h30+h31-h4+h61+h64

            C15[i * (N / 5) + j] = h15[i * (N / 5) + j] + h24[i * (N / 5) + j] + h25[i * (N / 5) + j] + h27[i * (N / 5) + j] - h28[i * (N / 5) + j] + h30[i * (N / 5) + j] + h31[i * (N / 5) + j] - h4[i * (N / 5) + j] + h61[i * (N / 5) + j] + h64[i * (N / 5) + j];

            // C25=-h10-h18-h2-h30-h38+h42-h43+h46+h67+h74

            C25[i * (N / 5) + j] = -h10[i * (N / 5) + j] - h18[i * (N / 5) + j] - h2[i * (N / 5) + j] - h30[i * (N / 5) + j] - h38[i * (N / 5) + j] + h42[i * (N / 5) + j] - h43[i * (N / 5) + j] + h46[i * (N / 5) + j] + h67[i * (N / 5) + j] + h74[i * (N / 5) + j];

            // C35=-h10+h12-h15+h28+h29-h2-h30-h3+h46+h4-h75

            C35[i * (N / 5) + j] = -h10[i * (N / 5) + j] + h12[i * (N / 5) + j] - h15[i * (N / 5) + j] + h28[i * (N / 5) + j] + h29[i * (N / 5) + j] - h2[i * (N / 5) + j] - h30[i * (N / 5) + j] - h3[i * (N / 5) + j] + h46[i * (N / 5) + j] + h4[i * (N / 5) + j] - h75[i * (N / 5) + j];

            // C45=-h12-h29+h30-h34+h35+h39+h3-h45+h57+h59

            C45[i * (N / 5) + j] = -h12[i * (N / 5) + j] - h29[i * (N / 5) + j] + h30[i * (N / 5) + j] - h34[i * (N / 5) + j] + h35[i * (N / 5) + j] + h39[i * (N / 5) + j] + h3[i * (N / 5) + j] - h45[i * (N / 5) + j] + h57[i * (N / 5) + j] + h59[i * (N / 5) + j];
        }
    }
    dx = 0, dy = 0;
    for (int i = 0; i < M / 4; i++)
    {
        for (int j = 0; j < N / 5; j++)
        {
            dx = 0, dy = 0;
            matC[(i + dx) * N + j + dy] = C11[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C12[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C13[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C14[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C15[i * (N / 5) + j];
            dx += (M / 4);
            dy = 0;
            matC[(i + dx) * N + j + dy] = C21[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C22[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C23[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C24[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C25[i * (N / 5) + j];
            dx += (M / 4);
            dy = 0;
            matC[(i + dx) * N + j + dy] = C31[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C32[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C33[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C34[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C35[i * (N / 5) + j];
            dx += (M / 4);
            dy = 0;
            matC[(i + dx) * N + j + dy] = C41[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C42[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C43[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C44[i * (N / 5) + j];
            dy += (N / 5);
            matC[(i + dx) * N + j + dy] = C45[i * (N / 5) + j];
        }
    }
}
static void printMatrix(float *C, int M, int N)
{
    printf("\n**********************\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f ", C[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}
static void matrix_multiplication_test(int rowsA, int colsB, int commonDimension, int maxValue)
{
    unsigned seed = time(0);
    srand(seed);
    clock_t start, end;

    for (int i = 0; i < 2; i++)
    {
        float *matrixA = (float *)malloc(rowsA * commonDimension * sizeof(float));
        float *matrixB = (float *)malloc(commonDimension * colsB * sizeof(float));
        float *resultDeepMind = (float *)malloc(rowsA * colsB * sizeof(float));
        float *resultNaive = (float *)malloc(rowsA * colsB * sizeof(float));

        for (int j = 0; j < rowsA * commonDimension; j++)
        {
            matrixA[j] = rand() % maxValue;
        }
        for (int j = 0; j < commonDimension * colsB; j++)
        {
            matrixB[j] = rand() % maxValue;
        }

        start = clock();
        naive_multiplication(matrixA, matrixB, resultNaive, rowsA, commonDimension, colsB);
        end = clock();
        double executionTimeStrassen = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "Naive multiplication execution time: " << executionTimeStrassen * 1000 << "ms" << endl;

        start = clock();
        deepmind(matrixA, matrixB, resultDeepMind, rowsA, commonDimension, colsB);
        end = clock();
        double executionTimeNaive = (double)(end - start) / CLOCKS_PER_SEC;
        cout << "deepmind execution time: " << executionTimeNaive * 1000 << "ms" << endl;


        for (int j = 0; j < rowsA * colsB; j++)
        {
            if (resultDeepMind[j] != resultNaive[j])
            {
                cout << "Matrix A:" << endl;
                printMatrix(matrixA, rowsA, commonDimension);
                cout << "Matrix B:" << endl;
                printMatrix(matrixB, commonDimension, colsB);
                cout << "DeepMind result:" << endl;
                printMatrix(resultDeepMind, rowsA, colsB);
                cout << "Naive result:" << endl;
                printMatrix(resultNaive, rowsA, colsB);
                return;
            }
        }
        cout << endl;
    }
}

signed main()
{
    int M, N, K, maxval;
    // rowsA
    M = 1000;
    // colsB
    N = 3000;
    // common dimension
    K = 3000;
    // maximum value in input
    maxval = 10;
    matrix_multiplication_test(M, N, K, maxval);
    return 0;
}
