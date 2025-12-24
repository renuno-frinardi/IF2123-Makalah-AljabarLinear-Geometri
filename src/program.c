#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>

// Deklarasi variabel global untuk parameter kurva eliptik
int64_t p;
int64_t a;
int64_t b;

// Deginisi struktur titik dalam koordinat affine
typedef struct {
    int64_t x;
    int64_t y;
    int is_inf;
} AffinePoint;

// Deginisi struktur titik dalam koordinat proyeksi berbobot
typedef struct {
    int64_t X;
    int64_t Y;
    int64_t Z;
    int is_inf;
} WeightedPoint;

// Fungsi modular aritmatika penjumlahan
int64_t modAdd(int64_t x, int64_t y) {
    return (x + y) % p;
}

// Fungsi modular aritmatika pengurangan
int64_t modSub(int64_t x, int64_t y) {
    int64_t res = (x - y) % p;
    if (res < 0) res += p;
    return res;
}

// Fungsi modular aritmatika perkalian
int64_t modMul(int64_t x, int64_t y) {
    __int128 mul = (__int128)x * y;
    return (int64_t)(mul % p);
}

// Fungsi modular aritmatika pemangkatan
int64_t modSqr(int64_t x) {
    return modMul(x, x);
}

// Fungsi modular aritmatika pemangkatan dengan eksponen
int64_t modPow(int64_t base, int64_t exp) {
    int64_t res = 1;
    base = base % p;
    while (exp > 0) {
        if (exp % 2 == 1) res = modMul(res, base);
        base = modSqr(base);
        exp /= 2;
    }
    return res;
}

// Fungsi modular aritmatika invers
int64_t modInv(int64_t n) {
    return modPow(n, p - 2);
}

// Operasi perkalian double titik pada kurva eliptik
AffinePoint affineDouble(AffinePoint P) {
    if (P.is_inf) return P;

    int64_t num = modAdd(modMul(3, modSqr(P.x)), a);
    int64_t den = modMul(2, P.y);

    if (den == 0) {
        return (AffinePoint){0, 0, 1};
    }
    
    int64_t inv = modInv(den); 
    int64_t lambda = modMul(num, inv);

    int64_t x3 = modSub(modSqr(lambda), modMul(2, P.x));
    int64_t y3 = modSub(modMul(lambda, modSub(P.x, x3)), P.y);

    AffinePoint R = {x3, y3, 0};
    return R;
}

// Operasi penjumlahan titik pada kurva eliptik
AffinePoint affineAdd(AffinePoint P, AffinePoint Q) {
    if (P.is_inf) return Q;
    if (Q.is_inf) return P;
    if (P.x == Q.x && P.y == Q.y) return affineDouble(P);

    int64_t num = modSub(Q.y, P.y);
    int64_t den = modSub(Q.x, P.x);

    if (den == 0) {
        return (AffinePoint){0, 0, 1};
    }

    int64_t inv = modInv(den);
    int64_t lambda = modMul(num, inv);

    int64_t x3 = modSub(modSub(modSqr(lambda), P.x), Q.x);
    int64_t y3 = modSub(modMul(lambda, modSub(P.x, x3)), P.y);

    AffinePoint R = {x3, y3, 0};
    return R;
}

// Operasi perkalian skalar titik pada kurva eliptik
AffinePoint affineScalarMul(AffinePoint P, int64_t k) {
    AffinePoint R = {0, 0, 1};
    AffinePoint Temp = P;

    while (k > 0) {
        if (k % 2 == 1) {
            R = affineAdd(R, Temp);
        }
        Temp = affineDouble(Temp);
        k /= 2;
    }
    return R;
}

// Operasi perkalian double titik pada kurva eliptik dalam proyeksi berbobot
WeightedPoint weightedDouble(WeightedPoint P) {
    if (P.is_inf) return P;

    int64_t X1 = P.X; int64_t Y1 = P.Y; int64_t Z1 = P.Z;

    int64_t A = modSqr(X1);
    int64_t B = modSqr(Y1);
    int64_t C = modSqr(B);
    
    int64_t Z1_sq = modSqr(Z1);
    int64_t Z1_4 = modSqr(Z1_sq);
    int64_t M = modAdd(modMul(3, A), modMul(a, Z1_4));

    int64_t Z3 = modMul(2, modMul(Y1, Z1));
    int64_t S = modMul(4, modMul(X1, B));
    int64_t X3 = modSub(modSqr(M), modMul(2, S));
    int64_t Y3 = modSub(modMul(M, modSub(S, X3)), modMul(8, C));

    WeightedPoint R = {X3, Y3, Z3, 0};
    return R;
}

// Operasi penjumlahan titik pada kurva eliptik dalam proyeksi berbobot
WeightedPoint weightedAdd(WeightedPoint P, WeightedPoint Q) {
    if (P.is_inf) return Q;
    if (Q.is_inf) return P;

    int64_t Z1_sq = modSqr(P.Z);
    int64_t Z2_sq = modSqr(Q.Z);

    int64_t U1 = modMul(P.X, Z2_sq);
    int64_t U2 = modMul(Q.X, Z1_sq);

    int64_t S1 = modMul(P.Y, modMul(Q.Z, Z2_sq)); 
    int64_t S2 = modMul(Q.Y, modMul(P.Z, Z1_sq)); 

    if (U1 == U2) {
        if (S1 != S2) return (WeightedPoint){0,0,0,1};
        return weightedDouble(P);
    }

    int64_t H = modSub(U2, U1);
    int64_t R = modSub(S2, S1);

    int64_t Z3 = modMul(modMul(P.Z, Q.Z), H);

    int64_t H_sq = modSqr(H);
    int64_t H_cu = modMul(H, H_sq);
    int64_t term = modMul(2, modMul(U1, H_sq));
    int64_t X3 = modSub(modSub(modSqr(R), H_cu), term);

    int64_t Y3 = modSub(modMul(R, modSub(modMul(U1, H_sq), X3)), modMul(S1, H_cu));

    WeightedPoint Res = {X3, Y3, Z3, 0};
    return Res;
}

// Operasi perkalian skalar titik pada kurva eliptik dengan
WeightedPoint weightedScalarMul(WeightedPoint P, int64_t k) {
    WeightedPoint R = {0, 1, 0, 1};
    WeightedPoint Temp = P;

    while (k > 0) {
        if (k % 2 == 1) {
            R = weightedAdd(R, Temp);
        }
        Temp = weightedDouble(Temp);
        k /= 2;
    }
    return R;
}

// Konversi antara koordinat affine dan proyeksi berbobot
WeightedPoint toWeighted(AffinePoint P) {
    return (WeightedPoint){P.x, P.y, 1, P.is_inf};
}

// Konversi antara koordinat proyeksi berbobot dan affine
AffinePoint toAffine(WeightedPoint P) {
    if (P.is_inf || P.Z == 0) return (AffinePoint){0, 0, 1};
    
    int64_t Z_inv = modInv(P.Z);
    int64_t Z_inv_sq = modSqr(Z_inv);
    int64_t Z_inv_cu = modMul(Z_inv, Z_inv_sq);

    int64_t x = modMul(P.X, Z_inv_sq);
    int64_t y = modMul(P.Y, Z_inv_cu);

    return (AffinePoint){x, y, 0};
}

// Fungsi utama dari program
int main() {
    AffinePoint G;
    int64_t k;
    int iterasi = 2000;

    printf("Masukkan Bilangan Prima (p): ");
    scanf("%lld", &p);
    
    printf("Masukkan Parameter Kurva a: ");
    scanf("%lld", &a);

    printf("Masukkan Parameter Kurva b: ");
    scanf("%lld", &b);

    printf("Masukkan Titik Generator x: ");
    scanf("%lld", &G.x);

    printf("Masukkan Titik Generator y: ");
    scanf("%lld", &G.y);
    G.is_inf = 0;

    printf("Masukkan Nilai Skalar k: ");
    scanf("%lld", &k);

    printf("\nMelakukan benchmark %d iterasi...\n", iterasi);

    clock_t start = clock();
    AffinePoint ResultA;
    for(int i=0; i<iterasi; i++) {
        ResultA = affineScalarMul(G, k);
    }
    clock_t end = clock();
    double time_affine = ((double)(end - start)) / CLOCKS_PER_SEC;

    clock_t start_w = clock();
    AffinePoint ResultB;
    for(int i=0; i<iterasi; i++) {
        WeightedPoint W = toWeighted(G);
        WeightedPoint ResW = weightedScalarMul(W, k);
        ResultB = toAffine(ResW); 
    }
    clock_t end_w = clock();
    double time_weighted = ((double)(end_w - start_w)) / CLOCKS_PER_SEC;

    printf("\n=== HASIL PERHITUNGAN ===\n");
    printf("Metode Affine              : (%lld, %lld)\n", ResultA.x, ResultA.y);
    printf("Metode Proyeksi Berbobot   : (%lld, %lld)\n", ResultB.x, ResultB.y);

    printf("\n=== PERBANDINGAN WAKTU ===\n");
    printf("Waktu Affine               : %.6f detik\n", time_affine);
    printf("Waktu Proyeksi Berbobot    : %.6f detik\n", time_weighted);
    
    if (time_weighted > 0) {
        printf("Peningkatan Kecepatan      : %.2fx\n", time_affine / time_weighted);
    } else if (time_affine > 0 && time_weighted == 0) {
        printf("Peningkatan Kecepatan      : > 100x (Waktu Weighted mendekati 0)\n");
    }

    return 0;
}