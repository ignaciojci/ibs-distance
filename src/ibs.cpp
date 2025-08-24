// [[Rcpp::plugins(cpp11)]]
// Enable OpenMP if available; harmless if not supported.
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// popcount for 64-bit words
static inline uint32_t popcnt64(uint64_t x) {
  #if defined(__GNUG__) || defined(__clang__)
    return static_cast<uint32_t>(__builtin_popcountll(x));
  #else
  // portable fallback (Hacker's Delight)
  x = x - ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  return static_cast<uint32_t>((((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL) >> 56);
#endif
}

// Pack an integer column (markers) for one sample into bitplanes:
// b0: genotype==0, b1: genotype==1 (het), b2: genotype==2, bM: is.na
static void pack_sample_bitplanes(const int* col, int p,
                                  std::vector<uint64_t>& b0,
                                  std::vector<uint64_t>& b1,
                                  std::vector<uint64_t>& b2,
                                  std::vector<uint64_t>& bM) {
  const int W = (p + 63) / 64;
  b0.assign(W, 0ULL);
  b1.assign(W, 0ULL);
  b2.assign(W, 0ULL);
  bM.assign(W, 0ULL);

  for (int site = 0; site < p; ++site) {
    int word = site >> 6;
    int bit  = site & 63;
    int g = col[site];
    uint64_t mask = (1ULL << bit);

    if (g == 0)      b0[word] |= mask;
    else if (g == 1) b1[word] |= mask;
    else if (g == 2) b2[word] |= mask;
    else             bM[word] |= mask; // NA/other
  }
}

// Compute IBS for one pair (i,j) using bitplanes, returns {dist, sites}
// IBS identity = (same - 0.5*hets) / sites, dist = 1 - identity.
// same counts (0-0), (2-2), and (1-1)
// diff counts (0-2), (2-0), (het vs hom), and (het-het)
// hets counts (1-1).
static inline std::pair<double,double>
pair_ibs(const std::vector<uint64_t>& b0i,
         const std::vector<uint64_t>& b1i,
         const std::vector<uint64_t>& b2i,
         const std::vector<uint64_t>& bMi,
         const std::vector<uint64_t>& b0j,
         const std::vector<uint64_t>& b1j,
         const std::vector<uint64_t>& b2j,
         const std::vector<uint64_t>& bMj) {

  const int W = static_cast<int>(b0i.size());
  uint64_t same_bits, diff_bits, half_bits, het_bits, valid;
  long long sameCnt = 0, diffCnt = 0, halfCnt = 0, hetCnt = 0, sites = 0;

  for (int w = 0; w < W; ++w) {
    // sites valid in both
    valid = ~(bMi[w] | bMj[w]);

    // same genotype (includes het-het)
    same_bits = ((b0i[w] & b0j[w]) | (b2i[w] & b2j[w]) | (b1i[w] & b1j[w])) & valid;

    // different genotype:
    // hom-hom different, het-hom, hom-het, and also include het-het so that
    // same + diff - het counts one site for het-het
    diff_bits = (
      (b0i[w] & b2j[w]) | (b2i[w] & b0j[w]) |
      (b1i[w] & b1j[w])
    ) & valid;

    half_bits = (
      (b1i[w] & (b0j[w] | b2j[w])) |
      (b1j[w] & (b0i[w] | b2i[w]))
    ) & valid;
    
    // het-het
    het_bits = (b1i[w] & b1j[w]) & valid;

    sameCnt += popcnt64(same_bits);
    diffCnt += popcnt64(diff_bits);
    halfCnt += popcnt64(half_bits);
    hetCnt  += popcnt64(het_bits);
  }

  sites = sameCnt + diffCnt + halfCnt - hetCnt; // ensure het-het counts once
  if (sites <= 0) return {NAN, 0.0};

  double identity = ( (double)sameCnt - 0.5 * (double)halfCnt ) / (double)sites;
  double dist = 1.0 - identity;
  return {dist, (double)sites};
}

// [[Rcpp::export]]
NumericMatrix fast_ibs(SEXP M_,
                                 int block_size = 256,
                                 bool diagonal_trueIBS = false,
                                 int min_sites = 0,
                                 int n_threads = 0) {
  // Expect integer matrix with values 0/1/2 and NA = missing
  IntegerMatrix M(M_);
  const int n = M.nrow();
  const int p = M.ncol();
  if (n == 0 || p == 0) stop("Empty matrix.");

  // Prepack all samples into bitplanes
  const int W = (p + 63) / 64;
  std::vector< std::vector<uint64_t> > B0(n), B1(n), B2(n), BM(n);

  // Pack columns by sample (rows are samples; transpose view to get a sample pointer)
  for (int i = 0; i < n; ++i) {
    // M(i, :) is contiguous with stride n: we copy into a temp column buffer
    std::vector<int> buf(p);
    for (int s = 0; s < p; ++s) buf[s] = M(i, s);
    B0[i].reserve(W); B1[i].reserve(W); B2[i].reserve(W); BM[i].reserve(W);
    pack_sample_bitplanes(buf.data(), p, B0[i], B1[i], B2[i], BM[i]);
  }

  NumericMatrix D(n, n);

  // Optionally set OpenMP threads
  #ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
  #endif

  // Block over i to improve cache locality
  const int B = std::max(1, block_size);

  for (int ib = 0; ib < n; ib += B) {
    int i_end = std::min(n, ib + B);

    // Parallelize the inner i-loop; each (i, *) row is independent
    // Parallelize only if OpenMP is available
    #ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = ib; i < i_end; ++i) {
      for (int j = i; j < n; ++j) {
        if (i == j && !diagonal_trueIBS) {
          D(i,i) = 0.0;
          continue;
        }
        auto pr = pair_ibs(B0[i], B1[i], B2[i], BM[i],
                           B0[j], B1[j], B2[j], BM[j]);
        double dist = pr.first;
        double sites = pr.second;
        if (sites < (double)min_sites) dist = NA_REAL;
        D(i,j) = dist;
        D(j,i) = dist;
      }
    }
  }
  return D;
}
