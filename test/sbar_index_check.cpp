// Standalone check of the safe_core_per_buffer indexing in buildLocalSbarMap.
//
// The ParticleBalancer constructor fills safe_core_per_buffer[core_nents * nbuffers] from a
// per-buffer safe-zone exchange, and buildLocalSbarMap is its only reader. Each cell below
// carries a unique (buffer, element) tag, so any wrong index is detected.
//
//   master (strided store):      element e of buffer b at e*nbuffers + b, read as i*nbuffers + j
//   8d6acb3 (contiguous store):  element e of buffer b at b*core_nents + e,
//                                read as j*nbuffers + i   <- the defect
//   this PR:                     read as j*nelms + i      (nelms == core_nents)
#include <vector>
#include <cstdio>
#include <cstdint>
int main() {
  long shapes[][2] = {{7,3},{698305,5},{1,4},{40,40},{1000,1},{13,7},{731000,8},{592923,16}};
  int nshapes = sizeof(shapes)/sizeof(shapes[0]);
  auto tag = [](int b, long e) -> uint64_t { return (uint64_t)(b+1) * 100000007ULL + (uint64_t)(e+1); };
  long total_corr = 0, total_bug = 0;
  for (int s = 0; s < nshapes; ++s) {
    long core_nents = shapes[s][0];
    int  nbuffers   = (int)shapes[s][1];
    std::vector<uint64_t> A_strided((size_t)core_nents*nbuffers, 0);
    std::vector<uint64_t> A_contig ((size_t)core_nents*nbuffers, 0);
    for (int b = 0; b < nbuffers; ++b)
      for (long e = 0; e < core_nents; ++e) {
        uint64_t v = tag(b, e);
        A_strided[(size_t)e*nbuffers + b]  = v;   // master strided store
        A_contig [(size_t)b*core_nents + e] = v;  // 8d6acb3 contiguous store
      }
    long mism_corr = 0, mism_bug = 0;
    for (long i = 0; i < core_nents; ++i)
      for (int j = 0; j < nbuffers; ++j) {
        uint64_t base = A_strided[(size_t)i*nbuffers + j];   // master consumer reads cell (j,i)
        uint64_t corr = A_contig [(size_t)j*core_nents + i]; // corrected index (this PR)
        size_t bug_index = (size_t)j*nbuffers + i;           // 8d6acb3 index
        if (corr != base) ++mism_corr;
        // An out-of-bounds 8d6acb3 index counts as a mismatch without being dereferenced.
        if (bug_index >= A_contig.size() || A_contig[bug_index] != base) ++mism_bug;
      }
    printf("core_nents=%-7ld nbuffers=%-3d : corrected_mismatches=%-9ld 8d6acb3_mismatches=%ld\n",
           core_nents, nbuffers, mism_corr, mism_bug);
    total_corr += mism_corr; total_bug += mism_bug;
  }
  printf("\nTOTAL corrected mismatches = %ld\n", total_corr);
  printf("TOTAL 8d6acb3 mismatches   = %ld\n", total_bug);
  return total_corr == 0 ? 0 : 1;
}
