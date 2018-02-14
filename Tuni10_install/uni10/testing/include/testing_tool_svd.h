#ifndef UNI10_TESTING_TOOL_SVD_H_
#define UNI10_TESTING_TOOL_SVD_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_svd(bool inplace)
  {
    uni10_uint64 MNsumstart = 2000, MNsumend = 10000, MNsumstp = 1000;
    uni10_uint64 Mstart     = 1000,                   Mstp     = 1000;
    int          ntest      = 10;

    cout << "          |               Running Time             |             Error            \n"
         << "    M    N|     WTime   CPUTime   GPUTime   MPITime| |I-U*U^T| |I-V*V^T||M-U*S*VT|\n"
         << "==================================================================================\n";

    for (uni10_uint64 MNsum = MNsumstart; MNsum <= MNsumend; MNsum += MNsumstp) {
      for (uni10_uint64 M = Mstart, N = MNsum - M; N > 0; M += Mstp, N -= Mstp) {
        vector<double> wtimes;
        vector<double> cputimes;
        vector<vector<double>> errors(3);

        for (int itest = 0; itest < ntest; ++itest) {
          // Assign random elements to the pending matrix.
          Matrix<uni10_type> matrix(M, N);
          uni10_rand(matrix, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

          // Declare results.
          vector<Matrix<uni10_type>> SVD(3);

          // Start the timer.
          wtimes.push_back(uni10_wtime());
          cputimes.push_back(uni10_cputime());

          // Perform SVD decomposition using uni10.
          if (inplace) svd(matrix, SVD.at(0), SVD.at(1), SVD.at(2), INPLACE);
          else         SVD = svd(matrix);

          // End the timer.
          wtimes.back() = uni10_wtime() - wtimes.back();
          cputimes.back() = uni10_cputime() - cputimes.back();

          // Error1 = |I - U*U^T|
          Matrix<uni10_type> UxUT =
            SVD[0].row() < SVD[0].col() ? dot(SVD[0], dagger(SVD[0]))
                                        : dot(dagger(SVD[0]), SVD[0]);
          Matrix<uni10_type> I(min(M, N), min(M, N));
          I.identity();
          errors.at(0).push_back(norm(I - UxUT));

          // Error2 = |I - V*V^T|
          Matrix<uni10_type> VxVT =
            SVD[2].row() < SVD[2].col() ? dot(SVD[2], dagger(SVD[2]))
                                        : dot(dagger(SVD[2]), SVD[2]);
          errors.at(1).push_back(norm(I - VxVT));

          // Error3 = |M - U*S*V^T|
          errors.at(2).push_back(norm(matrix - dot(dot(SVD[0], SVD[1]), SVD[2])));
        }

        // Calculate the averages.
        double wtime_avg =
          accumulate(wtimes.begin(), wtimes.end(), 0.0) / wtimes.size();
        double cputime_avg =
          accumulate(cputimes.begin(), cputimes.end(), 0.0) / cputimes.size();
        vector<double> error_avg;
        for (auto &error : errors)
          error_avg.push_back(accumulate(error.begin(), error.end(), 0.0) / error.size());

        // Print performance and error infos.
        cout << setw(5) << M << setw(5) << N << "|"
             << setw(10) << setprecision(3) << fixed << wtime_avg
             << setw(10) << setprecision(3) << fixed << cputime_avg
             << setw(10) << "--- "
             << setw(10) << "--- " << "|"
             << setw(10) << scientific << setprecision(3) << error_avg.at(0)
             << setw(10) << scientific << setprecision(3) << error_avg.at(1)
             << setw(10) << scientific << setprecision(3) << error_avg.at(2)
             << endl;
      }
    }

    return 0;
  }

#endif
