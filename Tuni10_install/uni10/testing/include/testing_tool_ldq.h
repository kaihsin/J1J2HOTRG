#ifndef UNI10_TESTING_TOOL_LDQ_H_
#define UNI10_TESTING_TOOL_LDQ_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_ldq(bool inplace)
  {
    uni10_uint64 Mstart = 1000,               Mstp = 1000;
    uni10_uint64 Nstart = 1000, Nend = 10000, Nstp = 1000;
    int          ntest = 10;

    cout << "          |               Running Time             |             Error            \n"
         << "    M    N|     WTime   CPUTime   GPUTime   MPITime|     |L*D| |I-Q*Q^T| |M-L*D*Q|\n"
         << "==================================================================================\n";

    for (uni10_uint64 N = Nstart; N <= Nend; N += Nstp) {
      for (uni10_uint64 M = Mstart; M <= N; M += Mstp) {
        vector<double> wtimes;
        vector<double> cputimes;
        vector<vector<double>> errors(3);

        for (int itest = 0; itest < ntest; ++itest) {
          // Assign random elements to the pending matrix.
          Matrix<uni10_type> matrix(M, N);
          uni10_rand(matrix, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

          // Declare results.
          vector<Matrix<uni10_type>> LDQ(3);

          // Start the timer.
          wtimes.push_back(uni10_wtime());
          cputimes.push_back(uni10_cputime());

          // Perform LDQ decomposition using uni10.
          if (inplace) ldq(matrix, LDQ.at(0), LDQ.at(1), LDQ.at(2), INPLACE);
          else         LDQ = ldq(matrix);

          // End the timer.
          wtimes.back() = uni10_wtime() - wtimes.back();
          cputimes.back() = uni10_cputime() - cputimes.back();

          // Error1 = |utLxD| where utL is the upper triangular part of L*D
          Matrix<uni10_type> utLxD(dot(LDQ[0], LDQ[1]));
          for (decltype(utLxD.row()) i = 0; i < utLxD.row(); ++i) {
            for (decltype(utLxD.col()) j = 0; j <= i; ++j) {
              utLxD[i * utLxD.row() + j] = 0;
            }
          }
          errors.at(0).push_back(norm(utLxD));

          // Error2 = |I - Q*Q^T|
          Matrix<uni10_type> QxQT = dot(LDQ[2], dagger(LDQ[2]));
          Matrix<uni10_type> I(M, M);
          I.identity();
          errors.at(1).push_back(norm(I - QxQT));

          // Error3 = |M-L*D*Q|
          errors.at(2).push_back(norm(matrix - dot(dot(LDQ[0], LDQ[1]), LDQ[2])));
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
