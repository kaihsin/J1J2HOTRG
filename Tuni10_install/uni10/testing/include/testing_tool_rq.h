#ifndef UNI10_TESTING_TOOL_RQ_H_
#define UNI10_TESTING_TOOL_RQ_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_rq(bool inplace)
  {
    uni10_uint64 Mstart = 1000,               Mstp = 1000;
    uni10_uint64 Nstart = 1000, Nend = 10000, Nstp = 1000;
    int          ntest = 10;

    cout << "          |               Running Time             |             Error            \n"
         << "    M    N|     WTime   CPUTime   GPUTime   MPITime|       |R| |I-Q*Q^T|   |M-R*Q|\n"
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
          vector<Matrix<uni10_type>> RQ;

          // Start the timer.
          wtimes.push_back(uni10_wtime());
          cputimes.push_back(uni10_cputime());

          // Perform RQ decomposition using uni10.
          if (inplace) rq(matrix, RQ.at(0), RQ.at(1), INPLACE);
          else         RQ = rq(matrix);

          // End the timer.
          wtimes.back() = uni10_wtime() - wtimes.back();
          cputimes.back() = uni10_cputime() - cputimes.back();

          // Error1 = |ltR| where ltR is the lower triangular part of R
          Matrix<uni10_type> ltR(RQ[0]);
          for (decltype(ltR.row()) i = 0; i < ltR.row(); ++i) {
            for (decltype(ltR.col()) j = i; j < ltR.col(); ++j) {
              ltR[i * ltR.row() + j] = 0;
            }
          }
          errors.at(0).push_back(norm(ltR));

          // Error2 = |I - Q*Q^T|
          Matrix<uni10_type> QxQT = dot(RQ[1], dagger(RQ[1]));
          Matrix<uni10_type> I(M, M);
          I.identity();
          errors.at(1).push_back(norm(I - QxQT));


          // Error3 = |M - R*Q|
          errors.at(2).push_back(norm(matrix - dot(RQ[0], RQ[1])));
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
