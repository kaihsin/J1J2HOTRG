#ifndef UNI10_TESTING_TOOL_QDR_H_
#define UNI10_TESTING_TOOL_QDR_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_qdr(bool inplace)
  {
    uni10_uint64 Mstart = 1000, Mend = 10000, Mstp = 1000;
    uni10_uint64 Nstart = 1000,               Nstp = 1000;
    int          ntest = 10;

    cout << "          |               Running Time             |             Error            \n"
         << "    M    N|     WTime   CPUTime   GPUTime   MPITime| |I-Q*Q^T|     |D*R| |M-Q*D*R|\n"
         << "==================================================================================\n";

    for (uni10_uint64 M = Mstart; M <= Mend; M += Mstp) {
      for (uni10_uint64 N = Nstart; N <= M; N += Nstp) {
        vector<double> wtimes;
        vector<double> cputimes;
        vector<vector<double>> errors(3);

        for (int itest = 0; itest < ntest; ++itest) {
          // Assign random elements to the pending matrix.
          Matrix<uni10_type> matrix(M, N);
          uni10_rand(matrix, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

          // Declare result.
          vector<Matrix<uni10_type>> QDR(3);

          // Start the timer.
          wtimes.push_back(uni10_wtime());
          cputimes.push_back(uni10_cputime());

          // Perform QDR decomposition using uni10.
          if (inplace) qdr(matrix, QDR.at(0), QDR.at(1), QDR.at(2), INPLACE);
          else         QDR = qdr(matrix);

          // End the timer.
          wtimes.back() = uni10_wtime() - wtimes.back();
          cputimes.back() = uni10_cputime() - cputimes.back();

          // Error1 = |I - Q*Q^T|
          Matrix<uni10_type> QxQT = dot(dagger(QDR[0]), QDR[0]);
          Matrix<uni10_type> I(N, N);
          I.identity();
          errors.at(0).push_back(norm(I - QxQT));

          // Error2 = |ltDR| where ltR is the lower triangular part of D * R
          Matrix<uni10_type> ltDxR(dot(QDR[1], QDR[2]));
          for (decltype(ltDxR.row()) i = 0; i < ltDxR.row(); ++i) {
            for (decltype(ltDxR.col()) j = i; j < ltDxR.col(); ++j) {
              ltDxR[i * ltDxR.row() + j] = 0;
            }
          }
          errors.at(1).push_back(norm(ltDxR));

          // Error3 = |M-Q*D*R|
          errors.at(2).push_back(norm(matrix - dot(dot(QDR[0], QDR[1]), QDR[2])));
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
