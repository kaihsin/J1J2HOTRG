#ifndef UNI10_TESTING_TOOL_QR_H_
#define UNI10_TESTING_TOOL_QR_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_qr(bool inplace)
  {
    uni10_uint64 Mstart = 1000, Mend = 10000, Mstp = 1000;
    uni10_uint64 Nstart = 1000,               Nstp = 1000;
    int          ntest = 10;

    cout << "          |               Running Time             |             Error            \n"
         << "    M    N|     WTime   CPUTime   GPUTime   MPITime| |I-Q*Q^T|       |R|   |M-Q*R|\n"
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
          vector<Matrix<uni10_type>> QR(2);

          // Start the timer.
          wtimes.push_back(uni10_wtime());
          cputimes.push_back(uni10_cputime());

          // Perform QR decomposition using uni10.
          if (inplace) qr(matrix, QR.at(0), QR.at(1), INPLACE);
          else         QR = qr(matrix);

          // End the timer.
          wtimes.back() = uni10_wtime() - wtimes.back();
          cputimes.back() = uni10_cputime() - cputimes.back();

          // Error1 = |I - Q*Q^T|
          Matrix<uni10_type> QxQT = dot(dagger(QR[0]), QR[0]);
          Matrix<uni10_type> I(N, N);
          I.identity();
          errors.at(0).push_back(norm(I - QxQT));

          // Error2 = |ltR| where ltR is the lower triangular part of R
          Matrix<uni10_type> ltR(QR[1]);
          for (decltype(ltR.row()) i = 0; i < ltR.row(); ++i) {
            for (decltype(ltR.col()) j = i; j < ltR.col(); ++j) {
              ltR[i * ltR.row() + j] = 0;
            }
          }
          errors.at(1).push_back(norm(ltR));

          // Error3 = |M - Q*R|
          errors.at(2).push_back(norm(matrix - dot(QR[0], QR[1])));
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
