#ifndef UNI10_TESTING_TOOL_INVERSE_H_
#define UNI10_TESTING_TOOL_INVERSE_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_inverse(bool inplace)
  {
    uni10_uint64 Nstart = 1000, Nend = 10000, Nstp = 1000;
    int          ntest = 10;

    cout << "          |               Running Time             |             Error            \n"
         << "         N|     WTime   CPUTime   GPUTime   MPITime|          |I-M*M^-1|          \n"
         << "==================================================================================\n";

    for (uni10_uint64 N = Nstart; N <= Nend; N += Nstp) {
      vector<double> wtimes;
      vector<double> cputimes;
      vector<double> errors;

      for (int itest = 0; itest < ntest; ++itest) {
        // Assign random elements to the pending matrix.
        Matrix<uni10_type> matrix(N, N);
        uni10_rand(matrix, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

        // Declare Result.
        Matrix<uni10_type> INV;
        // inplace version need to be initialized first
        if (inplace) INV = matrix;

        // Start the timer.
        wtimes.push_back(uni10_wtime());
        cputimes.push_back(uni10_cputime());

        // Perform matrix inversion using uni10.
        if (inplace) inverse(INV, INPLACE);
        else         INV = inverse(matrix);

        // End the timer.
        wtimes.back() = uni10_wtime() - wtimes.back();
        cputimes.back() = uni10_cputime() - cputimes.back();

        // Error = |I - M*M^-1|
        Matrix<uni10_type> MxM_1 = dot(matrix, INV);
        Matrix<uni10_type> I(N, N);
        I.identity();
        errors.push_back(norm(I - MxM_1));
      }

      // Calculate the averages.
      double wtime_avg =
        accumulate(wtimes.begin(), wtimes.end(), 0.0) / wtimes.size();
      double cputime_avg =
        accumulate(cputimes.begin(), cputimes.end(), 0.0) / cputimes.size();
      double error_avg =
        accumulate(errors.begin(), errors.end(), 0.0) / errors.size();

      // Print performance and error infos.
      cout << setw(5) << " " << setw(5) << N << "|"
           << setw(10) << setprecision(3) << fixed << wtime_avg
           << setw(10) << setprecision(3) << fixed << cputime_avg
           << setw(10) << "--- "
           << setw(10) << "--- " << "|"
           << setw(10) << " "
           << setw(10) << scientific << setprecision(3) << error_avg
           << setw(10) << " "
           << endl;
    }

    return 0;
  }

#endif
