#ifndef UNI10_TESTING_TOOL_QDR_CPIVOT_H_
#define UNI10_TESTING_TOOL_QDR_CPIVOT_H_

#include "uni10/uni10.hpp"

#include <iostream>
#include <iomanip>

#include "testing_timer.h"

using namespace std;
using namespace uni10;

template <typename uni10_type>
  int testing_qdr_cpivot(bool inplace)
  {
    uni10_uint64 Nstart = 100, Nend = 1000, Nstp = 100;
    int          ntest = 1;

    cout << "        |               Running Time             |             Error              \n"
         << "       N|     WTime   CPUTime   GPUTime   MPITime| |I-Q*Q^T|     D     R |M-Q*D*R|\n"
         << "==================================================================================\n";

    for (uni10_uint64 N = Nstart; N <= Nend; N += Nstp) {
      vector<double> wtimes;
      vector<double> cputimes;
      vector<vector<double>> errors(2);
      vector<bool> isCorrects(2);

      for (int itest = 0; itest < ntest; ++itest) {
        // Assign random elements to the pending matrix.
        Matrix<uni10_type> matrix(N, N);
        uni10_rand(matrix, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);

        // Declare result.
        vector<Matrix<uni10_type>> QDR(3);

        // Start the timer.
        wtimes.push_back(uni10_wtime());
        cputimes.push_back(uni10_cputime());

        // Perform QDR_CPIVOT decomposition using uni10.
        if (inplace) qdr_cpivot(matrix, QDR.at(0), QDR.at(1), QDR.at(2), INPLACE);
        else         QDR = qdr_cpivot(matrix);

        // End the timer.
        wtimes.back() = uni10_wtime() - wtimes.back();
        cputimes.back() = uni10_cputime() - cputimes.back();

        // Error1 = |I - Q*Q^T|
        Matrix<uni10_type> QxQT = dot(dagger(QDR[0]), QDR[0]);
        Matrix<uni10_type> I(N, N);
        I.identity();
        errors.at(0).push_back(norm(I - QxQT));

        // Check if D is Diagnal and the singular value is decreasing
        isCorrects.at(0) = true;
        for (decltype(QDR.at(1).row()) i = 1;
              i < min(QDR.at(1).row(), QDR.at(1).col()); ++i)
          if (abs(QDR.at(1)[i]) > abs(QDR.at(1)[i - 1])) isCorrects.at(0) = false;

        // Check if R can be upper triangular matrix
        isCorrects.at(1) = true;
        vector<bool> zeroCnt(QDR.at(2).col() + 1, false);
        for (decltype(QDR.at(2).col()) i = 0;
            i < QDR.at(2).col() && isCorrects.at(1); ++i) {
          int count = 0;
          auto elem_ind = i;
          while (count < QDR.at(2).row() && QDR.at(2)[elem_ind] != 0.0) {
            elem_ind += QDR.at(2).col();
            count++;
          }
          while (count < QDR.at(2).row() && QDR.at(2)[elem_ind] == 0.0) {
            elem_ind += QDR.at(2).col();
            count++;
          }
          if (count != QDR.at(2).row()) isCorrects.at(1) = false;
        }

        // Error2 = |M-Q*D*R|
        errors.at(1).push_back(norm(matrix - dot(dot(QDR[0], QDR[1]), QDR[2])));
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
      cout << setw(8) << N << "|"
           << setw(10) << setprecision(3) << fixed << wtime_avg
           << setw(10) << setprecision(3) << fixed << cputime_avg
           << setw(10) << "--- "
           << setw(10) << "--- " << "|"
           << setw(10) << scientific << setprecision(3) << error_avg.at(0)
           << setw(6) << (isCorrects.at(0) ? "OK" : "failed")
           << setw(6) << (isCorrects.at(1) ? "OK" : "failed")
           << setw(10) << scientific << setprecision(3) << error_avg.at(1)
           << endl;
    }

    return 0;
  }

#endif
