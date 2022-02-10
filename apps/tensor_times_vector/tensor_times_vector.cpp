#include "taco.h"
#include <iostream>
#include <cmath>
#include <string>
#include <chrono>
#include <iostream>

bool debug = true;
class AutoProfiler {
 public:
  AutoProfiler(std::string name)
      : m_name(std::move(name)),
        m_beg(std::chrono::high_resolution_clock::now()) { }
  ~AutoProfiler() {
    auto end = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - m_beg);
    std::cout << m_name << " : " << dur.count() << " musec\n";
  }
 private:
  std::string m_name;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_beg;
};

using namespace taco;

int main(int argc, char *argv[]) {
  Format csr({Dense, Sparse});
  Format csf({Sparse, Sparse, Sparse});
  Format sv({Sparse});

  Tensor<double> A("A", {2, 3}, csr);
  Tensor<double> B("B", {2, 3, 4}, csf);
  Tensor<double> c("c", {4}, sv);

  // Insert data into B and c
  B(0,0,0) = 1.0;
  B(1,2,0) = 2.0;
  B(1,2,1) = 3.0;
  c(0) = 4.0;
  c(1) = 5.0;

  IndexVar i, j, k;
  A(i, j) = B(i, j, k) * c(k);
  AutoProfiler ap("stuff");
  A.compile();
  A.evaluate();
  std::cout << A << std::endl;
}
