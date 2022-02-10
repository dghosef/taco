#include "pochivm.h"
#include <map>
#include <string>
#include <functional>
std::map<std::string, std::function<PochiVM::ValueVT(std::vector<PochiVM::ValueVT>&)>> fns = {
  std::make_pair("taco_binarySearchAfter", [](std::vector<PochiVM::ValueVT> &params) -> PochiVM::ValueVT {
  return PochiVM::CallFreeFn::taco_binarySearchAfter_wrapper(PochiVM::StaticCast<int *>(params[0]), PochiVM::StaticCast<int>(params[1]), PochiVM::StaticCast<int>(params[2]), PochiVM::StaticCast<int>(params[3]));
}),
  std::make_pair("taco_binarySearchBefore", [](std::vector<PochiVM::ValueVT> &params) -> PochiVM::ValueVT {
  return PochiVM::CallFreeFn::taco_binarySearchBefore_wrapper(PochiVM::StaticCast<int *>(params[0]), PochiVM::StaticCast<int>(params[1]), PochiVM::StaticCast<int>(params[2]), PochiVM::StaticCast<int>(params[3]));
})
};
