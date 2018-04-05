#include "test.h"
#include "horizon_tracking_factors.hxx"

using namespace LP_MP;

int main()
{
  int numNodes = 3;
  std::vector<INDEX> numLabels = {3, 2, 3};

  std::vector<MaxPairwisePotential> maxPairwisePotentials;
  std::vector<matrix<REAL>> LinearPairwisePotentials(numNodes - 1);
  for (int i = 0; i < numNodes - 1; i++)
  {
    LinearPairwisePotentials[i] = matrix<REAL>(numLabels[i], numLabels[i + 1]);
  }
  
  MaxPairwisePotential pot;
  
  pot.n1 = 0;
  pot.n2 = 1;
  pot.l1 = 0;
  pot.l2 = 0;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 1;
  pot.value = 2;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 0;
  pot.n2 = 1;
  pot.l1 = 0;
  pot.l2 = 1;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 2;
  pot.value = 3;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 0;
  pot.n2 = 1;
  pot.l1 = 1;
  pot.l2 = 0;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 3;
  pot.value = 2;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 0;
  pot.n2 = 1;
  pot.l1 = 1;
  pot.l2 = 1;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 4;
  pot.value = 1;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 0;
  pot.n2 = 1;
  pot.l1 = 2;
  pot.l2 = 0;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 1;
  pot.value = 1;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 0;
  pot.n2 = 1;
  pot.l1 = 2;
  pot.l2 = 1;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 2;
  pot.value = 2;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 1;
  pot.n2 = 2;
  pot.l1 = 0;
  pot.l2 = 0;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 1.5;
  pot.value = 2;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 1;
  pot.n2 = 2;
  pot.l1 = 0;
  pot.l2 = 1;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 5;
  pot.value = 3;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 1;
  pot.n2 = 2;
  pot.l1 = 0;
  pot.l2 = 2;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 1;
  pot.value = 3;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 1;
  pot.n2 = 2;
  pot.l1 = 1;
  pot.l2 = 0;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 3;
  pot.value = 6;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 1;
  pot.n2 = 2;
  pot.l1 = 1;
  pot.l2 = 1;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 1;
  pot.value = 1;

  maxPairwisePotentials.push_back(pot);

  pot.n1 = 1;
  pot.n2 = 2;
  pot.l1 = 1;
  pot.l2 = 2;
  LinearPairwisePotentials[pot.n1](pot.l1, pot.l2) = 6;
  pot.value = 0;

  maxPairwisePotentials.push_back(pot);

  max_potential_on_chain chain = max_potential_on_chain(maxPairwisePotentials, LinearPairwisePotentials, numLabels, numNodes);
  REAL objective = chain.LowerBound();
  test(objective == 4.5);
}
