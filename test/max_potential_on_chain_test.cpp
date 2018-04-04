#include "test.h"
#include "horizon_tracking_factors.hxx"

using namespace LP_MP;

int main()
{
  std::vector<max_linear_pairwise_pot> allPairwisePots;
  max_linear_pairwise_pot pot;

  pot.i = 0;
  pot.j = 1;
  pot.index_mini_node1 = 0;
  pot.index_mini_node2 = 3;
  pot.l1 = 0;
  pot.l2 = 0;
  pot.linear_pot = 1;
  pot.max_pot = 2;

  allPairwisePots.push_back(pot);

  pot.i = 0;
  pot.j = 1;
  pot.index_mini_node1 = 0;
  pot.index_mini_node2 = 4;
  pot.l1 = 0;
  pot.l2 = 1;
  pot.linear_pot = 2;
  pot.max_pot = 3;

  allPairwisePots.push_back(pot);

  pot.i = 0;
  pot.j = 1;
  pot.index_mini_node1 = 1;
  pot.index_mini_node2 = 3;
  pot.l1 = 1;
  pot.l2 = 0;
  pot.linear_pot = 3;
  pot.max_pot = 2;

  allPairwisePots.push_back(pot);

  pot.i = 0;
  pot.j = 1;
  pot.index_mini_node1 = 1;
  pot.index_mini_node2 = 4;
  pot.l1 = 1;
  pot.l2 = 1;
  pot.linear_pot = 4;
  pot.max_pot = 1;

  allPairwisePots.push_back(pot);

  pot.i = 0;
  pot.j = 1;
  pot.index_mini_node1 = 2;
  pot.index_mini_node2 = 3;
  pot.l1 = 2;
  pot.l2 = 0;
  pot.linear_pot = 1;
  pot.max_pot = 1;

  allPairwisePots.push_back(pot);

  pot.i = 0;
  pot.j = 1;
  pot.index_mini_node1 = 2;
  pot.index_mini_node2 = 4;
  pot.l1 = 2;
  pot.l2 = 1;
  pot.linear_pot = 2;
  pot.max_pot = 2;

  allPairwisePots.push_back(pot);

  pot.i = 1;
  pot.j = 2;
  pot.index_mini_node1 = 3;
  pot.index_mini_node2 = 5;
  pot.l1 = 0;
  pot.l2 = 0;
  pot.linear_pot = 2;
  pot.max_pot = 2;

  allPairwisePots.push_back(pot);

  pot.i = 1;
  pot.j = 2;
  pot.index_mini_node1 = 3;
  pot.index_mini_node2 = 6;
  pot.l1 = 0;
  pot.l2 = 1;
  pot.linear_pot = 5;
  pot.max_pot = 3;

  allPairwisePots.push_back(pot);

  pot.i = 1;
  pot.j = 2;
  pot.index_mini_node1 = 3;
  pot.index_mini_node2 = 7;
  pot.l1 = 0;
  pot.l2 = 2;
  pot.linear_pot = 1;
  pot.max_pot = 3;

  allPairwisePots.push_back(pot);

  pot.i = 1;
  pot.j = 2;
  pot.index_mini_node1 = 4;
  pot.index_mini_node2 = 5;
  pot.l1 = 1;
  pot.l2 = 0;
  pot.linear_pot = 3;
  pot.max_pot = 6;

  allPairwisePots.push_back(pot);

  pot.i = 1;
  pot.j = 2;
  pot.index_mini_node1 = 4;
  pot.index_mini_node2 = 6;
  pot.l1 = 1;
  pot.l2 = 1;
  pot.linear_pot = 1;
  pot.max_pot = 1;

  allPairwisePots.push_back(pot);

  pot.i = 1;
  pot.j = 2;
  pot.index_mini_node1 = 4;
  pot.index_mini_node2 = 7;
  pot.l1 = 1;
  pot.l2 = 2;
  pot.linear_pot = 6;
  pot.max_pot = 0;

  allPairwisePots.push_back(pot);

  max_potential_on_chain chain = max_potential_on_chain(allPairwisePots, 3, 8);
  chain.Solve();
  REAL objective = chain.GetSolutionObjective();
  test(objective == 5);
}
