add_executable(max_potential_on_chain_test max_potential_on_chain_test.cpp)
target_link_libraries(max_potential_on_chain_test LP_MP-horizon-tracking DD_ILP lingeling m stdc++ pthread )

add_executable(max_potential_on_tree_test max_potential_on_tree_test.cpp)
target_link_libraries(max_potential_on_tree_test LP_MP-horizon-tracking DD_ILP lingeling m stdc++ pthread )

add_executable(max_potential_uai_chain_test max_potential_uai_chain_test.cpp)
target_link_libraries(max_potential_uai_chain_test LP_MP-horizon-tracking DD_ILP FW-MAP arboricity lingeling m stdc++ pthread  )

add_executable(max_potential_uai_grid_test max_potential_uai_grid_test.cpp)
target_link_libraries(max_potential_uai_grid_test LP_MP-horizon-tracking DD_ILP FW-MAP arboricity lingeling m stdc++ pthread  )
