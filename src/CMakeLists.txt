add_executable(horizon_tracking_fwmap horizon_tracking_fwmap.cpp)
target_link_libraries(horizon_tracking_fwmap LP_MP-horizon-tracking DD_ILP FW-MAP arboricity lingeling m stdc++ pthread )

add_executable(horizon_tracking_conic_bundle horizon_tracking_conic_bundle.cpp)
target_link_libraries(horizon_tracking_conic_bundle LP_MP-horizon-tracking DD_ILP CONIC_BUNDLE arboricity lingeling m stdc++ pthread )

add_executable(horizon_tracking_subgradient horizon_tracking_subgradient.cpp)
target_link_libraries(horizon_tracking_subgradient LP_MP-horizon-tracking DD_ILP arboricity lingeling m stdc++ pthread )
