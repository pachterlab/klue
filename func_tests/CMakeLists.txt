project(FuncTests)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtests.sh
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_1.fq.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_2.fq.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_w.fq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_x.fq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_y.fq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/snp_x.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/snp_y.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/snp_z.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/snp_insertion_x.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/snp_insertion_y.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/extendedvariation_x.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/extendedvariation_y.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1_V.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1revcomp_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2_V.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2revcomp_V.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/expected_output)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/extendedvariation_expected_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/extendedvariation_expected_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/extendedvariation_expected_V0.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/extendedvariation_expected_V1.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_insertion_expected_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_insertion_expected_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_insertion_expected_V0.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_insertion_expected_V1.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xy_expected_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xy_expected_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xy_expected_V0.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xy_expected_V1.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xyz_expected_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xyz_expected_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xyz_expected_V0.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xyz_expected_V1.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/klue_output/snp_xyz_expected_V2.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/klue_output)