include_directories(gtest/include ../src)
add_executable(unittest mmtest.cpp utiltest.cpp alignmenttest.cpp scaffoldtest.cpp readfiletest.cpp
        nucleotidetest.cpp kmermdtest.cpp kmertest.cpp tstringtest.cpp
        ../src/tstring.cpp ../src/nucleotide.cpp ../src/kmeroverlap.cpp ../src/alignment.cpp
        ../src/util.cpp)

target_link_libraries(unittest readfile gtest essaMEM
                      gtest_main ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} pthread)

add_subdirectory(gtest)
