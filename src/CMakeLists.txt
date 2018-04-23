add_executable(brownie graphaln.cpp nodechain.cpp kmercounttable.cpp kmeroverlaptable.cpp readcorrection.cpp alignment.cpp  library.cpp kmernpp.cpp kmertable.cpp  dsnode.cpp nucleotide.cpp nodeendstable.cpp settings.cpp util.cpp tstring.cpp kmeroverlap.cpp graph.cpp brownie.cpp suffix_tree.c   markovChain.cpp )

add_executable(compare2 comparison.cpp alignment.cpp util.cpp library.cpp )
target_link_libraries(compare2 readfile pthread)
target_link_libraries(brownie readfile essaMEM pthread)

add_executable(sam2Alignment samToAlignment.cpp alignment.cpp util.cpp library.cpp  )
target_link_libraries(sam2Alignment readfile pthread)

if (ZLIB_FOUND)
   target_link_libraries(brownie ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)


if (ZLIB_FOUND)
   target_link_libraries(compare2 ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)


if (ZLIB_FOUND)
   target_link_libraries(sam2Alignment ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

install(TARGETS brownie RUNTIME DESTINATION bin)
add_subdirectory(readfile)
add_subdirectory(essaMEM-master)
