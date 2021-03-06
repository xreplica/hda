include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/core" "${PROJECT_SOURCE_DIR}/src/algorithms" "${CMAKE_SOURCE_DIR}/algorithms")

#add_executable(clustering clustering.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(efc efc.cpp $<TARGET_OBJECTS:fuzzyentropy> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(fcm fcm.cpp $<TARGET_OBJECTS:fuzzycmeans> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(fcmp fcmp.cpp $<TARGET_OBJECTS:fuzzycmeans> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(hda hda.cpp $<TARGET_OBJECTS:hdiana> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance> $<TARGET_OBJECTS:params>)
#add_executable(km km.cpp $<TARGET_OBJECTS:kmeans> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(kmp kmp.cpp $<TARGET_OBJECTS:kmeansplus> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(rac rac.cpp $<TARGET_OBJECTS:racamaic> $<TARGET_OBJECTS:fuzzycmeans> $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
#add_executable(validate validate.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)

add_executable(clustering clustering.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(efc efc.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(fcm fcm.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(fcmp fcmp.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(hda hda.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance> $<TARGET_OBJECTS:params>)
add_executable(km km.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(kmp kmp.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(rac rac.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)
add_executable(validate validate.cpp $<TARGET_OBJECTS:point> $<TARGET_OBJECTS:cluster> $<TARGET_OBJECTS:matrix> $<TARGET_OBJECTS:validation> $<TARGET_OBJECTS:distance>)

target_link_libraries(efc fuzzyentropy_static)
target_link_libraries(fcm fuzzycmeans_static)
target_link_libraries(fcmp fuzzycmeans_static)
target_link_libraries(hda hdiana_static)
target_link_libraries(km kmeans_static)
target_link_libraries(kmp kmeansplus_static)
target_link_libraries(rac racamaic_static fuzzycmeans_static)
