cmake_minimum_required(VERSION 2.8)

set(SRCS
    Edge.cc
    FloatImage.cc
    Gaussian.cc
    GLine2D.cc
    GLineSegment2D.cc
    GrayModel.cc
    Homography33.cc
    MathUtil.cc
    Quad.cc
    Segment.cc
    TagDetection.cc
    TagDetector.cc
    TagFamily.cc
    UnionFindSimple.cc)

add_library(apriltags SHARED ${SRCS})

file(GLOB HEADERS "${PROJECT_SOURCE_DIR}/include/apriltags/*.h"
    "${PROJECT_SOURCE_DIR}/include/apriltags/*.hpp")

set_target_properties(apriltags PROPERTIES
    OUTPUT_NAME apriltags
    PUBLIC_HEADER "${HEADERS}")

install(TARGETS apriltags
    EXPORT AprilTagsTargets
    LIBRARY DESTINATION "${LIBRARY_INSTALL_DIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}"
    COMPONENT dev)

target_link_libraries(apriltags ${AprilTags_LIBS})
