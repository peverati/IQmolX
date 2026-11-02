find_package(Qt5 COMPONENTS Core REQUIRED)

add_executable(XtbProcessTest EXCLUDE_FROM_ALL
   XtbProcessTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
)
target_include_directories(XtbProcessTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(XtbProcessTest PRIVATE Qt5::Core)

add_test(NAME XtbProcessTest COMMAND XtbProcessTest)
