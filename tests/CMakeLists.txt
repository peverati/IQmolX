find_package(Qt5 COMPONENTS Core REQUIRED)

add_executable(XtbProcessTest EXCLUDE_FROM_ALL
   XtbProcessTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
)
target_include_directories(XtbProcessTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(XtbProcessTest PRIVATE Qt5::Core)

add_test(NAME XtbProcessTest COMMAND XtbProcessTest)


# libxtb tests: only built when XTB_ROOT points to an xtb installation (see build-xtb-macos.sh)
if (XTB_ROOT AND EXISTS ${XTB_ROOT}/include/xtb.h)
   add_executable(XtbLibraryTest EXCLUDE_FROM_ALL XtbLibraryTest.C)
   target_include_directories(XtbLibraryTest PRIVATE ${XTB_ROOT}/include)
   target_link_libraries(XtbLibraryTest PRIVATE ${XTB_ROOT}/lib/libxtb.dylib)
   set_target_properties(XtbLibraryTest PROPERTIES BUILD_RPATH ${XTB_ROOT}/lib)
   if (OpenMP_CXX_FOUND)
      # Same OpenMP runtime as IQmol itself, to test that it coexists with libxtb's libgomp
      target_link_libraries(XtbLibraryTest PRIVATE OpenMP::OpenMP_CXX)
   endif()
   add_test(NAME XtbLibraryTest COMMAND XtbLibraryTest)
endif()
