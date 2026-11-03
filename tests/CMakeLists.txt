find_package(Qt5 COMPONENTS Core REQUIRED)

add_executable(XtbProcessTest EXCLUDE_FROM_ALL
   XtbProcessTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbOutput.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
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

   # Both backends side by side
   add_executable(XtbBackendTest EXCLUDE_FROM_ALL
      XtbBackendTest.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbLibrary.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbOutput.C
   )
   target_include_directories(XtbBackendTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src ${XTB_ROOT}/include)
   target_link_libraries(XtbBackendTest PRIVATE Qt5::Core ${XTB_ROOT}/lib/libxtb.dylib)
   set_target_properties(XtbBackendTest PROPERTIES BUILD_RPATH ${XTB_ROOT}/lib)
   if (OpenMP_CXX_FOUND)
      target_link_libraries(XtbBackendTest PRIVATE OpenMP::OpenMP_CXX)
   endif()
   add_test(NAME XtbBackendTest COMMAND XtbBackendTest ${XTB_ROOT}/bin/xtb ${CMAKE_CURRENT_SOURCE_DIR}/data)
endif()


# Optimizer: no dependencies
add_executable(OptimizerTest EXCLUDE_FROM_ALL OptimizerTest.C)
target_include_directories(OptimizerTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
add_test(NAME OptimizerTest COMMAND OptimizerTest)


# The xtb window (needs the Util library, runs on Qt's offscreen platform)
find_package(Qt5 COMPONENTS Widgets REQUIRED)
add_executable(XtbDialogTest EXCLUDE_FROM_ALL XtbDialogTest.C)
target_include_directories(XtbDialogTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(XtbDialogTest PRIVATE Util Qt5::Widgets Qt5::Gui Qt5::Core)
add_test(NAME XtbDialogTest COMMAND XtbDialogTest)


# Vibrational analysis: no dependencies
add_executable(VibrationsTest EXCLUDE_FROM_ALL
   VibrationsTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
)
target_include_directories(VibrationsTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
add_test(NAME VibrationsTest COMMAND VibrationsTest)


# Orbitals from xtb, checked on a grid with IQmol's own basis function code
find_package(Qt5 COMPONENTS Xml OpenGL Gui REQUIRED)
add_executable(XtbOrbitalsTest EXCLUDE_FROM_ALL
   XtbOrbitalsTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbOutput.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
)
target_include_directories(XtbOrbitalsTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(XtbOrbitalsTest PRIVATE Parser Data Util Math yaml-cpp openbabel ${OPENMESH_LIBRARIES}
   Qt5::Widgets Qt5::Core Qt5::Xml Qt5::OpenGL Qt5::Gui Qt5::Sql Qt5::PrintSupport)
add_test(NAME XtbOrbitalsTest COMMAND XtbOrbitalsTest ${XTB_ROOT}/bin/xtb)
