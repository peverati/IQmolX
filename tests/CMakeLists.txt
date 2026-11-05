find_package(Qt5 COMPONENTS Core REQUIRED)

add_executable(XtbProcessTest EXCLUDE_FROM_ALL
   XtbProcessTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/GxtbSetup.C
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
   target_link_libraries(XtbLibraryTest PRIVATE ${XTB_LIBRARY})
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
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/GxtbSetup.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbLibrary.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
      ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbOutput.C
   )
   target_include_directories(XtbBackendTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src ${XTB_ROOT}/include)
   target_link_libraries(XtbBackendTest PRIVATE Qt5::Core ${XTB_LIBRARY})
   set_target_properties(XtbBackendTest PROPERTIES BUILD_RPATH ${XTB_ROOT}/lib)
   if (OpenMP_CXX_FOUND)
      target_link_libraries(XtbBackendTest PRIVATE OpenMP::OpenMP_CXX)
   endif()
   add_test(NAME XtbBackendTest COMMAND XtbBackendTest ${XTB_ROOT}/bin/xtb${CMAKE_EXECUTABLE_SUFFIX} ${CMAKE_CURRENT_SOURCE_DIR}/data)
endif()


# Python engine helper: the C++ client, tested against the demo engine
# (needs a Python with numpy; the test is given the interpreter and the helper's source folder)
find_package(Python3 COMPONENTS Interpreter)
add_executable(PythonEngineTest EXCLUDE_FROM_ALL
   PythonEngineTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/PythonEngine.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
)
target_include_directories(PythonEngineTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(PythonEngineTest PRIVATE Qt5::Core)
if (Python3_EXECUTABLE)
   add_test(NAME PythonEngineTest COMMAND PythonEngineTest ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/../python/iqmolx-engines/src)
endif()


# The g-xTB backend: the download (a local archive), the command line, and a real g-xTB if given
add_executable(GxtbTest EXCLUDE_FROM_ALL
   GxtbTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbProcess.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/GxtbSetup.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbOutput.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
)
target_include_directories(GxtbTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(GxtbTest PRIVATE Qt5::Core)
add_test(NAME GxtbTest COMMAND GxtbTest)


# The creation of the Python environment with uv, tested with a fake uv
add_executable(PythonSetupTest EXCLUDE_FROM_ALL
   PythonSetupTest.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/PythonSetup.C
)
target_include_directories(PythonSetupTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(PythonSetupTest PRIVATE Qt5::Core)
add_test(NAME PythonSetupTest COMMAND PythonSetupTest)


# The Python engines window, against the demo engine (offscreen Qt)
add_executable(PythonDialogTest EXCLUDE_FROM_ALL PythonDialogTest.C)
target_include_directories(PythonDialogTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(PythonDialogTest PRIVATE Util Qt5::Widgets Qt5::Gui Qt5::Core)
if (Python3_EXECUTABLE)
   add_test(NAME PythonDialogTest COMMAND PythonDialogTest ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/../python/iqmolx-engines/src)
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
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/GxtbSetup.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/XtbOutput.C
   ${CMAKE_CURRENT_SOURCE_DIR}/../src/Util/Vibrations.C
)
target_include_directories(XtbOrbitalsTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
target_link_libraries(XtbOrbitalsTest PRIVATE Parser Data Util Math yaml-cpp openbabel ${OPENMESH_LIBRARIES}
   Qt5::Widgets Qt5::Core Qt5::Xml Qt5::OpenGL Qt5::Gui Qt5::Sql Qt5::PrintSupport)
add_test(NAME XtbOrbitalsTest COMMAND XtbOrbitalsTest ${XTB_ROOT}/bin/xtb${CMAKE_EXECUTABLE_SUFFIX})
