cmake_minimum_required(VERSION 3.10)

# Set policies first to avoid warnings
cmake_policy(SET CMP0156 NEW)  # Use modern library de-duplication

project(QteExternals C)

# Fix Xcode SDK warning by explicitly setting the SDK path
set(CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk")

# Build for both Intel (x86_64) and Apple Silicon (arm64)
set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")

# Set the build directory explicitly
set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build")

# Force all paths to be within the project
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(MAX_INSTALL_DIR "${CMAKE_BINARY_DIR}")
set(MAX_BUILD_DIR "${CMAKE_BINARY_DIR}")

# Override environment variables that might be used by scripts
set(ENV{MAX_BUILD_DIR} "${CMAKE_BINARY_DIR}")
set(ENV{MAX_INSTALL_DIR} "${CMAKE_BINARY_DIR}")

# Override any default paths that might be causing issues
set(ZERO_CHECK_SOURCE_DIR "${CMAKE_BINARY_DIR}")
set(ZERO_CHECK_BINARY_DIR "${CMAKE_BINARY_DIR}")
set(CMAKE_CACHEFILE_DIR "${CMAKE_BINARY_DIR}")

# Include the max-sdk pretarget script using your absolute path.
# TEMPORARILY COMMENT THIS OUT TO TEST if it's the source of the problem
# include("/Users/pedrogonzalez/Documents/Max9/Packages/max-sdk/source/max-sdk-base/script/max-pretarget.cmake")

# Instead, define the minimum needed from the max-pretarget.cmake
set(MAX_SDK_INCLUDES "/Users/pedrogonzalez/Documents/Max9/Packages/max-sdk/source/c74support/max-includes")
set(MAX_SDK_MSP_INCLUDES "/Users/pedrogonzalez/Documents/Max9/Packages/max-sdk/source/c74support/msp-includes")
set(MAX_SDK_JIT_INCLUDES "/Users/pedrogonzalez/Documents/Max9/Packages/max-sdk/source/c74support/jit-includes")

#############################################################
# MAX EXTERNALS
#############################################################

# Include Max SDK header directories.
include_directories(
    "${MAX_SDK_INCLUDES}"
    "${MAX_SDK_MSP_INCLUDES}"
    "${MAX_SDK_JIT_INCLUDES}"
)

# Set the C standard to C99 for <complex.h> support.
set(CMAKE_C_STANDARD 99)

# Path to your minimal Info.plist file.
set(MACOSX_BUNDLE_INFO_PLIST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Info.plist")

# Create a custom function to reduce repetition
function(add_max_external target_name source_file)
    add_library(${target_name} MODULE ${source_file})
    set_target_properties(${target_name} PROPERTIES
        BUNDLE True
        BUNDLE_EXTENSION "mxo"
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        MACOSX_BUNDLE_INFO_PLIST "${MACOSX_BUNDLE_INFO_PLIST_FILE}"
    )
endfunction()

# Quantum Harmonic Oscillator (qte.quantumho)
add_max_external(qte.quantumho quantum_ho.c)

# Hermitian Matrix Maker (qte.hermitmaker)
add_max_external(qte.hermitmaker herm_maker.c)

# Random Hermitian Matrix Generator (qte.randherm)
add_max_external(qte.randherm rand_herm.c)

# Hermitian Matrix Combiner (qte.hermcombiner)
add_max_external(qte.hermcombiner herm_combiner.c)

# Eigenbasis Calculator (qte.eigencalc) - requires Accelerate framework
add_max_external(qte.eigencalc eigen_calc.c)
target_link_libraries(qte.eigencalc PUBLIC "-framework Accelerate")

# InitialState_CoefficientCalculator (qte.initstatecalc)
add_max_external(qte.initstatecalc initstate_calc.c)

# Time Developer (qte.timedev)
add_max_external(qte.timedev time_dev.c)

# Add this to ensure we're not trying to use /Users/externals
set_directory_properties(PROPERTIES
    ADDITIONAL_CLEAN_FILES ""
    ADDITIONAL_MAKE_CLEAN_FILES ""
)

# Explicitly disable the ZERO_CHECK target location to avoid permission issues
set(CMAKE_SUPPRESS_REGENERATION TRUE)