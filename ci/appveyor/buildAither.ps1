function BuildAither() {
    # Go to build directory and build
    # Build Aither with cmake
    md build
    cd build
    cmake -G "Visual Studio 14 2015 Win64" -DMPI_DIR="C:\Program Files (x86)\Microsoft SDKs\MPI" ..
    cmake --build . --target INSTALL --config release
    cd ..
}

function main() {
    BuildAither
}

main