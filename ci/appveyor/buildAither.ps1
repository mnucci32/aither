function BuildAither() {
    # Go to build directory and build
    # Build Aither with cmake
    md build
    cd build
    cmake -G "Visual Studio 14 2015 Win64" ..
    cmake --build . --target INSTALL --config release
    cd ..
}

function main() {
    BuildAither
}

main