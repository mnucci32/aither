function InstallMPI() {
    md mpi
    cd mpi
    # install MPI SDK and Runtime
    Write-Host "Installing Microsoft MPI SDK..."
    appveyor DownloadFile http://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/msmpisdk.msi
    Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i msmpisdk.msi" -Wait
    Write-Host "Microsoft MPI SDK installation complete"
    Write-Host "Installing Microsoft MPI Runtime..."
    appveyor DownloadFile http://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/MSMpiSetup.exe
    Start-Process -FilePath MSMpiSetup.exe -ArgumentList -unattend -Wait
    Write-Host "Microsoft MPI Runtime installation complete..."
    cd ..
}

function InstallCmake() {
    md cmake-38
    cd cmake-38
    # install cmake
    Write-Host "Installing CMake 3.8.0 ..."
    appveyor DownloadFile https://cmake.org/files/v3.8/cmake-3.8.0-win64-x64.msi
    Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i cmake-3.8.0-win64-x64.msi" -Wait
    Write-Host "CMake 3.8.0 installed..."
    cd ..
}

function main() {
    InstallMPI
    # InstallCmake
}

main