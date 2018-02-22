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
    Write-Host "Installing CMake 3.4.0 ..." -ForegroundColor Cyan
    $exePath = "$($env:USERPROFILE)\cmake-3.4.0-rc2-win32-x86.exe"
    Write-Host "Downloading..."
    (New-Object Net.WebClient).DownloadFile('https://cmake.org/files/v3.4/cmake-3.4.0-rc2-win32-x86.exe', $exePath)
    Write-Host "Installing..."
    cmd /c start /wait $exePath /S
    cmake --version
    Write-Host "CMake 3.4.0 installed" -ForegroundColor Green
}

function main() {
    InstallMPI
    InstallCmake
}

main