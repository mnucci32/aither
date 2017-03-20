function InstallMPI($downloads) {
    md mpi
    cd mpi
    # install MPI SDK and Runtime
    Write-Host "Installing Microsoft MPI SDK..."
    appveyor DownloadFile http://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/msmpisdk.msi -FileName $downloads
    Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i $downloads\msmpisdk.msi" -Wait
    Write-Host "Microsoft MPI SDK installation complete"
    Write-Host "Installing Microsoft MPI Runtime..."
    appveyor DownloadFile http://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/MSMpiSetup.exe -FileName $downloads
    Start-Process -FilePath $downloads\MSMpiSetup.exe -ArgumentList -unattend -Wait
    # Get MPI environment variables
    $envfile = "MPI_Env.cmd"
    Write-Host "Saving Microsoft MPI environment variables to" $envfile
    $envlist = @("MSMPI_BIN", "MSMPI_INC", "MSMPI_LIB32", "MSMPI_LIB64")
    $stream = [IO.StreamWriter] $envfile
    foreach ($variable in $envlist) {
        $value = [Environment]::GetEnvironmentVariable($variable, "Machine")
        if ($value) { $stream.WriteLine("SET $variable=$value") }
        if ($value) { Write-Host "$variable=$value" }
    }
    $stream.Close()
    cd ..
}

function main() {
    $downloads = "C:\Downloads\MSMPI"
    InstallMPI($downloads)
}

main