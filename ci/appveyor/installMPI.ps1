md mpi
cd mpi
# install MPI SDK and Runtime
Remove-item alias:curl
Write-Host "Installing Microsoft MPI SDK..."
curl -fsS -o msmpisdk.msi http://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/msmpisdk.msi
Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i msmpisdk.msi" -Wait
Write-Host "Microsoft MPI SDK installation complete"
Write-Host "Installing Microsoft MPI Runtime..."
curl -fsS -o MSMpiSetup.exe http://download.microsoft.com/download/B/2/E/B2EB83FE-98C2-4156-834A-E1711E6884FB/MSMpiSetup.exe
Start-Process -FilePath MSMpiSetup.exe -ArgumentList -unattend -Wait
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
Start-Process -FilePath $envfile
cd ..