param(
    [int]$NP = 1,          # MPI processes
    [int]$OMP = 1,         # OpenMP threads per process
    [string]$Dest = "L:\coding projects\AMCFD project\AMCFD-FORTRAN new\AMCFD\timing_results"
)

Write-Host "Running: MPI np=$NP, OMP threads=$OMP  (total cores used: $($NP * $OMP))"

wsl -- bash -c "cd '/mnt/l/coding projects/AMCFD project/AMCFD-FORTRAN new/AMCFD/fortran_new' && export OMP_NUM_THREADS=$OMP && mpirun --allow-run-as-root -np $NP ./cluster_main"

New-Item -ItemType Directory -Force -Path $Dest | Out-Null
$n = (Get-ChildItem "$Dest\timing_report_*.txt" -ErrorAction SilentlyContinue).Count + 1
$tag = "np${NP}_omp${OMP}"
$outName = "timing_report_${n}_${tag}.txt"
Copy-Item "result\timing_report.txt" "$Dest\$outName"
Write-Host "Saved as $outName"
