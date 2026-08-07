param(
  [string]$InputRoot = "outputs/tcga_gtex_integration/recount3_objects",
  [string]$Rscript = "Rscript",
  [string]$Projects = "",
  [string]$K = "1,2,3",
  [int]$HvgN = 5000,
  [switch]$Force
)

$ErrorActionPreference = "Stop"

$normalizer = Join-Path $PSScriptRoot "normalize_tcga_gtex_ruvprps_unsupervised.R"
if (!(Test-Path -LiteralPath $normalizer)) {
  throw "Normalizer script not found: $normalizer"
}
if (!(Test-Path -LiteralPath $InputRoot)) {
  throw "Input root not found: $InputRoot"
}

$selected = @()
if ($Projects.Trim().Length -gt 0) {
  $selected = $Projects.Split(",") | ForEach-Object { $_.Trim().ToUpperInvariant() } | Where-Object { $_ }
}

$projectDirs = Get-ChildItem -LiteralPath $InputRoot -Directory | Sort-Object Name
foreach ($projectDir in $projectDirs) {
  $project = $projectDir.Name.ToUpperInvariant()
  if ($selected.Count -gt 0 -and $project -notin $selected) {
    continue
  }

  $inputRds = Join-Path $projectDir.FullName "tcga_gtex_integrated_se.rds"
  $outputRds = Join-Path $projectDir.FullName ("{0}_all_samples_ruvprps_unsupervised_cca.rds" -f $project)

  if (!(Test-Path -LiteralPath $inputRds)) {
    Write-Host "SKIP ${project}: tcga_gtex_integrated_se.rds not found"
    continue
  }
  if ((Test-Path -LiteralPath $outputRds) -and !$Force) {
    Write-Host "EXISTS ${project}: $outputRds"
    continue
  }

  Write-Host "START $project"
  & $Rscript $normalizer `
    --input-rds $inputRds `
    --output-rds $outputRds `
    --samples-for-prps all `
    --approach cca `
    --hvg-n $HvgN `
    --k $K

  if ($LASTEXITCODE -ne 0) {
    throw "PRPS/RUV failed for $project with exit code $LASTEXITCODE"
  }
  Write-Host "DONE $project"
}
