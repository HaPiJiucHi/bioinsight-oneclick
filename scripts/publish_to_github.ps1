param(
  [string]$RepoName = "differential-analysis-software",
  [string]$Visibility = "public",
  [string]$Version = "v1.0.0"
)

$ErrorActionPreference = "Stop"

$gh = Get-Command gh -ErrorAction SilentlyContinue
if (-not $gh) {
  $localGh = Join-Path (Split-Path -Parent $PSScriptRoot) "..\tools\gh\gh.exe"
  if (Test-Path $localGh) {
    $gh = Get-Item $localGh
  }
}

if (-not $gh) {
  throw "找不到 GitHub CLI。请先安装 gh，或使用本机 pipeline/tools/gh/gh.exe。"
}

& $gh.Source auth status *> $null
if ($LASTEXITCODE -ne 0) {
  & $gh.Source auth login --hostname github.com --web --git-protocol https --scopes repo
}

$login = (& $gh.Source api user --jq ".login").Trim()
git config user.name $login
git config user.email "$login@users.noreply.github.com"

git remote remove origin 2>$null
& $gh.Source repo create "$login/$RepoName" "--$Visibility" --source . --remote origin --description "Windows desktop app for differential expression analysis with Shiny and limma." --push

& $gh.Source release create $Version "dist/差异分析软件.zip" --title "$Version" --notes-file RELEASE_NOTES.md

Write-Host "Published: https://github.com/$login/$RepoName"
