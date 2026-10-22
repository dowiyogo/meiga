#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${repo_root}/build"
app_dir="${build_dir}/Applications/G4HodoscopeSimulator"
log_file="${repo_root}/diagnostics/hodoscope_geometry_check.log"
output_file="${app_dir}/output.json"

cd "${repo_root}"
cmake -S src -B "${build_dir}"
cmake --build "${build_dir}" --target G4HodoscopeSimulator -j2
cp src/Applications/G4HodoscopeSimulator/G4HodoscopeSimulator.json "${app_dir}/G4HodoscopeSimulator.json"
cp src/Applications/G4HodoscopeSimulator/DetectorList.xml "${app_dir}/DetectorList.xml"

mkdir -p diagnostics
(
	cd "${app_dir}"
	./G4HodoscopeSimulator -c ./G4HodoscopeSimulator.json
) >"${log_file}" 2>&1

check_overlaps_result="PASS: no severe overlap message found"
if grep -Eiq "overlap|GeomVol1002|GeomVol1003" "${log_file}"; then
	if grep -Eiq "GeomVol1002|GeomVol1003|overlap.*detected" "${log_file}"; then
		check_overlaps_result="CHECK LOG: possible overlap warning in diagnostics/hodoscope_geometry_check.log"
	fi
fi

python3 diagnostics/analyze_hodoscope_output.py "${output_file}" \
	--outdir diagnostics \
	--check-overlaps-result "${check_overlaps_result}"

echo "Geometry check log: ${log_file}"
echo "Output JSON: ${output_file}"
echo "CheckOverlaps: ${check_overlaps_result}"
