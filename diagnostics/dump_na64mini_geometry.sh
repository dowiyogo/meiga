#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 G4HodoscopeSimulator OUTPUT.csv" >&2
  exit 2
fi

binary=$(realpath "$1")
dump=$(realpath -m "$2")
mkdir -p "$(dirname "$dump")"
run_dir=$(mktemp -d "${TMPDIR:-/tmp}/meiga-geom-dump.XXXXXX")
trap 'rm -rf "$run_dir"' EXIT
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

python3 - "$here" "$run_dir" "$dump" <<'PY'
import json
import sys
from pathlib import Path

here, run_dir, dump = map(Path, sys.argv[1:])
app = here.parent / 'src' / 'Applications' / 'G4HodoscopeSimulator'
source = app / 'G4HodoscopeSimulator.json'
config = json.loads(source.read_text())
config['Input']['InputNParticles'] = 1
config['Output']['OutputFile'] = str(Path(str(dump) + '.event.json'))
config['DetectorList'] = str(app / 'DetectorList.xml')
config['DetectorProperties'] = str(here.parent / 'src' / 'Config' / 'DetectorProperties.xml')
target = run_dir / 'geometry_dump.json'
target.write_text(json.dumps(config, indent=2) + '\n')
PY

GEOM_DUMP="$dump" "$binary" -c "$run_dir/geometry_dump.json"
test -s "$dump"
