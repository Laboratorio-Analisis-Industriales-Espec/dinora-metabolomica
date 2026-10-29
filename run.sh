#!/usr/bin/env bash
# Ejecuta el pipeline completo.
#   ./run.sh                     usa config/config.yml
#   ./run.sh config/otro.yml     usa otra configuración
#   MAMBA_BIN=/ruta/micromamba ./run.sh   si micromamba no está en el PATH
set -euo pipefail
cd "$(dirname "$0")"

CONFIG="${1:-config/config.yml}"

MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/share/micromamba}"
export MAMBA_ROOT_PREFIX

MAMBA_BIN="${MAMBA_BIN:-}"
if [ -z "$MAMBA_BIN" ]; then
    for cand in "$HOME/bin/micromamba" "$(command -v micromamba 2>/dev/null || true)" \
                "/tmp/opencode/bin/micromamba"; do
        [ -n "$cand" ] && [ -x "$cand" ] && { MAMBA_BIN="$cand"; break; }
    done
fi
[ -n "$MAMBA_BIN" ] || { echo "No encuentro micromamba. Corre env/setup.sh"; exit 1; }

run() { "$MAMBA_BIN" run -n dinora Rscript "$@"; }

echo ">> Config: $CONFIG"
echo ">> 1/2 Preprocesamiento"
run R/01_preprocess.R "$CONFIG"

echo ">> 2/2 Anotación"
run R/02_annotate.R "$CONFIG"

echo ">> Hecho. Resultados en results/"
