#!/usr/bin/env bash
# Crea el entorno de R con xcms + ecosistema RforMassSpectrometry.
# No requiere root: todo se instala en el home del usuario.
#
#   bash env/setup.sh
#
# Después de esto:
#   eval "$(~/bin/micromamba shell hook -s bash)"
#   micromamba activate dinora
#   Rscript R/01_preprocess.R
set -euo pipefail

MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/share/micromamba}"
MAMBA_BIN="$HOME/bin/micromamba"
ENV_NAME="dinora"

if [ ! -x "$MAMBA_BIN" ]; then
    echo ">> Descargando micromamba..."
    mkdir -p "$HOME/bin"
    tmp="$(mktemp -d)"
    curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest -o "$tmp/mm.tar.bz2"
    tar -xjf "$tmp/mm.tar.bz2" -C "$tmp" bin/micromamba
    mv "$tmp/bin/micromamba" "$MAMBA_BIN"
    rm -rf "$tmp"
fi

export MAMBA_ROOT_PREFIX

echo ">> Creando entorno '$ENV_NAME' (esto descarga ~2-3 GB)..."
"$MAMBA_BIN" create -y -n "$ENV_NAME" -f env/environment.yml

echo
echo ">> Listo. Para activarlo:"
echo "   eval \"\$($MAMBA_BIN shell hook -s bash)\""
echo "   micromamba activate $ENV_NAME"
echo
echo ">> Comprobación rápida:"
echo "   micromamba run -n $ENV_NAME Rscript -e 'library(xcms); packageVersion(\"xcms\")'"
