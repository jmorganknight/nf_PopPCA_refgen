#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  download_models.sh [--outdir <dir>] [--bundle-url <url>]

Options:
  --outdir       Destination directory (default: example/models)
  --bundle-url   URL to model bundle .tar.gz (default: GitHub Releases v1.0.0 asset URL)
  -h, --help     Show this help message

Environment:
  MODEL_BUNDLE_URL  Optional default URL override
EOF
}

OUTDIR="example/models"
DEFAULT_URL="https://github.com/jmorganknight/nf_PopPCA_refgen/releases/download/v1.0.0/nf_PopPCA_refgen_v1.0.0_models.tar.gz"
BUNDLE_URL="${MODEL_BUNDLE_URL:-$DEFAULT_URL}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;
    --bundle-url)
      BUNDLE_URL="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
  echo "ERROR: curl or wget is required." >&2
  exit 2
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
ARCHIVE="$TMP_DIR/model_bundle.tar.gz"
EXTRACT_DIR="$TMP_DIR/extracted"
mkdir -p "$EXTRACT_DIR" "$OUTDIR"

fetch() {
  local url="$1"
  local dest="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 3 --retry-delay 2 -o "$dest" "$url"
  else
    wget -O "$dest" "$url"
  fi
}

echo "Downloading model bundle from: $BUNDLE_URL"
fetch "$BUNDLE_URL" "$ARCHIVE"

echo "Extracting model bundle into: $EXTRACT_DIR"
tar -xf "$ARCHIVE" -C "$EXTRACT_DIR"

# Locate the models/ folder inside extracted archive
SRC_MODELS_DIR="$(find "$EXTRACT_DIR" -type d -name "models" | head -n 1 || true)"
if [[ -z "$SRC_MODELS_DIR" ]]; then
  SRC_MODELS_DIR="$EXTRACT_DIR"
fi

echo "Staging files to destination: $OUTDIR"
cp -r "$SRC_MODELS_DIR/"* "$OUTDIR/"

echo "Model bundle successfully unpacked and staged into: $OUTDIR"