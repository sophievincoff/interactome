#!/usr/bin/env bash
set -euo pipefail

# Usage: ./download_unzip.sh URL DEST_DIR FILENAME [UNZIP] [DELETE_ZIP] [INSECURE]
# Example:
#   ./download_unzip.sh "https://..." "/tmp/negatome" "file.txt" false false

URL="${1:?need URL}"
DEST_DIR="${2:?need DEST_DIR}"
FILENAME="${3:?need FILENAME (e.g., file.txt or file.gz)}"
UNZIP="${4:-true}"              # "true" or "false"
DELETE_ZIP="${5:-false}"        # "true" or "false"
INSECURE="${6:-false}"          # "true" to allow insecure SSL (skip verification)

FILE="$DEST_DIR/$FILENAME"
mkdir -p "$DEST_DIR"

echo "Starting download from $URL..."

# Try to find a CA bundle via certifi (works on most conda/py installs)
CA_FILE=""
if command -v python >/dev/null 2>&1; then
  CA_FILE="$(python - <<'PY' || true
import sys
try:
    import certifi
    print(certifi.where())
except Exception:
    sys.exit(1)
PY
  )"
fi

# Build wget command
wget_cmd=(wget -O "$FILE" --timeout=60 --tries=3 --no-verbose)

if [[ "$INSECURE" == "true" ]]; then
  wget_cmd+=(--no-check-certificate)
elif [[ -n "$CA_FILE" && -f "$CA_FILE" ]]; then
  # Use a trusted CA bundle if we found one
  wget_cmd+=(--ca-certificate="$CA_FILE")
else
  # Last resort: use the system default CA store (usually works if installed)
  # You could also set --ca-directory=/etc/ssl/certs if needed.
  :
fi

wget_cmd+=("$URL")

# Run download and capture errors
if ! "${wget_cmd[@]}"; then
  echo "ERROR: Download failed for $URL" >&2
  echo "Tip: If this is a trusted link and the CA store is missing, retry with INSECURE=true:" >&2
  echo "     $0 \"$URL\" \"$DEST_DIR\" \"$FILENAME\" \"$UNZIP\" \"$DELETE_ZIP\" true" >&2
  exit 1
fi

if [[ "$UNZIP" == "true" ]]; then
  # Post-download handling (same logic you had)
  if [[ "$FILENAME" == *.zip ]]; then
    echo "File is a .zip archive. Unzipping to $DEST_DIR in 10s..."
    sleep 10
    unzip -o "$FILE" -d "$DEST_DIR"

    if [[ "$DELETE_ZIP" == "true" ]]; then
      echo "delete_zip=true: removing $FILE in 30s"
      sleep 30
      rm -f "$FILE"
    fi

  elif [[ "$FILENAME" == *.gz ]]; then
    echo "File is a .gz archive. Extracting in 10s..."
    sleep 10
    gunzip -c "$FILE" > "${FILE%.gz}"

    if [[ "$DELETE_ZIP" == "true" ]]; then
      echo "delete_zip=true: removing $FILE in 30s"
      sleep 30
      rm -f "$FILE"
    fi

  else
    echo "File is not a .zip or .gz archive. Skipping extraction and deletion."
  fi
fi


echo "Done."
