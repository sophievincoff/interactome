#!/bin/bash

URL="$1"
DEST_DIR="$2"
FILENAME="$3"  # optional, e.g., intact.zip or biogrid.tab3.zip

FILE="$DEST_DIR/$FILENAME"

mkdir -p "$DEST_DIR"

echo "Starting download from $URL..."
wget "$URL" -O "$FILE"

# Only unzip if the file ends in .zip
if [[ "$FILENAME" == *.zip ]]; then
  echo "Download complete. Unzipping to $DEST_DIR in 10s..."
  sleep 10
  unzip "$FILE" -d "$DEST_DIR"
else
  echo "Download complete. File is not a zip archive, skipping unzip."
fi

echo "Done."
