#!/bin/bash

URL="$1"
DEST_DIR="$2"
FILENAME="$3"  # optional, e.g., intact.zip or biogrid.tab3.zip

FILE="$DEST_DIR/$FILENAME"

mkdir -p "$DEST_DIR"

echo "Starting download from $URL..."
wget "$URL" -O "$FILE"

echo "Waiting 30s..."
sleep 30

echo "Download complete. Unzipping to $DEST_DIR..."
unzip "$FILE" -d "$DEST_DIR"

echo "Done."
