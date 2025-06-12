#!/bin/bash

URL="https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.246/BIOGRID-ALL-4.4.246.tab3.zip"
DEST_DIR="$(git rev-parse --show-toplevel)/interactome/data_files/raw/biogrid"
FILE="$DEST_DIR/BIOGRID-ALL-4.4.246.tab3.zip"

mkdir -p "$DEST_DIR"  # create directory if it doesn't exist

echo "Starting download..."
wget "$URL" -O "$FILE"

echo "Waiting 30s"
sleep 30

echo "Download complete. Unzipping..."
unzip "$FILE" -d "$DEST_DIR"

echo "Done."
