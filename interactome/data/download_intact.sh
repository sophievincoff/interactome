#!/bin/bash

URL="ftp://ftp.ebi.ac.uk/pub/databases/intact/current/all.zip"
DEST_DIR="$(git rev-parse --show-toplevel)/interactome/data_files/raw/intact"
FILE="$DEST_DIR/all_release250.zip"

mkdir -p "$DEST_DIR"  # create directory if it doesn't exist

echo "Starting download..."
wget "$URL" -O "$FILE"

echo "Waiting 30s"
sleep 30

echo "Download complete. Unzipping..."
unzip "$FILE" -d "$DEST_DIR"

echo "Done."
