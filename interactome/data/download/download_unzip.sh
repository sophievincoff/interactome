#!/bin/bash

URL="$1"
DEST_DIR="$2"
FILENAME="$3"        # e.g., intact.zip or biogrid.tab3.gz
DELETE_ZIP="$4"      # "true" or "false"

FILE="$DEST_DIR/$FILENAME"

mkdir -p "$DEST_DIR"

echo "Starting download from $URL..."
wget "$URL" -O "$FILE"

# Handle .zip files
if [[ "$FILENAME" == *.zip ]]; then
  echo "File is a .zip archive. Unzipping to $DEST_DIR in 10s..."
  sleep 10
  unzip "$FILE" -d "$DEST_DIR"

  if [[ "$DELETE_ZIP" == "true" ]]; then
    echo "delete_zip=true: removing $FILE in 30s"
    sleep 30
    rm -f "$FILE"
  fi

# Handle .gz files
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

echo "Done."
