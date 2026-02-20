#!/usr/bin/env bash
shopt -s globstar

# set absolute path !
files=(~/BicWin26/terrain/**/Mat_portatif/Trisonica/*.txt)
# destination folder
DEST="data_trisonica_portable"

mkdir -p "$DEST"

for f in "${files[@]}"; do
  echo $f
  [[ -f "$f" ]] && cp -- "$f" "$DEST" && echo "Copied: $f"
done

echo "Done."
