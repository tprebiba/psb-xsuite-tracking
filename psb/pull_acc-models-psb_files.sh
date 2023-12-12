#!/bin/bash

read -p "gitlab.cern.ch username: " GITLAB_USERNAME
read -s -p "gitlab.cern.ch Password: " GITLAB_PASSWORD
echo

GITLAB_REPO_URL="https://gitlab.cern.ch/acc-models/acc-models-psb"
FILES=("psb.seq" "psb_aperture.dbx" "scenarios/lhc/1_flat_bottom/psb_fb_lhc.str")
LOCAL_DESTINATION_DIRECTORY="."

for FILE_PATH in "${FILES[@]}"; do
    FILE_URL="$GITLAB_REPO_URL/raw/2021/$FILE_PATH"
    OUTPUT_FILE="$LOCAL_DESTINATION_DIRECTORY/$(basename $FILE_PATH)"
    curl --user "$GITLAB_USERNAME:$GITLAB_PASSWORD" "$FILE_URL" -o "$OUTPUT_FILE"

    echo "File '$FILE_PATH' downloaded to '$OUTPUT_FILE'."
done