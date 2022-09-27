#!/bin/bash

# location of run script to be created
TARGET_FILE="run_pipeline.sh"

echo -e "\n~~~ yEvo Pipeline v1.0 - Generate Run Script ~~~\n"
echo -e "\nCreating run script...\n"
sed 's:{{BASE_DIR}}:'`pwd`':' scripts/templates/_run_pipeline.sh > $TARGET_FILE

echo -e "\nCreated $TARGET_FILE\n"
echo -e "\nChanging file permissions to make run script executable...\n"
chmod +x $TARGET_FILE

# success
echo -e "\nDONE!\n"
