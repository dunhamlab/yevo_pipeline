#!/bin/bash

# test data URL
DATA_URL="https://shared-misc.s3.us-east-2.amazonaws.com/yevo_pipeline/v1.0/test_data.zip"

echo -e "\n~~~ yEvo Pipeline v1.0 - Test Data Download Script ~~~\n"
echo -e "\nDownloading compressed test data...\n"
curl $DATA_URL -o test_data.zip

echo -e "\nUnzipping downloaded test data...\n"
unzip test_data.zip

echo -e "\nCopying downloaded files to destination...\n"
cp -r test_data/* data

echo -e "\nCleaning up...\n"
rm test_data.zip
rm -rf test_data/

# success
echo -e "\nDONE!\n"
