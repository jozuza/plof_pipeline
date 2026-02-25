

# ### CLOUD
# !/bin/bash

############################################################################
# Optional: Install Google Cloud SDK (includes gsutil)
# This section allows the script to copy files from Google Cloud Storage
# directly into the container when an input VCF is missing.
############################################################################
## new Oct 23
# Update package lists and install basic utilities required for setup.
# - curl: used to fetch files from the internet (e.g., the Google signing key)
# - lsb-release: provides Linux distribution information (needed by some installers)
# - gnupg: allows managing encryption keys and verifying package authenticity
apt-get update && apt-get install -y curl lsb-release gnupg

# Add the official Google Cloud SDK repository to apt sources.
# This repository provides the 'google-cloud-sdk' package, which includes 'gsutil'.
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

# Download and add Google's public signing key to verify the authenticity
# of packages downloaded from the Cloud SDK repository.
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -

# Update package lists again (to include the new Google repo),
# and install the full Google Cloud SDK, which contains 'gsutil'.
# 'gsutil' is used to copy files from Google Cloud Storage (gs://... paths).
apt-get update && apt-get install -y google-cloud-sdk
##
