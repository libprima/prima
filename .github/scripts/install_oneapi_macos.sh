#!/bin/bash

# SPDX-FileCopyrightText: 2020 Intel Corporation
#
# SPDX-License-Identifier: MIT

URL=$1
COMPONENTS=$2

curl --output webimage.dmg --url "$URL" --retry 5 --retry-delay 5
hdiutil attach webimage.dmg
if [ -z "$COMPONENTS" ]; then
  sudo /Volumes/"$(basename "$URL" .dmg)"/bootstrapper.app/Contents/MacOS/bootstrapper -s --action install --eula=accept --continue-with-optional-error=yes --log-dir=.
  installer_exit_code=$?
else
  sudo /Volumes/"$(basename "$URL" .dmg)"/bootstrapper.app/Contents/MacOS/bootstrapper -s --action install --components="$COMPONENTS" --eula=accept --log-dir=.
  installer_exit_code=$?
fi
hdiutil detach /Volumes/"$(basename "$URL" .dmg)" -quiet
exit $installer_exit_code
