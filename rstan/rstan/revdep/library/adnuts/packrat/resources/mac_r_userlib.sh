#!/bin/bash

set -e

# R system library to user library migration script for Mac OS X
#
# Date: January 14, 2014
# Author: Joe Cheng <joe@rstudio.com>
#
# From https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html:
#   The official CRAN binaries come pre-packaged in such a way that
#   administrator have sufficient privileges to update R and install
#   packages system-wide.
#
# This means that any install.packages() call, or using Install Package
# from RStudio, causes packages to be installed in the system library
# (e.g. /Library/Frameworks/R.framework/Versions/3.0/Resources/library).
# The system library contains base and recommended packages as well.
#
# We believe it's more hygienic to keep base/recommended packages
# separate from user-installed packages, and this separation is
# necessary for the Packrat[0] dependency management system to provide
# isolation benefits.
#
# This script creates a personal library directory, and migrates any
# non-base, non-recommended packages from the system library into
# it. It then sets the permissions on the system library to only be
# writable by root. This will ensure that future install.packages calls
# will not add more packages to the system library.
#
# [0] https://rstudio.github.io/packrat/


# The system-wide library
RLIBDIR=`R --vanilla --slave -e "cat(tail(.libPaths(), 1))"`

# The user library (which might not exist yet)
RLIBSUSER=`R --vanilla --slave -e "cat(path.expand(head(Sys.getenv('R_LIBS_USER'), 1)))"`

# The list of non-base, non-recommended packages in the system-wide library
PKGS=`R --vanilla --slave -e "cat(with(as.data.frame(installed.packages(tail(.libPaths(), 1))), paste(Package[is.na(Priority)])))"`

if [ "$RLIBDIR" == "" ]; then
  echo "ERROR: Couldn't detect system library directory, aborting" >&2
  exit 1
fi

if [ "$RLIBSUSER" == "" ]; then
  echo "ERROR: Couldn't detect R_LIBS_USER directory, aborting" >&2
  exit 1
fi

echo "Saving backup of $RLIBDIR to ./SystemLibBackup.tar.gz"
if [ -f ./SystemLibBackup.tar.gz ]; then
  echo "SystemLibBackup.tar.gz exists. Press Enter to overwrite, or Ctrl-C to abort:" >&2
  read -s < /dev/tty
  echo "Backing up..."
fi
tar -czPf SystemLibBackup.tar.gz "$RLIBDIR"
#tar -czf SystemLibBackup.tar.gz -C "$RLIBDIR" $(ls $RLIBDIR)
echo "Backup successful."

echo "Migrating user-installed packages to $RLIBSUSER"
echo "Press Enter to continue, or Ctrl-C to abort"
read -s < /dev/tty

mkdir -p "$RLIBSUSER"

for pkg in $PKGS; do
  echo "Moving $pkg"
  if [ -d "$RLIBSUSER/$pkg" ]; then
    echo "ERROR: The directory $RLIBSUSER/$pkg already exists, aborting" >&2
    echo "Please delete the package $pkg from either $RLIBDIR or $RLIBSUSER."
    exit 3
  fi
  # Do a copy to get default permissions
  cp -R "$RLIBDIR/$pkg" "$RLIBSUSER"
  sudo rm -rf "$RLIBDIR/$pkg"
done

echo
echo "Making $RLIBDIR writable only by root (chmod 755)"
sudo chmod -R 755 "$RLIBDIR"

echo
echo Success!
