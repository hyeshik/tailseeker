#!/bin/sh
TAILSEEKER_VERSION=$(env PYTHONPATH=.. python -c 'from tailseeker import __version__; print(__version__)')
PACKAGE_NAME=tailseeker-${TAILSEEKER_VERSION}

if [ -d "$PACKAGE_NAME" ]; then
  rm -rf "$PACKAGE_NAME"
fi

umask 022

git clone https://github.com/hyeshik/tailseeker.git "${PACKAGE_NAME}"
(cd "${PACKAGE_NAME}" && git checkout -q "tags/v${TAILSEEKER_VERSION}")

# Remove unnessary files from the tree
(cd "$PACKAGE_NAME" && rm -rf .git .gitignore)

# Make a tarball
tar -czf "${PACKAGE_NAME}.tar.gz" "${PACKAGE_NAME}"
rm -rf "${PACKAGE_NAME}"

# Sign the package
gpg -u EBB09D8E --output ${PACKAGE_NAME}.tar.gz.sig --detach-sig \
  ${PACKAGE_NAME}.tar.gz
