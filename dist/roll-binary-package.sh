#!/bin/sh
TAILSEEKER_VERSION=$(env PYTHONPATH=.. python -c 'from tailseeker import __version__; print(__version__)')
PLATFORM=ubuntu_xenial
PACKAGE_NAME=tailseeker-${TAILSEEKER_VERSION}-bundle-${PLATFORM}

if [ -d "$PACKAGE_NAME" ]; then
  rm -rf "$PACKAGE_NAME"
fi

# Extract files from the latest docker image
CONTAINER_ID=$(docker create hyeshik/tailseeker:latest)
docker cp $CONTAINER_ID:/opt/tailseeker $PACKAGE_NAME
docker rm -v $CONTAINER_ID

# Remove unnessary files from the tree
(cd "$PACKAGE_NAME" && \
  rm -rf .git .gitignore include share refdb/level2/.snakemake \
    refdb/level3/.snakemake lib/pkgconfig lib/libhts.a \
    tailseeker/__pycache__ tailseeker/.snakemake \
)
find "$PACKAGE_NAME" -name '*.o' -delete

# Convert absolute paths to tailseeker binaries to relative
perl -pi -e 's,/opt/tailseeker/,,' "$PACKAGE_NAME/conf/paths.conf"

# Replace the entry script to use other binaries from the relative
# paths
cp -pf "$PACKAGE_NAME/install/tailseeker.in" \
       "$PACKAGE_NAME/bin/tseek"

# Overlay additional documents for binary releases
tar -C files -cf - . | tar -C "$PACKAGE_NAME" -xf -

# Fix permissions
chmod go+r "$PACKAGE_NAME"/LICENSE*
chmod go+rX "$PACKAGE_NAME"

# Make a tarball
tar -czf "${PACKAGE_NAME}.tar.gz" "${PACKAGE_NAME}"
rm -rf "${PACKAGE_NAME}"

# Sign the package
gpg -u EBB09D8E --output ${PACKAGE_NAME}.tar.gz.sig --detach-sig \
  ${PACKAGE_NAME}.tar.gz
