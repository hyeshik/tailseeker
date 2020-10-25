export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig

make -C $SRC_DIR/src

# Copy tailseeker dir
INSTALL_SUBDIRS="bin conf docs refdb scripts tailseeker templates"
TARGET_DIR="$PREFIX/share/tailseeker"
mkdir -p m 755 $TARGET_DIR

for subdir in $INSTALL_SUBDIRS; do
  cp -Rp $SRC_DIR/$subdir $TARGET_DIR/$subdir
done
rm -f $TARGET_DIR/bin/tailseq-docker-wrap

# Install entrypoint script.
sed -e "s|%%PREFIX%%|$PREFIX|g" \
  -e "s|%%VERSION%%|$PKG_VERSION conda:$PKG_VERSION-$PKG_BUILDNUM|g" \
  $SRC_DIR/conda/tseek.in > \
  $PREFIX/bin/tseek
chmod 755 $PREFIX/bin/tseek

# Set paths
sed -e "s|%%PREFIX%%|$PREFIX|g" $SRC_DIR/conda/paths.conf.in > \
  $TARGET_DIR/conf/paths.conf
chmod 644 $TARGET_DIR/conf/paths.conf

