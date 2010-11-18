root_dir=`pwd`
for d in '.' 'contrib/aracne' 'contrib/tinyxml' 'conan' 'test' 'python'; do
  cd ${root_dir}/${d}
  aclocal && autoconf && automake --add-missing --copy && autoconf
  rm -rf autom4te.cache
done
