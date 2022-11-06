# Need update `exports.h` by running `bash update.sh` manually.

c = run_command('bash', '-c', 'cd "$MESON_SUBDIR" ; ls *.cpp', check: true)
cqlat_utils_sources = files(c.stdout().strip().split('\n'))
