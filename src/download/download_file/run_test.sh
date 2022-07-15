

./$meta_functionality_name \
  --input https://raw.githubusercontent.com/scala/scala/2.13.x/NOTICE \
  --output NOTICE

[ ! -f NOTICE ] && echo Output file could not be found && exit 1

if ! grep -q 'Licensed under the Apache License' NOTICE; then
  echo Could not find content
  exit 1
fi

echo Test succeeded!