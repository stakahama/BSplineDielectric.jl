#/bin/sh

[ -f ./README.md ] && rm -f ../README.md && mv ./README.md ..; \
    [ -d ./README_files ] && rm -rf ../README_files && mv ./README_files ..
