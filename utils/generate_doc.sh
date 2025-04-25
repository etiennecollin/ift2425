#!/usr/bin/env bash

domain="ift2425.etiennecollin.com"
docs_path="../docs"
crate_name="utils"

# Check if arg -o is passed
if [ "$1" == "-o" ]; then
    cargo doc --workspace --document-private-items --no-deps --open
else
    cargo doc --workspace --document-private-items --no-deps
fi

# Delete old docs and copy new ones
rm -rf "$docs_path"
mkdir -p "$docs_path"
cp -r ./target/doc/* "$docs_path"

# Remove source code from the docs
# rm -rf "$docs_path/src"

# Add a redirect to the index page
echo "<meta http-equiv=\"refresh\" content=\"0; url=$crate_name/index.html\">" >"$docs_path/index.html"

# Add a CNAME file for GitHub Pages
echo "$domain" >"$docs_path/CNAME"

git add "$docs_path"
git commit -m "Updated docs" "$docs_path"
