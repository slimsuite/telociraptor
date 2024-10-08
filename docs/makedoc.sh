#export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

python3 ../code/telociraptor.py dochtml newlog

python3 ../code/telociraptor.py --description | sed 's/^/# /' > ../Telociraptor.md
echo >> ../Telociraptor.md

echo '```' >> ../Telociraptor.md
python3 ../code/telociraptor.py --details >> ../Telociraptor.md
echo '```' >> ../Telociraptor.md
echo >> ../Telociraptor.md
echo 'For a better rendering and navigation of this document, please download and open [`./docs/telociraptor.docs.html`](./docs/telociraptor.docs.html), or visit <https://slimsuite.github.io/telociraptor/>.' >> ../Telociraptor.md
echo 'Documentation can also be generated by running Telociraptor with the `dochtml=T` option. (R and pandoc must be installed - see below.)' >> ../Telociraptor.md
echo >> ../Telociraptor.md
echo '## Introduction' >> ../Telociraptor.md
echo >> ../Telociraptor.md
grep -A 10000 "Telociraptor is a genome" telociraptor.docs.Rmd >> ../Telociraptor.md

cp telociraptor.docs.html ../index.html
cp ../Telociraptor.md ../README.md
