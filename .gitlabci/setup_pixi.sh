#!/bin/bash -e

echo "+ installing pixi..."
curl -fsSL https://minio.retd.edf.fr/codeaster/tools/pixi-install.sh | bash


echo "+ setting config.toml..."
mkdir -p ${PIXI_HOME}
cat << EOF > ${PIXI_HOME}/config.toml
[mirrors]
"https://conda.anaconda.org/conda-forge" = ["https://nexus.retd.edf.fr/repository/conda-forge/conda-forge/"]

[pypi-config]
index-url = "https://nexus.retd.edf.fr/repository/pypi-all/simple"
allow-insecure-host = ["nexus.retd.edf.fr"]
EOF

pixi --version
pixi config set --global run-post-link-scripts insecure
pixi global install git minio-client
