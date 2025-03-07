# Create mamba environment
mamba create --name cpcr \
    blast biopython gtfparse rich-argparse spyder primer3 primer3-py pysam tqdm \
    -c bioconda -y

# Install cPCR
git clone https://github.com/CherWeiYuan/Panthera.Git
pip install -e .

