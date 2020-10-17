Gene search and comparison engine
--- 

This repository contains a simple gene search and comparison engine. The engine parses DNA sequences and generates codon frequency comparison matrix in Phylip format. All six reading frames are being used.

##### Instructions

1. Clone repository and `cd` into it.
2. Create virtual environment `python3 -m venv .env`
3. Log into created virtual environment using `source ./.env/bin/activate`
4. Install dependencies with `pip3 install -r requirements.txt`
5. Upload DNA sequence files into `./dna-data` directory (default format is `fasta`)
6. Change engine config according to your needs in `./config/config.py`
7. Run `python3 ./main.py`.
8. Copy and paste generated Phylip matrix into [Tree Viewer](http://www.trex.uqam.ca/index.php) to inspect clusterization and compare results.
9. Use `deactivate` CLI command to exit virtual environment.
