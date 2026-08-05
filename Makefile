# MS2Bac Make pipeline
# Usage examples:
#   make help
#   make setup
#   make database ANI_FILE=00_metafiles/ANI_report_saccharomycotina.txt
#   make mapping
#   make prepare-search
#   make search RAW=/absolute/path/sample.raw
#   make report RAW=/absolute/path/sample.raw EXPERIMENT=sample

SHELL := /bin/bash
.DELETE_ON_ERROR:

ROOT ?= $(CURDIR)
PYTHON ?= python
CONDA_ENV ?= MS2Bac
CONDA_RUN ?= conda run -n $(CONDA_ENV)
DOCKER ?= docker
DOCKER_IMAGE ?= nava20ir/fragpipebbm:latest
N_CHUNKS ?= 10
JAVA_MEM ?= 32G

META_DIR := $(ROOT)/00_metafiles
FASTA_DIR := $(ROOT)/04_fasta
SEARCH_DIR := $(FASTA_DIR)/search_space
RESULTS_DIR ?= $(ROOT)/results
IDB ?= $(ROOT)/idb.pkl.gz
BACCI ?= $(ROOT)/bacci_finder.pkl.gz
ANI_FILE ?= $(META_DIR)/ANI_report_prokaryotes.txt
ASSEMBLY_SUMMARY ?= $(META_DIR)/assembly_summary.txt

SCRIPTS := download_database_MS2Bac.py mapping_file_generator.py \
           make_fragger_parameter_file.py chunk_fasta_file.py \
           combine_fragger_results.py identification.py visualize_output.py

RAW ?=
EXPERIMENT ?= $(if $(RAW),$(basename $(notdir $(RAW))),experiment)
RAW_BASE = $(basename $(notdir $(RAW)))
RUN_DIR = $(RESULTS_DIR)/$(EXPERIMENT)
PSM = $(RUN_DIR)/psm.csv
OUTPUT_TABLE = $(RUN_DIR)/output_table.csv
SUMMARY = $(RUN_DIR)/summary_output.csv
PDF = $(RUN_DIR)/$(EXPERIMENT)_identification.pdf

.PHONY: help setup pull-image metadata database mapping prepare-search \
        check-raw search combine identify visualise report clean-search clean-results

help:
	@echo "MS2Bac targets:"
	@echo "  make setup"
	@echo "  make database ANI_FILE=00_metafiles/ANI_report_saccharomycotina.txt"
	@echo "  make mapping"
	@echo "  make prepare-search N_CHUNKS=10"
	@echo "  make search RAW=/absolute/path/sample.raw EXPERIMENT=sample"
	@echo "  make report RAW=/absolute/path/sample.raw EXPERIMENT=sample"
	@echo "  make visualise EXPERIMENT=sample"

setup: linux_conda.yml
	conda env create -f $<

pull-image:
	$(DOCKER) pull $(DOCKER_IMAGE)

$(META_DIR):
	mkdir -p $@

metadata: $(ASSEMBLY_SUMMARY)

$(ASSEMBLY_SUMMARY): | $(META_DIR)
	curl -L --fail --retry 3 -o $@ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# ANI_FILE is intentionally not downloaded automatically because fungal and
# prokaryotic workflows use different ANI reports. Put the selected ANI report
# in 00_metafiles and pass ANI_FILE=... on the command line.
database: $(ASSEMBLY_SUMMARY)
	@test -f "$(ANI_FILE)" || { echo "Missing ANI file: $(ANI_FILE)"; exit 2; }
	$(CONDA_RUN) $(PYTHON) download_database_MS2Bac.py "$(ROOT)" "$(ANI_FILE)"

mapping:
	$(CONDA_RUN) $(PYTHON) mapping_file_generator.py "$(ROOT)"
	@test -f "$(IDB)" -o -f "$(ROOT)/04_fasta/idb.pkl.gz" -o -f "$(ROOT)/03_digest/idb.pkl.gz" || \
	  echo "NOTE: locate idb.pkl.gz and set IDB=/absolute/path/idb.pkl.gz"
	@test -f "$(BACCI)" -o -f "$(ROOT)/04_fasta/bacci_finder.pkl.gz" -o -f "$(ROOT)/03_digest/bacci_finder.pkl.gz" || \
	  echo "NOTE: locate bacci_finder.pkl.gz and set BACCI=/absolute/path/bacci_finder.pkl.gz"

prepare-search:
	mkdir -p "$(SEARCH_DIR)"
	@test -n "$$(find "$(FASTA_DIR)" -maxdepth 1 -name '*.fasta' -print -quit)" || { echo "No FASTA files in $(FASTA_DIR)"; exit 2; }
	cat "$(FASTA_DIR)"/*.fasta > "$(SEARCH_DIR)/all_fasta.fasta"
	cp make_fragger_parameter_file.py chunk_fasta_file.py "$(SEARCH_DIR)/"
	cd "$(SEARCH_DIR)" && $(PYTHON) chunk_fasta_file.py all_fasta.fasta $(N_CHUNKS)
	rm -f "$(SEARCH_DIR)/all_fasta.fasta"

check-raw:
	@test -n "$(RAW)" || { echo "RAW is required, e.g. make report RAW=/path/sample.raw EXPERIMENT=sample"; exit 2; }
	@test -f "$(RAW)" || { echo "RAW file not found: $(RAW)"; exit 2; }
	@test -n "$$(find "$(SEARCH_DIR)" -maxdepth 1 -name '*.fasta' -print -quit)" || { echo "No search FASTAs. Run make prepare-search"; exit 2; }

# Run one independent MSFragger search for every FASTA chunk.
# Results go into results/<experiment>, so *.tsv from another experiment are
# not accidentally combined.
search: check-raw pull-image
	mkdir -p "$(RUN_DIR)"
	cp "$(RAW)" "$(RUN_DIR)/$(notdir $(RAW))"
	@set -euo pipefail; \
	for fasta in "$(SEARCH_DIR)"/*.fasta; do \
	  name=$$(basename "$$fasta"); \
	  cp "$$fasta" "$(RUN_DIR)/$$name"; \
	  $(PYTHON) make_fragger_parameter_file.py "$(RUN_DIR)/$$name" "$(RUN_DIR)/$$name.params"; \
	  echo "Running MSFragger with $$name"; \
	  $(DOCKER) run --rm \
	    -v "$(RUN_DIR):/fragger_projects" \
	    $(DOCKER_IMAGE) \
	    java -Xmx$(JAVA_MEM) -jar /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/MSFragger-4.4.1/MSFragger-4.4.1/MSFragger-4.4.1.jar \
	    "/fragger_projects/$$name.params" "/fragger_projects/$(notdir $(RAW))"; \
	  test -f "$(RUN_DIR)/$(RAW_BASE).tsv"; \
	  mv "$(RUN_DIR)/$(RAW_BASE).tsv" "$(RUN_DIR)/$(RAW_BASE).$$name.tsv"; \
	done

combine: search
	$(PYTHON) combine_fragger_results.py "$(RUN_DIR)" "$(EXPERIMENT)"
	@test -s "$(PSM)"

# identification.py writes fixed filenames in its current working directory,
# therefore run it from RUN_DIR using absolute paths.
identify: combine
	@test -f "$(IDB)" || { echo "Missing IDB: $(IDB)"; exit 2; }
	@test -f "$(BACCI)" || { echo "Missing bacci finder: $(BACCI)"; exit 2; }
	cd "$(RUN_DIR)" && $(PYTHON) "$(ROOT)/identification.py" "$(PSM)" "$(IDB)" "$(BACCI)"
	@test -s "$(OUTPUT_TABLE)"

visualise:
	@test -f "$(OUTPUT_TABLE)" || { echo "Missing $(OUTPUT_TABLE); run make identify ... first"; exit 2; }
	$(PYTHON) visualize_output.py -i "$(OUTPUT_TABLE)" -o "$(PDF)"
	@test -s "$(PDF)"

report: identify visualise
	@echo "Finished: $(PDF)"

clean-search:
	rm -rf "$(SEARCH_DIR)"

clean-results:
	@test -n "$(EXPERIMENT)" && rm -rf "$(RUN_DIR)"