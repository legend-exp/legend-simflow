SHELL := /bin/bash
SOURCEDIR = source
BUILDDIR = build

all: apidoc
	sphinx-build \
      -T \
      -b html \
      -d "$(BUILDDIR)"/doctrees \
      -D language=en \
      -W --keep-going \
      -n \
      "$(SOURCEDIR)" "$(BUILDDIR)"

apidoc: clean-apidoc
	sphinx-apidoc \
      --private \
      --module-first \
      --force \
      --output-dir "$(SOURCEDIR)/api" \
      ../workflow/src \
      ../workflow/src/legendsimflow/scripts
	python tools/generate_rule_docs.py

rulegraph:
	touch ../tests/dummyprod/inputs/simprod/l200cfg01-optmap-dummy.lh5
	python tools/make_rulegraph.py 2>/dev/null > "$(SOURCEDIR)/img/rulegraph.mmd"

clean-apidoc:
	rm -rf "$(SOURCEDIR)/api"

clean: clean-apidoc
	rm -rf "$(BUILDDIR)"
