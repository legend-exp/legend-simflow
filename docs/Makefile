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

clean-apidoc:
	rm -rf "$(SOURCEDIR)/api"

clean: clean-apidoc
	rm -rf "$(BUILDDIR)"
