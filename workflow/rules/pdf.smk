

rule gen_all_tier_pdf:
    """Aggregate and produce all the pdf tier files."""
    input:
        aggregate.gen_list_of_all_tier_pdf_outputs(config),


rule gen_pdf_release:
    """Generate a compressed archive with all the PDF files.

    Pack all PDF files into a `.tar.xz` archive using `tar --create --xz`,
    renaming them to a flat `{experiment}-pdfs/` directory structure.

    No wildcards are used.
    """
    message:
        "Generating pdf release"
    input:
        aggregate.gen_list_of_all_tier_pdf_outputs(config),
    params:
        exp=config["experiment"],
        ro_input=lambda wildcards, input: utils.as_ro(config, input),
    output:
        Path(config["paths"]["pdf_releases"]) / (config["experiment"] + "-pdfs.tar.xz"),
    shell:
        r"""
        tar --create --xz \
            --file {output} \
            --transform 's|.*/\({params.exp}-.*-tier_pdf\..*\)|{params.exp}-pdfs/\1|g' \
            {params.ro_input}
        """


rule build_tier_pdf:
    """Produce a `pdf` tier file.

    Run the `build-pdf` command, which reads `evt` tier data and bins it into
    histograms (the PDFs) according to the PDF configuration file.

    Uses wildcard `simid`.
    """
    message:
        "Producing output file for job pdf.{wildcards.simid}"
    input:
        cvt_file=patterns.output_tier_cvt_filename(config),
    output:
        patterns.output_tier_pdf_filename(config),
    log:
        patterns.log_tier_pdf_filename(config),
    benchmark:
        patterns.benchmark_tier_pdf_filename(config)
    script:
        "../src/legendsimflow/scripts/tier/pdf.py"
