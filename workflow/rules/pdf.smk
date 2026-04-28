from pathlib import Path

from legendsimflow import aggregate, patterns
from legendsimflow.metadata import get_tier_settings

_pdf_settings = get_tier_settings(config, "pdf")


rule gen_all_tier_pdf:
    """Aggregate and produce all the `pdf` tier files.

    No wildcards are used.
    """
    input:
        aggregate.gen_list_of_all_tier_pdf_outputs(config),


rule archive_pdfs:
    """Archive all pdf tier files into a single tarball.

    Must be triggered manually with ``snakemake archive_pdfs`` — it is not
    part of the default ``all`` target. Collects all LH5 files produced under
    ``tier/pdf/`` and packs them into ``tarballs/<cycle>-pdfs.tar.xz``,
    preserving the directory tree structure.

    No wildcards are used.
    """
    localrule: True
    input:
        aggregate.gen_list_of_all_tier_pdf_outputs(config),
    output:
        patterns.pdf_tarball_filename(config),
    run:
        from pathlib import Path
        from legendsimflow.archive import create_pdfs_tarball

        create_pdfs_tarball(
            pdf_dir=config.paths.tier.pdf,
            output=Path(output[0]),
            prefix=Path(output[0]).name.removesuffix(".tar.xz"),
        )


rule build_tier_pdf:
    """Produce a `pdf` tier file.

    Reads `cvt` tier data and bins it into histograms (the PDFs) according to
    the PDF configuration file.

    Uses wildcard `simid`.
    """
    message:
        "Producing output file for job pdf.{wildcards.simid}"
    input:
        cvt_file=patterns.output_tier_cvt_filename(config),
    params:
        # surfaced so Snakemake invalidates outputs when the regex map changes
        detector_groups=_pdf_settings.get("detector_groups", None),
    output:
        patterns.output_tier_pdf_filename(config),
    log:
        patterns.log_tier_pdf_filename(config),
    benchmark:
        patterns.benchmark_tier_pdf_filename(config)
    script:
        "../src/legendsimflow/scripts/tier/pdf.py"
