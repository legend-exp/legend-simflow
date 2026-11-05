from pathlib import Path

from snakemake.logging import logger

from legendsimflow import aggregate, patterns, utils
from legendsimflow.metadata import get_tier_settings

_pdf_settings = get_tier_settings(config, "pdf")

if _pdf_settings.build_l200data_pdfs and config.paths.get("l200data") is None:
    logger.warning(
        "build_l200data_pdfs is set but paths.l200data is not: "
        "the LEGEND-200 data histograms are not built"
    )


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
        detector_groups=_pdf_settings.detector_groups,
    output:
        patterns.output_tier_pdf_filename(config),
    log:
        patterns.log_tier_pdf_filename(config),
    benchmark:
        patterns.benchmark_tier_pdf_filename(config)
    script:
        "../src/legendsimflow/scripts/tier/pdf.py"


rule build_pdf_from_l200data:
    """Histogram the LEGEND-200 data of a run in the format of the `pdf` tier.

    Reads the `evt` (or `pet`) tier files of the run from `paths.l200data`,
    selects the events as `build_tier_pdf` does (see
    {mod}`legendsimflow.selections`) and writes the histograms of each
    category and detector group in the same layout, with the livetime of the
    run in place of the number of simulated events. The A/E of the detectors
    that the simulation does not model is invalid, as in the simulation.
    `merge_l200data_pdfs` sums the runs. Only with `paths.l200data` and the
    `build_l200data_pdfs` setting of the `pdf` tier.

    Uses wildcard `runid`.
    """
    message:
        "Histogramming the data of run {wildcards.runid} in the pdf tier format"
    input:
        evt_files=lambda wildcards: utils.lookup_evt_files(
            config.paths.l200data,
            wildcards.runid,
            utils.get_evt_tier_name(config.paths.l200data),
        ),
        # params.modeled_detectors reads the modelable-HPGe checkpoint
        _hpge_cache=smk_wait_for_hpge_cache,
    params:
        modeled_detectors=lambda wildcards: sorted(
            smk_load_hpge_cache()[wildcards.runid]
        ),
        detector_groups=_pdf_settings.detector_groups,
    output:
        patterns.output_l200data_pdf_filename(config),
    log:
        patterns.log_l200data_pdf_filename(config),
    script:
        "../src/legendsimflow/scripts/build_pdf_from_l200data.py"


rule merge_l200data_pdfs:
    """Sum the data histograms of all runs of the Simflow into one `pdf` tier file.

    Adds the histograms of `build_pdf_from_l200data` over the runs, and their
    livetimes, with {func}`legendsimflow.pdf.merge_pdf_files`. The output
    sits next to the simulated PDFs, which cover the same runs.

    No wildcards are used.
    """
    message:
        "Merging the data histograms of all runs"
    localrule: True
    input:
        aggregate.gen_list_of_all_l200data_pdf_outputs(config),
    output:
        patterns.output_merged_l200data_pdf_filename(config),
    run:
        from legendsimflow.pdf import merge_pdf_files

        merge_pdf_files(input, output[0])
