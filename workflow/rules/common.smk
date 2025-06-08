def gen_target_all():
    if config.get("simlist", "*") in ("all", "*"):
        if "pdf" in make_tiers:
            return rules.gen_pdf_release.output
        elif "evt" in make_tiers:
            return rules.gen_all_tier_evt.output
        elif "hit" in make_tiers:
            return rules.gen_all_tier_hit.output
        elif "stp" in make_tiers:
            return (
                rules.gen_all_tier_stp.output,
                aggregate.gen_list_of_all_plots_outputs(config, tier="stp"),
                aggregate.gen_list_of_all_plots_outputs(config, tier="ver"),
            )
    else:
        return aggregate.process_simlist(config)
