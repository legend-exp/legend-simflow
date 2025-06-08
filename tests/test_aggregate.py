from legendsimflow import aggregate as agg


def test_gen_lists(simflow_config):
    assert agg.gen_list_of_plots_outputs(simflow_config, "stp", "test-simid") == []
