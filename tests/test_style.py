"""save_figure must emit both formats and preserve dotted output stems."""


def _fig():
    from _style import apply_style
    apply_style()
    import matplotlib.pyplot as plt
    return plt


def test_dotted_stem_keeps_all_dots(tmp_path):
    from _style import save_figure
    plt = _fig()
    fig = plt.figure()
    written = save_figure(fig, str(tmp_path / "K562_DNASE.epochs.svg"))
    plt.close(fig)
    assert sorted(p.name for p in written) == [
        "K562_DNASE.epochs.png", "K562_DNASE.epochs.svg"]
    assert (tmp_path / "K562_DNASE.epochs.svg").exists()
    assert (tmp_path / "K562_DNASE.epochs.png").exists()


def test_suffixless_path_gets_both_formats(tmp_path):
    from _style import save_figure
    plt = _fig()
    fig = plt.figure()
    save_figure(fig, str(tmp_path / "peaks_overview"))
    plt.close(fig)
    assert (tmp_path / "peaks_overview.svg").exists()
    assert (tmp_path / "peaks_overview.png").exists()
