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


def test_despine_numeric_axis_ends_on_ticks():
    from _style import despine
    plt = _fig()
    fig, ax = plt.subplots()
    ax.plot([0.3, 9.2], [1.1, 7.7])
    despine(ax)
    for (lo, hi), ticks, (dmin, dmax) in [
        (ax.get_xlim(), ax.get_xticks(), (0.3, 9.2)),
        (ax.get_ylim(), ax.get_yticks(), (1.1, 7.7)),
    ]:
        assert (lo, hi) == (ticks[0], ticks[-1])
        assert lo <= dmin and hi >= dmax
    assert tuple(ax.spines["bottom"].get_bounds()) == ax.get_xlim()
    assert tuple(ax.spines["left"].get_bounds()) == ax.get_ylim()
    plt.close(fig)


def test_despine_date_axis():
    import datetime as dt

    from _style import despine
    plt = _fig()
    import matplotlib.dates as mdates
    fig, ax = plt.subplots()
    days = [dt.datetime(2026, 8, 27) + dt.timedelta(days=i) for i in range(30)]
    ax.plot(days, range(30))
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO))
    despine(ax)
    lo, hi = ax.get_xlim()
    assert lo <= mdates.date2num(days[0]) and hi >= mdates.date2num(days[-1])
    assert tuple(ax.spines["bottom"].get_bounds()) == (lo, hi)
    assert all(mdates.num2date(t).weekday() == 0 for t in ax.get_xticks())
    plt.close(fig)


def test_despine_unclips_fitted_axes():
    from _style import despine
    plt = _fig()
    fig, ax = plt.subplots()
    (line,) = ax.plot([0, 1, 2, 3], [0, 1, 4, 9], marker="o")
    despine(ax)
    assert line.get_clip_on() is False
    plt.close(fig)


def test_despine_keeps_clipping_with_explicit_limits():
    from _style import despine
    plt = _fig()
    fig, ax = plt.subplots()
    (line,) = ax.plot([0, 1, 2, 3], [0, 1, 4, 9])
    ax.set_xlim(1, 2)
    despine(ax)
    assert line.get_clip_on() is True
    plt.close(fig)


def test_despine_categorical_y_drops_left_spine():
    from _style import despine
    plt = _fig()
    fig, ax = plt.subplots()
    ax.barh(["a", "b", "c"], [3.2, 1.0, 7.5])
    despine(ax, categorical_y=True)
    assert not ax.spines["left"].get_visible()
    assert all(t.tick1line.get_markersize() == 0 for t in ax.yaxis.get_major_ticks())
    assert ax.get_xlim()[1] >= 7.5
    assert tuple(ax.spines["bottom"].get_bounds()) == ax.get_xlim()
    plt.close(fig)
