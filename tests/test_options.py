import sys

import pytest

from dsr_shelx.options import OptionsParser


def parse(monkeypatch, *args):
    monkeypatch.setattr(sys, 'argv', ['dsr', *args])
    return OptionsParser('')


@pytest.mark.parametrize('args, option, value', [
    (['-r', 'my', 'file.res'], 'res_file', 'my file.res'),
    (['-re', 'p21c.res'], 'external_restr', 'p21c.res'),
    (['-target', '1', '2', '3'], 'target_coords', [1.0, 2.0, 3.0]),
    (['-s', 'benz'], 'search_string', 'benz'),
    (['-x', 'tol', 'uene'], 'search_extern', 'toluene'),
    (['-ah', ' Toluene '], 'head_for_gui', 'toluene'),
    (['-c', ' benzene '], 'export_clip', 'benzene'),
    (['-e', 'benzene'], 'export_fragment', 'benzene'),
    (['-i', 'ALA.tgz'], 'import_grade', 'ALA.tgz'),
    (['-shx', '/does/not/exist/shelxl'], 'shelxl_ex', ''),
    (['-t'], 'invert', True),
    (['-g'], 'rigid_group', True),
    (['-n'], 'no_refine', True),
    (['-l'], 'list_db', True),
    (['-lc'], 'list_db_csv', True),
    (['-ea'], 'export_all', True),
    (['-u'], 'selfupdate', True),
    ([], 'res_file', False),
])
def test_options(monkeypatch, args, option, value):
    assert getattr(parse(monkeypatch, *args), option) == value


@pytest.mark.parametrize('args, option', [
    (['-target', '1', '2'], 'target_coords'),
    (['-s', 'benz;rm'], 'search_string'),
    (['-x', 'a/b'], 'search_extern'),
])
def test_invalid_options_exit(monkeypatch, args, option):
    options = parse(monkeypatch, *args)
    with pytest.raises(SystemExit):
        getattr(options, option)


def test_error_prints_help(monkeypatch, capsys):
    with pytest.raises(SystemExit):
        parse(monkeypatch).error()
    assert 'Disordered Structure Refinement' in capsys.readouterr().out
