import re
import shutil
import unittest
from pathlib import Path
from types import SimpleNamespace

import pytest

from dsr_shelx.atomhandling import get_atoms, get_atomtypes, Elem_2_Sfac, FindAtoms, remove_emptied_afix_blocks, \
    replace_after_fit, NumberScheme, SfacTable
from dsr_shelx.resfile import ResList

TESTS_DIR = Path(__file__).parent


def test_get_atoms():
    result = get_atoms(['F8    4    0.349210   0.073474   0.519443 -21.00000   0.03106', 'foo'])
    assert result == [['F8', '4', '0.349210', '0.073474', '0.519443']]


def test_get_atomtypes():
    atoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'],
             ['C1', 1, '-0.00146', '0.26814', '0.06351'],
             ['Au15', 8, '0.00146', '0.26814', '0.06351']]
    assert get_atomtypes(atoms) == ['O', 'C', 'AU']


def test_get_atomtypes_invalid_element():
    atoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'],
             ['C1', 1, '-0.00146', '0.26814', '0.06351'],
             ['Ax15', 8, '0.00146', '0.26814', '0.06351']]
    with pytest.raises(KeyError):
        get_atomtypes(atoms)


def test_is_atom():
    result = FindAtoms.is_atom(atomline='O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
    assert result == ['O1', '3', '0.120080', '0.336659', '0.494426']


def test_is_atom_without_sfac_number():
    assert FindAtoms.is_atom(atomline='O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...') == []


def test_get_partnumber():
    assert FindAtoms.get_partnumber(partstring='PART 2 -21') == 2


@pytest.mark.parametrize('element, sfac', [('Br', 4), ('c', 3), ('Ar', None)])
def test_elem_2_sfac(element, sfac):
    assert Elem_2_Sfac(['Ca', 'H', 'C', 'BR']).elem_2_sfac(element) == sfac


@pytest.mark.parametrize('sfac, element', [(1, 'CA'), (2, 'H'), (3, 'C'), (4, 'BR'), (5, None), (0, None)])
def test_sfac_2_elem(sfac, element):
    assert Elem_2_Sfac(['Ca', 'H', 'C', 'BR']).sfac_2_elem(sfac) == element


def test_replace_after_fit_removes_hydrogens_and_their_afix(tmp_path, capsys):
    res_file = tmp_path / 'p21c.res'
    shutil.copy(TESTS_DIR / 'p21c.res', res_file)
    lines = res_file.read_text(encoding='latin1').splitlines(keepends=True)
    coords = {}
    for line in lines:
        match = re.match(r'^(C3[0-5])\s+1\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if match:
            coords[match.group(1)] = match.group(2, 3, 4)
    # Put a fitted "benzene" exactly onto C30 to C35:
    hklf = [num for num, line in enumerate(lines) if line.startswith('HKLF')][0]
    lines[hklf:hklf] = ['C{} 1 {} {} {} 11.0 0.03\n'.format(n + 1, *coords['C3{}'.format(n)]) for n in range(6)]
    res_file.write_text(''.join(lines), encoding='latin1')
    rl = ResList(str(res_file))
    reslist = rl.get_res_list()
    cell = [float(x) for x in reslist[1].split()[2:8]]
    resi = SimpleNamespace(dsrp=SimpleNamespace(resiflag=False))
    reslist, _ = replace_after_fit(rl, reslist, resi, ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'], cell)
    assert 'H30 H32 H34' in capsys.readouterr().out
    atom_names = [line.split()[0] for line in reslist if line.strip()]
    for hydrogen in ('H30', 'H32', 'H34'):
        assert hydrogen not in atom_names
    assert len([line for line in reslist if re.match(r'AFIX\s+43', line)]) == 3
    # No empty AFIX mn / AFIX 0 pairs are left behind:
    for num, name in enumerate(atom_names[:-1]):
        if name == 'AFIX' and atom_names[num + 1] == 'AFIX':
            pytest.fail('Empty AFIX block found.')


def test_remove_emptied_afix_blocks_removes_empty_afix_and_afix0():
    res = ['C34 1 0.1 0.2 0.3 11 0.02\n', 'AFIX  43\n', '', 'AFIX   0\n', 'C33 1 0.1 0.2 0.3 11 0.02\n']
    assert remove_emptied_afix_blocks(res) == ['C34 1 0.1 0.2 0.3 11 0.02\n', '', '', '',
                                               'C33 1 0.1 0.2 0.3 11 0.02\n']


def test_remove_emptied_afix_blocks_keeps_afix_with_atoms():
    res = ['AFIX  43\n', 'H34 2 0.1 0.2 0.3 11 -1.2\n', 'AFIX   0\n']
    assert remove_emptied_afix_blocks(res) == ['AFIX  43\n', 'H34 2 0.1 0.2 0.3 11 -1.2\n', 'AFIX   0\n']


def test_remove_emptied_afix_blocks_followed_by_other_afix():
    res = ['AFIX 137\n', '', '', 'AFIX  43\n', 'H1 2 0.1 0.2 0.3 11 -1.2\n', 'AFIX 0\n']
    assert remove_emptied_afix_blocks(res) == ['', '', '', 'AFIX  43\n', 'H1 2 0.1 0.2 0.3 11 -1.2\n',
                                               'AFIX 0\n']


def test_find_atoms_to_replace_keeps_fragment_atoms_in_residue_0():
    # A fitted acetonitrile in PART 0 without residue. N1A-C1A is shorter than the replace distance.
    reslist = ['SFAC C N\n', 'UNIT 1 1\n',
               'C5    1    0.300000    0.410000    0.300000    11.00000    0.04\n',
               'C6    1    0.800000    0.800000    0.800000    11.00000    0.04\n',
               'N1A   2    0.300000    0.300000    0.300000    11.00000    0.04\n',
               'C1A   1    0.415000    0.300000    0.300000    11.00000    0.04\n',
               'C2A   1    0.561000    0.300000    0.300000    11.00000    0.04\n',
               'HKLF 4\n']
    fa = FindAtoms(reslist)
    to_delete = fa.find_atoms_to_replace(['N1A', 'C1A', 'C2A'], [10, 10, 10, 90, 90, 90], remdist=1.3)
    assert to_delete == ['C5']


def test_number_scheme_respects_lowercase_atom_names():
    reslist = ['SFAC C H O\n', 'UNIT 1 1 1\n',
               'o1a   3    0.100000    0.200000    0.300000    11.00000    0.04\n',
               'c1a   1    0.200000    0.200000    0.300000    11.00000    0.04\n',
               'HKLF 4\n']
    dbatoms = [['O1', 3, 0.0, 0.0, 0.0], ['C1', 1, 1.0, 0.0, 0.0]]
    names = NumberScheme(reslist, dbatoms, SimpleNamespace(resiflag=False)).get_fragment_number_scheme()
    assert len(names) == 2
    assert not {'O1A', 'C1A'} & set(names)


def test_remove_adjacent_hydrogens_isotropic_atom_with_afix():
    reslist = ['SFAC C H\n', 'UNIT 1 1\n',
               'C1    1    0.100000    0.200000    0.300000    11.00000    0.04\n',
               'AFIX 43\n',
               'H1    2    0.100000    0.250000    0.300000    11.00000   -1.20000\n',
               'AFIX 0\n',
               'C2    1    0.200000    0.200000    0.300000    11.00000    0.04\n',
               'HKLF 4\n']
    fa = FindAtoms(reslist)
    assert fa.remove_adjacent_hydrogens(['C1'], ['C', 'H']) == ['H1']
    assert reslist[3:6] == ['', '', '']
    assert reslist[2].startswith('C1') and reslist[6].startswith('C2')


def test_remove_adjacent_hydrogens_anisotropic_atom_without_afix():
    reslist = ['SFAC C H\n', 'UNIT 1 1\n',
               'C1    1    0.100000    0.200000    0.300000    11.00000    0.01    0.02 =\n',
               '      0.03    0.00    0.00    0.00\n',
               'H1    2    0.100000    0.250000    0.300000    11.00000   -1.20000\n',
               'C2    1    0.200000    0.200000    0.300000    11.00000    0.04\n',
               'HKLF 4\n']
    fa = FindAtoms(reslist)
    assert fa.remove_adjacent_hydrogens(['C1'], ['C', 'H']) == ['H1']
    assert reslist[3] == '      0.03    0.00    0.00    0.00\n'
    assert reslist[4] == ''
    assert reslist[5].startswith('C2')


def test_sfac_table_without_unit_exits_with_proper_message(capsys):
    reslist = ['SFAC C H\n', 'C1    1    0.100000    0.200000    0.300000    11.00000    0.04\n', 'HKLF 4\n']
    with pytest.raises(SystemExit):
        SfacTable(reslist, ['C']).set_sfac_table()
    assert 'No UNIT instruction' in capsys.readouterr().out


class TestFindAtoms(unittest.TestCase):

    def test_get_resi_definition_dict_sorted_keys(self):
        result = FindAtoms.get_resi_definition_dict('RESI 1 TOL')
        sorted_keys = sorted(list(result.keys()))
        self.assertEqual(sorted_keys, ['class', 'number'])

    def test_get_resi_definition_dict_sorted_values(self):
        result = FindAtoms.get_resi_definition_dict('RESI 1 TOL')
        sorted_values = sorted(result.values())
        self.assertEqual(sorted_values, ['1', 'TOL'])

    def test_get_resi_definition_dict(self):
        result = FindAtoms.get_resi_definition_dict('RESI 1 TOL')
        self.assertEqual(result, {'class': 'TOL', 'number': '1'})

    def test_get_resi_definition_dict_empty(self):
        result = FindAtoms.get_resi_definition_dict('RESI')
        self.assertEqual(result, {'class': None, 'number': None})

    def test_floating_ponit_alias(self):
        with self.assertRaises(SystemExit):
            result = FindAtoms.get_resi_definition_dict('RESI 5 TOL 21.0')
            self.assertEqual(['5', 'TOL', '21.0'], result)


if __name__ == '__main__':
    unittest.main()
