import unittest
from io import StringIO
from unittest.mock import patch

from atomhandling import FindAtoms


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
