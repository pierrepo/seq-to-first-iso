
"""Test files provided."""

import filecmp
from pathlib import Path

from pyteomics import mass
import pytest

import seq_to_first_iso as stfi

data_dir = Path(__file__).parent.resolve().joinpath("test_data")
REL = 1e-6

# Tests for pyteomics.
@pytest.fixture(scope="session")
def get_composition():
    return [("ACDE", {'H': 24, 'C': 15, 'O': 9, 'N': 4, 'S': 1}),
            ("VPKER", {'H': 49, 'C': 27, 'O': 8, 'N': 9}),
            ("ACDEFGHIKLMNPQRSTVWY", {'H': 159, 'C': 107, 'O': 30, 'N': 29, 'S': 2}),
            ("", {}),
            ("ACTG", {'H': 22, 'C': 12, 'O': 6, 'N': 4, 'S': 1}),
            ("LILI", {'H': 46, 'C': 24, 'O': 5, 'N': 4}),
            ("LLLL", {'H': 46, 'C': 24, 'O': 5, 'N': 4})
            ]


def test_composition(get_composition):
    assert mass.Composition("ACDE")
    assert mass.Composition("A") + mass.Composition("C")
    assert mass.Composition(parsed_sequence="ACDE") == {'H': 22, 'C': 15, 'O': 8, 'N': 4, 'S': 1}

    for data in get_composition:
        sequence = data[0]
        expected = data[1]
        assert mass.Composition(sequence) == expected


@pytest.fixture(scope="session")
def get_mass():
    return [("", 0),
            ("ACDE", 436.12639936491996),
            ("VPKER", 627.37040957111),
            ("ACDEFGHIKLMNPQRSTVWY", 2394.1249068251295),
            ("ACTG", 350.1260054421),
            ("LILI", 470.34682059221996),
            ("LLLL", 470.34682059221996),
            ("AAAAAALAAAAAIAAAAAAAAAAAAAGAAAAAAAAAAVIYWWSSEPED", 4108.065244683701),
            ]


def test_calculate_mass(get_mass):

    assert mass.calculate_mass("ACDE") == pytest.approx(436.12639936, REL)
    assert mass.calculate_mass(mass.Composition("ACDE")) == pytest.approx(436.12639936, REL)
    assert mass.calculate_mass(parsed_sequence="ACDE") == pytest.approx(418.115834, REL)
    assert mass.calculate_mass("A") == pytest.approx(89.04767846841, REL)

    for data in get_mass:
        sequence = data[0]
        expected = data[1]
        assert mass.calculate_mass(sequence) == pytest.approx(expected, REL)


# def test_fast_mass(get_mass):
#    assert mass.fast_mass("A")
#    assert mass.fast_mass2("A")
#
#    # fast_mass has a problem with empty sequences.
#    for data in get_mass[1:]:
#        sequence = data[0]
#        expected = data[1]
#        assert mass.fast_mass(sequence) == pytest.approx(expected, REL)
#        assert mass.fast_mass2(sequence) == pytest.approx(expected, REL)


# Tests for seq_to_first_iso.
def test_constants():
    assert stfi.AMINO_ACIDS == set("ACDEFGHIKLMNPQRSTVWY")
    assert stfi.isotopic_abundance == {"H[1]": 0.999885, "H[2]": 0.000115,
                                       "C[12]": 0.9893, "C[13]": 0.0107,
                                       "X[12]": 0.9893, "X[13]": 0.0107,
                                       "N[14]": 0.99632, "N[15]": 0.00368,
                                       "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                                       "S[32]": 0.9493, "S[33]": 0.0076, "S[34]": 0.0429}

    assert stfi.C12_abundance["C[12]"] == pytest.approx(0.9999, rel=REL)
    assert stfi.C12_abundance["C[13]"] == pytest.approx(0.0001, rel=REL)
    assert {k: v for k, v in stfi.C12_abundance.items()
            if k not in ["C[12]", "C[13]"]} == {"H[1]": 0.999885, "H[2]": 0.000115,
                                                "X[12]": 0.9893, "X[13]": 0.0107,
                                                "N[14]": 0.99632, "N[15]": 0.00368,
                                                "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                                                "S[32]": 0.9493, "S[33]": 0.0076, "S[34]": 0.0429}


def test_parser():
    test_file = data_dir.joinpath("sample_sequence.txt")
    bad_file = data_dir.joinpath("bad_sample_sequence.zip")
    expected_sequences = ["VPKER", "LLIDRI", "FHNK", "NEAT", "SACFTK", "NA"]

    with pytest.raises(FileNotFoundError):
        stfi.sequence_parser("not_a_file")
    with pytest.raises(UnicodeDecodeError):
        stfi.sequence_parser(bad_file)

    assert stfi.sequence_parser(test_file)

    parser_output = stfi.sequence_parser(test_file)
    assert type(parser_output) is tuple
    assert type(parser_output[0]) is list
    assert type(parser_output[1]) is int
    assert parser_output[0] == expected_sequences
    assert parser_output[1] == 4


def test_separation():
    test_sequence = "ABCDEFG"
    test_sequence_long = "ABCDEFGCDCDDC"
    test_unlabelled = "CDW"

    assert stfi.separate_labelled(test_sequence, test_unlabelled)
    assert stfi.separate_labelled(test_sequence_long, test_unlabelled)
    assert stfi.separate_labelled("", "")
    test_small = stfi.separate_labelled(test_sequence, test_unlabelled)
    test_long = stfi.separate_labelled(test_sequence_long, test_unlabelled)
    test_empty = stfi.separate_labelled("", "")

    assert type(test_small) is tuple
    assert type(test_long) is tuple
    assert type(test_empty) is tuple

    assert test_small == ("ABEFG", "CD")
    assert test_long == ("ABEFG", "CDCDCDDC")
    assert test_empty == ("", "")


def test_computation_isotopologue():
    test_composition = mass.Composition("ACDE")
    assert stfi.compute_M0_nl(test_composition, stfi.isotopic_abundance) == pytest.approx(0.77662382, REL)
    assert stfi.compute_M0_nl(test_composition, stfi.C12_abundance) == pytest.approx(0.911253268, REL)

    assert stfi.compute_M1_nl(test_composition, stfi.isotopic_abundance) == pytest.approx(0.1484942353, REL)
    assert stfi.compute_M1_nl(test_composition, stfi.C12_abundance) == pytest.approx(0.0277650369575, REL)


def test_formula_X():
    assert stfi.seq_to_midas("ACDE", "") == {'H': 24, 'C': 15, 'O': 9, 'N': 4, 'S': 1}
    assert stfi.seq_to_midas("ACDE", "FGH") == {'H': 43, 'C': 15, 'O': 12, 'N': 9, 'S': 1, 'X': 17}


def test_string_casting():
    assert stfi.formula_to_str(mass.Composition("ACDE")) == "C15H24O9N4S1"

    test_formula_X = stfi.seq_to_midas("ACDE", "FGH")
    assert stfi.formula_to_str(test_formula_X) == "C15H43O12N9S1X17"


def test_seq_to_tsv():
    sequences_given = ["VPKER", "LLIDRI", "FHNK", "NEAT", "SACFTK", "NA"]
    output_file = data_dir.joinpath("output.tsv")
    unlabelled_output_file = data_dir.joinpath("unlabelled_output.tsv")
    assert stfi.seq_to_tsv(sequences_given, output_file)
    assert output_file.is_file()
    assert stfi.seq_to_tsv(sequences_given, unlabelled_output_file, "AT")
    assert unlabelled_output_file.is_file()

    assert data_dir.joinpath("reference_sequence.tsv")
    assert filecmp.cmp(output_file, data_dir.joinpath("reference_sequence.tsv"), shallow=False)
    assert data_dir.joinpath("reference_sequence_AT.tsv")
    assert filecmp.cmp(unlabelled_output_file, data_dir.joinpath("reference_sequence_AT.tsv"), shallow=False)
