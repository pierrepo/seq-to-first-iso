
"""Test files provided."""

import filecmp
from pathlib import Path

import pandas as pd
from pyteomics import mass
from pyteomics.auxiliary.structures import PyteomicsError
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


@pytest.fixture(scope="session")
def get_separated_sequence():
    """For C,D and W unlabelled."""
    return [("ABCDEFG", ("ABEFG", "CD")),
            ("ABCDEFGCDCDDC", ("ABEFG", "CDCDCDDC")),
            ("", ("", "")),
            ]


def test_separation(get_separated_sequence):
    test_unlabelled = "CDW"
    for data in get_separated_sequence:
        sequence = data[0]
        expected = data[1]
        assert stfi.separate_labelled(sequence, test_unlabelled)
        result = stfi.separate_labelled(sequence, test_unlabelled)
        assert type(result) is tuple
        assert result == expected


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
#    assert type(parser_output) is tuple
#    assert type(parser_output[0]) is list
#    assert type(parser_output[1]) is list
#    assert type(parser_output[2]) is int
#    assert parser_output[1] == expected_sequences
#    assert parser_output[2] == 4
    assert type(parser_output) is dict
    assert type(parser_output.get("annotations")) is list
    assert type(parser_output.get("sequences")) is list
    assert type(parser_output.get("ignored_lines")) is int
    assert parser_output.get("sequences") == expected_sequences
    assert parser_output.get("ignored_lines") == 4


def test_parser_annotation(caplog):
    empty_file = data_dir.joinpath("sample_empty_file.tsv")
    test_file = data_dir.joinpath("sample_sequence.tsv")
    expected_sequences = ["VPKER", "LLIDRI", "FHNK", "NEAT", "SACFTK", "NA"]

    empty_output = stfi.sequence_parser(empty_file, sep="")
    assert "file is empty" in caplog.text
    assert "separator is empty" in caplog.text
    # No annotations or sequences were taken.
    assert empty_output.get("annotations") == empty_output.get("sequences") == []

    parser_output = stfi.sequence_parser(test_file)
    assert type(parser_output) is dict
    assert type(parser_output.get("annotations")) is list
    assert type(parser_output.get("sequences")) is list
    assert type(parser_output.get("ignored_lines")) is int
    assert parser_output.get("sequences") == expected_sequences
    assert parser_output.get("ignored_lines") == 3


def test_deprecated_computation_isotopologue():
    test_composition = mass.Composition("ACDE")
    stfi.compute_M0 = stfi.seq_to_first_iso.compute_M0
    stfi.compute_M1 = stfi.seq_to_first_iso.compute_M1
    assert stfi.compute_M0(test_composition, stfi.isotopic_abundance) == pytest.approx(0.77662382, REL)
    assert stfi.compute_M0(test_composition, stfi.C12_abundance) == pytest.approx(0.911253268, REL)

    assert stfi.compute_M1(test_composition, stfi.isotopic_abundance) == pytest.approx(0.1484942353, REL)
    assert stfi.compute_M1(test_composition, stfi.C12_abundance) == pytest.approx(0.0277650369575, REL)


def test_formula_X():
    assert stfi.seq_to_midas("ACDE", "") == {'H': 24, 'C': 15, 'O': 9, 'N': 4, 'S': 1}
    assert stfi.seq_to_midas("ACDE", "FGH") == {'H': 43, 'C': 15, 'O': 12, 'N': 9, 'S': 1, 'X': 17}


def test_computation_isotopologue():
    # Standard formula.
    test_composition = mass.Composition("ACDE")
    assert stfi.compute_M0_nl(test_composition, stfi.isotopic_abundance) == pytest.approx(0.77662382, REL)
    assert stfi.compute_M0_nl(test_composition, stfi.C12_abundance) == pytest.approx(0.911253268, REL)

    assert stfi.compute_M1_nl(test_composition, stfi.isotopic_abundance) == pytest.approx(0.1484942353, REL)
    assert stfi.compute_M1_nl(test_composition, stfi.C12_abundance) == pytest.approx(0.0277650369575, REL)

    unlabelled_composition = stfi.seq_to_midas("AC", "DE")
    assert stfi.compute_M0_nl(unlabelled_composition, stfi.isotopic_abundance) == pytest.approx(0.77662382, REL)
    assert stfi.compute_M0_nl(unlabelled_composition, stfi.C12_abundance) == pytest.approx(0.8279079739944033, REL)

    assert stfi.compute_M1_nl(unlabelled_composition, stfi.isotopic_abundance) == pytest.approx(0.1484942353, REL)
    assert stfi.compute_M1_nl(unlabelled_composition, stfi.C12_abundance) == pytest.approx(0.10507024116588572, REL)


def test_string_casting():
    assert stfi.formula_to_str(mass.Composition("ACDE")) == "C15H24O9N4S1"

    test_formula_X = stfi.seq_to_midas("ACDE", "FGH")
    assert stfi.formula_to_str(test_formula_X) == "C15H43O12N9S1X17"


@pytest.fixture(scope="session")
def get_mods():
    return [(["Oxidation"], mass.Composition({"O":1})),
            (["Acetyl", "Phospho"], mass.Composition({'H':3, 'C':2, 'O':4, "P":1})),
            (["Acetyl", "Phospho", "not_mod"], mass.Composition({'H':3, 'C':2, 'O':4, "P":1})),
            ([], mass.Composition()),
            ]


def test_get_mods_composition(get_mods, caplog):
    get_mods_composition = stfi.seq_to_first_iso.get_mods_composition
    get_mods_composition(["not_mod"])
    assert "entry not found" in caplog.text
    # Mod with element not in CHONPSX.
    get_mods_composition(["Heme"])
    assert "Fe is not supported" in caplog.text
    for data in get_mods:
        modifications = data[0]
        expected = data[1]
        assert get_mods_composition(modifications) == expected


def test_seq_to_tsv(caplog):
#    sequences_given = ["VPKER", "LLIDRI", "FHNK", "NEAT", "SACFTK", "NA"]
#    output_file = data_dir.joinpath("output.tsv")
#    unlabelled_output_file = data_dir.joinpath("unlabelled_output.tsv")
    # TODO: add dataframe comparison.
#    assert stfi.seq_to_tsv(sequences_given, unlabelled_aa=[])
#    assert output_file.is_file()
#    assert stfi.seq_to_tsv(sequences_given, unlabelled_output_file, unlabelled_aa="AT")
#    assert unlabelled_output_file.is_file()
#
#    assert data_dir.joinpath("reference_sequence.tsv")
#    assert filecmp.cmp(output_file, data_dir.joinpath("reference_sequence.tsv"), shallow=False)
#    assert data_dir.joinpath("reference_sequence_AT.tsv")
#    assert filecmp.cmp(unlabelled_output_file, data_dir.joinpath("reference_sequence_AT.tsv"), shallow=False)
    # Parse non valid sequence.
    with pytest.raises(PyteomicsError):
        stfi.seq_to_tsv(["b"], [])
    # Testing  annotations and sequences with different lengths.
    df = stfi.seq_to_tsv(["AC"], [], annotations=["id1", "id2"])
    assert "different lengths" in caplog.text
    assert type(df) is pd.DataFrame
    # Annotations and sequences with same length.
    df = stfi.seq_to_tsv(["AC"], [], annotations=["id1"])
    assert type(df) is pd.DataFrame

    # Test verification with modifications
    no_mod_input = {"sequences": ["A"], "unlabelled_aa":[],
                    "raw_sequences":["A"], "not_modifications":[]}
    no_rseq_input = {"sequences": ["A"], "unlabelled_aa":[],
                     "modifications":["Oxidation"]}
    stfi.seq_to_tsv(**no_mod_input)
    assert "not_modifications not recognized" in caplog.text
    assert "raw_sequences and modifications have different" in caplog.text
    stfi.seq_to_tsv(**no_rseq_input)
    assert "raw_sequences and sequences have different" in caplog.text


def test_cli_parser(caplog):
    cli_parser = stfi.seq_to_first_iso.user_input
    test_file = data_dir.joinpath("sample_sequence.txt")

    assert cli_parser([str(test_file)])
    with pytest.raises(SystemExit):
        cli_parser(["not_a_file"])

    minimal_cli = cli_parser([str(test_file)])
    assert minimal_cli.non_labelled_aa == []
    # It seems like a space between -n and A here won't be escaped.
    cli_labelled = cli_parser([str(test_file), "-nA,c"])
    assert cli_labelled.non_labelled_aa == ["A", "C"]
    cli_parser([str(test_file), "-nerror"])
    assert "WARNING" in caplog.text


def test_main(caplog):
    main_cli = stfi.seq_to_first_iso.cli
    test_file = data_dir.joinpath("sample_sequence.txt")
    bad_file = data_dir.joinpath("sample_bad_sequence.txt")

    main_cli([str(test_file), "-nA"])
    assert "Amino acid" in caplog.text
    assert "lines ignored" in caplog.text
    assert Path("sample_sequence_stfi.tsv").is_file()

    main_cli([str(test_file), "-ooutput"])
    assert Path("output.tsv").is_file()

    with pytest.raises(SystemExit):
        main_cli([str(bad_file)])
    with pytest.raises(SystemExit):
        # Error raised due to not having command line arguments.
        main_cli([])

