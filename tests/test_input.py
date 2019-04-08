
"""Test files provided."""

import seq_to_first_iso as stfi


def test_constants():
    assert stfi.AMINO_ACIDS == set("ACDEFGHIKLMNPQRSTVWY")
    assert stfi.isotopic_abundance ==  {"H[1]": 0.999885,
                                        "H[2]": 0.000115,
                                        "C[12]": 0.9893,
                                        "C[13]": 0.0107,
                                        "X[12]": 0.9893,
                                        "X[13]": 0.0107,
                                        "N[14]": 0.99632,
                                        "N[15]": 0.00368,
                                        "O[16]": 0.99757,
                                        "O[17]": 0.00038,
                                        "O[18]": 0.00205,
                                        "S[32]": 0.9493,
                                        "S[33]": 0.0076,
                                        "S[34]": 0.0429}


