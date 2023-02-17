import mdss.parser as p
import pytest


def test_parser_with_all_arguments():
    OUT_DIR = "testing"
    PROPERTY = "RMSD"
    SELECTION = "name CA"
    SAMPLER = "RandomSampler"
    SIZE = "100"
    DISSIMILARITY = "Dissimilarity"
    TRJ_FILENAME = "input.gro"
    TOP_FILENAME = "input.xtc"
    OUT_PREFIX = "UNITTEST"

    args_list = [
        "--traj",
        TRJ_FILENAME,
        "--top",
        TOP_FILENAME,
        "--prefix",
        OUT_PREFIX,
        "--output-folder",
        OUT_DIR,
        "--property",
        PROPERTY,
        "--atom-selection",
        SELECTION,
        "--sampler",
        SAMPLER,
        "--seed-number",
        "1999",
        "--size",
        str(SIZE),
        "--dissimilarity",
        DISSIMILARITY,
    ]

    args = p.parse_args(args_list)


def test_parser_with_all_arguments():
    OUT_DIR = "testing"
    PROPERTY = "RMSD"
    SELECTION = "name CA"
    SAMPLER = "RandomSampler"
    SIZE = "100"
    DISSIMILARITY = "Dissimilarity"
    TRJ_FILENAME = "input.gro"
    TOP_FILENAME = "input.xtc"
    OUT_PREFIX = "UNITTEST"

    args_list = [
        "--traj",
        TRJ_FILENAME,
        "--top",
        TOP_FILENAME,
        "--prefix",
        OUT_PREFIX,
        "--output-folder",
        OUT_DIR,
        "--property",
        PROPERTY,
        "--atom-selection",
        SELECTION,
        "--sampler",
        SAMPLER,
        "--seed-number",
        "1999",
        "--size",
        str(SIZE),
        "--dissimilarity",
        DISSIMILARITY,
    ]

    args = p.parse_args(args_list)


def test_parser_with_missing_traj_file():
    OUT_DIR = "testing"
    PROPERTY = "RMSD"
    SELECTION = "name CA"
    SAMPLER = "RandomSampler"
    SIZE = "100"
    DISSIMILARITY = "Dissimilarity"
    TOP_FILENAME = "input.xtc"
    OUT_PREFIX = "UNITTEST"

    args_list = [
        "--traj",
        "--top",
        TOP_FILENAME,
        "--prefix",
        OUT_PREFIX,
        "--output-folder",
        OUT_DIR,
        "--property",
        PROPERTY,
        "--atom-selection",
        SELECTION,
        "--sampler",
        SAMPLER,
        "--seed-number",
        "1999",
        "--size",
        str(SIZE),
        "--dissimilarity",
        DISSIMILARITY,
    ]
    with pytest.raises(SystemExit):
        p.parse_args(args_list)


def test_parser_with_missing_seeds_number():
    OUT_DIR = "testing"
    PROPERTY = "RMSD"
    SELECTION = "name CA"
    SAMPLER = "RandomSampler"
    SIZE = "100"
    DISSIMILARITY = "Dissimilarity"
    TRJ_FILENAME = "input.gro"
    TOP_FILENAME = "input.xtc"
    OUT_PREFIX = "UNITTEST"

    args_list = [
        "--traj",
        TRJ_FILENAME,
        "--top",
        TOP_FILENAME,
        "--prefix",
        OUT_PREFIX,
        "--output-folder",
        OUT_DIR,
        "--property",
        PROPERTY,
        "--atom-selection",
        SELECTION,
        "--sampler",
        SAMPLER,
        "--seed-number",
        "--size",
        str(SIZE),
        "--dissimilarity",
        DISSIMILARITY,
    ]
    with pytest.raises(SystemExit):
        p.parse_args(args_list)
