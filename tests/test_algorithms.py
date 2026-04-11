"""Test algorithm classes with toy model."""

import pytest
import sys

sys.path.insert(0, ".")
from src.algorithms import (
    GimmeAlgorithm,
    TinitAlgorithm,
    CordaAlgorithm,
    RiptideAlgorithm,
    FastcoreAlgorithm,
    ImatAlgorithm,
)

MODEL_PATH = "tests/toy_model.xml"
EXPR_PATH = "tests/toy_expr.csv"
CONFIG_PATH = "tests/config.yaml"


class Args:
    def __init__(self, model=MODEL_PATH, expr=EXPR_PATH):
        self.model = model
        self.expr = expr


def run_algo(algo_class):
    args = Args()
    algo = algo_class(args=args, config_path=CONFIG_PATH)
    return algo.run()


def test_gimme():
    result = run_algo(GimmeAlgorithm)
    assert result is not None


@pytest.mark.xfail(reason="tINIT sets solver tolerances not supported by GLPK")
def test_tinit():
    result = run_algo(TinitAlgorithm)
    assert result is not None


@pytest.mark.xfail(reason="CORDA uses hardcoded solver")
def test_corda():
    result = run_algo(CordaAlgorithm)
    assert result is not None


def test_riptide():
    result = run_algo(RiptideAlgorithm)
    assert result is not None


def test_fastcore():
    result = run_algo(FastcoreAlgorithm)
    assert result is not None


def test_imat():
    result = run_algo(ImatAlgorithm)
    assert result is not None
