import pytest


def pytest_addoption(parser):
    parser.addoption("--uchime-ref-db-fp", action="store", required=True, help="path to PR2 database")


