import pytest
from typer.testing import CliRunner


@pytest.fixture
def runner():
    runner = CliRunner()
    return runner
