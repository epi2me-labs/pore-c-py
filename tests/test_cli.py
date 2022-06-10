import pytest
from typer.testing import CliRunner

from pore_c2.cli import app

runner = CliRunner()


@pytest.mark.parametrize("command", ["digest", "map", "merge", "haplotag", "to-cooler"])
def test_help(command: str):
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 0
