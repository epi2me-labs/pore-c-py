import pytest
from typer.testing import CliRunner

from pore_c2.cli import app

runner = CliRunner()


@pytest.mark.parametrize("command", ["index", "map", "merge"])
def test_help(command: str):
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 0
