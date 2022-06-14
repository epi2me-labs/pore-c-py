import pytest
from typer.testing import CliRunner

from pore_c2.cli import app


@pytest.mark.parametrize("command", ["index", "map", "merge"])
def test_help(runner: CliRunner, command: str):
    result = runner.invoke(app, [command, "--help"])
    print(result.stdout)
    assert result.exit_code == 0
