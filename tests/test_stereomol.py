"""stereomol tests."""

import stereomol


def test_stub() -> None:
    """Stub test to ensure the test suite runs."""
    print(stereomol.__version__)  # noqa: T201


def test__greet() -> None:
    """Test the greet function."""
    assert stereomol.greet("World") == "Hello, World!"


def test__greet_jim() -> None:
    """Test the greet_jim function."""
    assert stereomol.greet_jim() == "Hello, Jim!"
