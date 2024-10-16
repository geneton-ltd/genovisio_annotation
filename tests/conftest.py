import pytest

import annotation


@pytest.fixture(scope="session")
def annotation_cnv():
    yield annotation.Annotation.load_from_json("tests/data/cnv.json").cnv
