import annotation


def test_gnomad():
    annot = annotation.Annotation.load_from_json("tests/data/gnomad.json")

    var_regions = annot.get_common_variability_regions()
    assert len(var_regions) == 1
    assert var_regions[0]["frequency"] == 0.2
